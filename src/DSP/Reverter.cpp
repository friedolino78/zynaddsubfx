/*
  ZynAddSubFX - a software synthesizer

  Reverter.cpp - Reverse Delay
  Copyright (C) 2023-2024 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cstring>

#include "../Misc/Allocator.h"
#include "../Misc/Util.h"
#include "Reverter.h"

namespace zyn {

Reverter::Reverter(Allocator *alloc, float delay_, unsigned int srate, int bufsize, float tRef_, AbsTime *time_)
    : syncMode(AUTO), input(nullptr), gain(1.0f), delay(delay_), phase(0.0f), crossfade(0.16f),
      tRef(tRef_), buffer_offset(0), buffer_counter(0), reverse_index(0.0f), phase_offset_old(0.0f),
      phase_offset_fade(0.0f), fade_counter(0), rms_hist(999.9f), time(time_), memory(*alloc) 
{
    samplerate = srate;
    buffersize = bufsize;
    fading_samples = static_cast<int>(srate * crossfade);
    max_delay = srate * MAX_REV_DELAY_SECONDS;
    mem_size = static_cast<int>(ceilf(max_delay * 4.0f)) + static_cast<int>(1.27f * samplerate) + 2;

    input = static_cast<float *>(memory.alloc_mem(mem_size * sizeof(float)));
    reset();

    setdelay(delay);
    global_offset = fmodf(tRef, delay);
    pos_writer = 0;
    pos_reader = 0;
    pos_start = 0;
}

Reverter::~Reverter() {
    memory.dealloc(input);
}

inline float Reverter::sampleLerp(float *smp, float pos) {
    int poshi = static_cast<int>(pos);
    float poslo = pos - static_cast<float>(poshi);
    return smp[poshi] + poslo * (smp[poshi + 1] - smp[poshi]);
}

inline float hanningWindow(float x) {
    return 0.5f * (1.0f - cos(1.0f * M_PI * x));
}

inline void Reverter::switchBuffers(float offset) {
    reverse_index = 0;
    pos_start = pos_writer + buffersize + offset;
    float pos_next = fmodf(float(pos_start + mem_size) - (reverse_index + phase_offset), mem_size);
    delta_crossfade = pos_reader - 1.0f - pos_next;
    fade_counter = 0;
}

void Reverter::filterout(float *smp) {
    writeToRingBuffer(smp);
    processBuffer(smp);
}

void Reverter::writeToRingBuffer(float *smp) {
    int space_to_end = mem_size - pos_writer;
    if (buffersize <= space_to_end) {
        // No wrap around, copy in one go
        memcpy(&input[pos_writer], smp, buffersize * sizeof(float));
    } else {
        // Wrap around, split into two copies
        memcpy(&input[pos_writer], smp, space_to_end * sizeof(float));
        memcpy(&input[0], smp + space_to_end, (buffersize - space_to_end) * sizeof(float));
    }

// Update pos_writer
pos_writer = (pos_writer + buffersize) % mem_size;
}

void Reverter::processBuffer(float *smp) {
    float rms = 0.0f;
    for (int i = 0; i < buffersize; i++) {
        reverse_index++;
        rms += fabsf(smp[i]);
        phase_offset = phase_offset_old + static_cast<float>(i) * phase_offset_fade;

        handleSync(i);
        updateReaderPosition();
        crossfadeSamples(smp, i);
        applyGain(smp[i]);
    }

    rms_hist = rms / static_cast<float>(buffersize);
    phase_offset_old = phase_offset;
    phase_offset_fade = 0.0f;
}

void Reverter::handleSync(int i) {
    switch (syncMode) {
        case AUTO:
            if (reverse_index >= delay && state != IDLE) switchBuffers(delay / 2.0f);
            break;
        case HOST:
        case MIDI:
            if (doSync && reverse_index >= syncPos) switchBuffers(delay / 2.0f);
            break;
        case NOTEON:
        case NOTEONOFF:
            handleNoteSync(i);
            break;
    }
    if (reverse_index >= max_delay && state == PLAYING) switchBuffers(delay / 2.0f);
}

void Reverter::handleNoteSync(int i) {
    if (syncMode == NOTEON && reverse_index >= delay && state != IDLE) {
        handleStateChange("NOTEON", i);
    } else if (syncMode == NOTEONOFF && (
                   reverse_index >= recorded_samples && state == PLAYING ||
                   reverse_index >= max_delay && state == PLAYING ||
                   rms_hist < 0.001f && state == RECORDING)) {
        handleStateChange("NOTEONOFF", i);
    }
}

void Reverter::handleStateChange(const char *mode, int i) {
    if (state == RECORDING) {
        recorded_samples = reverse_index;
        state = PLAYING;
    } else if (state == PLAYING) {
        state = IDLE;
    }
    switchBuffers(delay / 2.0f);
    printf("syncMode: %s state change at index %d\n", mode, i);
}

void Reverter::updateReaderPosition() {
    pos_reader = fmodf(float(pos_start + mem_size) - (reverse_index + phase_offset), mem_size);
}

void Reverter::crossfadeSamples(float *smp, int i) {
    if (fade_counter < fading_samples) {
        float slope = static_cast<float>(fade_counter) / static_cast<float>(fading_samples);
        float fadein = hanningWindow(slope);
        float fadeout = 1.0f - fadein;
        fade_counter++;

        if (state != IDLE) {
            smp[i] = applyFade(fadein, fadeout);
        } else {
            smp[i] = fadeout * sampleLerp(input, fmodf(pos_reader + mem_size + delta_crossfade, mem_size));
        }
    } else {
        smp[i] = (state == PLAYING) ? sampleLerp(input, pos_reader) : 0.0f;
    }
}

float Reverter::applyFade(float fadein, float fadeout) {
    if (syncMode == NOTEON || syncMode == NOTEONOFF) {
        return fadein * sampleLerp(input, pos_reader);
    } else {
        return fadein * sampleLerp(input, pos_reader) + fadeout * sampleLerp(input, fmodf(pos_reader + mem_size + delta_crossfade, mem_size));
    }
}

void Reverter::applyGain(float &sample) {
    sample *= gain;
}

void Reverter::sync(float pos) {
    if (state == IDLE) {
        state = RECORDING;
        reverse_index = 0.0f;
        rms_hist = 1.0f;
    } else {
        syncPos = pos + reverse_index;
        doSync = true;
    }
}

void Reverter::update_phase(float value) {
    float phase_offset_new = value * delay;
    phase_offset_fade = (phase_offset_new - phase_offset_old) / static_cast<float>(buffersize);
}

void Reverter::setdelay(float value) {
    delay = value * static_cast<float>(samplerate);
    fading_samples = static_cast<int>(crossfade * static_cast<float>(samplerate));
    if (delay < 2.0f * static_cast<float>(fading_samples)) fading_samples = static_cast<int>(delay * 0.5f);

    global_offset = fmodf(tRef, delay);
    update_phase(phase);
}

void Reverter::setphase(float value) {
    phase = value;
    update_phase(phase);
}

void Reverter::setcrossfade(float value) {
    crossfade = value;
    fading_samples = static_cast<int>(crossfade * static_cast<float>(samplerate));
    if (delay < 2.0f * static_cast<float>(fading_samples)) fading_samples = static_cast<int>(delay * 0.5f);
}

void Reverter::setgain(float value) {
    gain = dB2rap(value);
}

void Reverter::setsyncMode(SyncMode value) {
    if (value != syncMode) {
        syncMode = value;
        state = (syncMode == NOTEON || syncMode == NOTEONOFF) ? IDLE : PLAYING;
        printf("setting syncMode: %d state: %d\n", syncMode, state);
    }
}

void Reverter::reset() {
    memset(input, 0, mem_size * sizeof(float));
}

};
