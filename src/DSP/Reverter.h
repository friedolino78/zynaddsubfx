/*
  ZynAddSubFX - a software synthesizer

  Reverter.h - Reverse Delay
  Copyright (C) 2023-2024 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#pragma once
#include "Filter.h"
#include "Value_Smoothing_Filter.h"

namespace zyn {

#define SYNCMODES   AUTO,\
                    MIDI,\
                    HOST,\
                    NOTEON,\
                    NOTEONOFF

enum SyncMode
{
    SYNCMODES
};

#define STATE   RECORDING,\
                PLAYING,\
                IDLE

enum State
{
    STATE
};

/**
 * @brief Reverse Delay effect class
 * 
 * Implements a reverse delay effect using a ring buffer. Supports synchronization
 * with host, MIDI, or note-on/off events.
 */
class Reverter
{
public:
    /**
     * @brief Constructor for the Reverter class.
     * 
     * @param alloc Allocator for memory management.
     * @param delay Initial delay time in seconds.
     * @param srate Sample rate in Hz.
     * @param bufsize Size of the audio buffer.
     * @param tRef Optional time reference (default is 0.0f).
     * @param time_ Optional pointer to an AbsTime object for timing.
     */
    Reverter(Allocator *alloc, float delay, unsigned int srate, int bufsize, float tRef = 0.0f, AbsTime *time_ = nullptr);

    /// Destructor
    ~Reverter();

    /**
     * @brief Process and filter the input audio samples.
     * 
     * @param smp Pointer to the buffer of input samples.
     */
    void filterout(float *smp);

    /**
     * @brief Set the delay time in seconds.
     * 
     * @param value Delay time.
     */
    void setdelay(float value);

    /**
     * @brief Set the phase of the delay.
     * 
     * @param value Phase value.
     */
    void setphase(float value);

    /**
     * @brief Set the crossfade time for the transition between buffers.
     * 
     * @param value Crossfade time.
     */
    void setcrossfade(float value);

    /**
     * @brief Set the output gain.
     * 
     * @param value Gain value in dB.
     */
    void setgain(float value);

    /**
     * @brief Set the synchronization mode for the reverter.
     * 
     * @param value Synchronization mode (e.g., AUTO, MIDI, HOST).
     */
    void setsyncMode(SyncMode value);

    /// Reset the internal state and buffer.
    void reset();

    /**
     * @brief Synchronize the delay based on an external position.
     * 
     * @param pos External position for syncing.
     */
    void sync(float pos);

private:
    /**
     * @brief Update the phase of the delay.
     * 
     * @param phase New phase value.
     */
    void update_phase(float phase);

    /**
     * @brief Switch between buffers when reverse playback is triggered.
     * 
     * @param offset Offset for buffer switching.
     */
    void switchBuffers(float offset);

    /**
     * @brief Perform linear interpolation between two samples.
     * 
     * @param smp Pointer to the sample buffer.
     * @param pos Position to interpolate.
     * @return Interpolated sample value.
     */
    float sampleLerp(float *smp, float pos);

    /**
     * @brief Write new input samples to the ring buffer.
     * 
     * @param smp Pointer to the buffer of input samples.
     */
    void writeToRingBuffer(float *smp);

    /**
     * @brief Process the samples and apply the reverse delay effect.
     * 
     * @param smp Pointer to the buffer of input samples.
     */
    void processBuffer(float *smp);

    /**
     * @brief Handle synchronization modes and state transitions.
     * 
     * @param i Index of the current sample in the buffer.
     */
    void handleSync(int i);

    /**
     * @brief Handle synchronization in NOTEON and NOTEONOFF modes.
     * 
     * @param i Index of the current sample in the buffer.
     */
    void handleNoteSync(int i);

    /**
     * @brief Handles state transitions in NOTEON and NOTEONOFF modes.
     * 
     * @param mode Name of the synchronization mode triggering the state change.
     * @param i Index of the current sample in the buffer.
     */
    void handleStateChange(const char *mode, int i);

    /**
     * @brief Update the read position in the ring buffer.
     */
    void updateReaderPosition();

    /**
     * @brief Apply crossfading between two buffer segments.
     * 
     * @param smp Pointer to the buffer of input samples.
     * @param i Index of the current sample in the buffer.
     */
    void crossfadeSamples(float *smp, int i);

    /**
     * @brief Apply fading for transitions between buffer segments.
     * 
     * @param fadein Fade-in factor.
     * @param fadeout Fade-out factor.
     * @return Crossfaded sample value.
     */
    float applyFade(float fadein, float fadeout);

    /**
     * @brief Apply the output gain to a sample.
     * 
     * @param sample Reference to the sample to modify.
     */
    void applyGain(float &sample);

    /// Current synchronization mode
    SyncMode syncMode;

    /// Current state of the reverter (RECORDING, PLAYING, IDLE)
    State state;

    /// Ring buffer for input samples
    float* input;

    /// Output gain value
    float gain;

    /// Delay time in samples
    float delay;

    /// Phase value
    float phase;

    /// Crossfade duration in samples
    float crossfade;

    /// Time reference for syncing
    float tRef;

    /// Offset in the buffer for reading/writing
    int buffer_offset;

    /// Counter for buffer processing
    int buffer_counter;

    /// Offset for global time synchronization
    float global_offset;

    /// Index for reverse playback
    float reverse_index;

    /// Phase offset for delay transitions
    float phase_offset_old;
    float phase_offset_fade;

    /// Number of samples for crossfading
    int fading_samples;

    /// Counter for crossfade processing
    int fade_counter;

    /// Root mean square (RMS) history for input volume analysis
    float rms_hist;

    /// Pointer to timing object for synchronization
    AbsTime *time;

    /// Memory allocator
    Allocator &memory;

    /// Total size of the memory buffer in samples
    unsigned int mem_size;

    /// Sample rate of the audio engine
    int samplerate;

    /// Size of the audio buffer in samples
    int buffersize;

    /// Maximum allowable delay in samples
    float max_delay;

    /// Start position for the read head in the ring buffer
    unsigned int pos_start;

    /// Offset for phase correction
    float phase_offset;

    /// Number of samples recorded for playback
    float recorded_samples;

    /// Position to synchronize playback
    float syncPos;

    /// Current position of the read head in the ring buffer
    float pos_reader;

    /// Delta between crossfaded buffer positions
    float delta_crossfade;

    /// Flag indicating whether synchronization is active
    bool doSync;

    /// Write head position in the ring buffer
    unsigned int pos_writer;
};

}
