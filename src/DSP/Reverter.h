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

#define MAX_CROSSFADE_SECONDS 1.27f

namespace zyn {

#define SYNCMODES   AUTO,\
                    AUTOFLIP,\
                    MIDI,\
                    HOST,\
                    HOSTFLIP,\
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

class Reverter
{
    public:

        SyncMode syncMode;

        //! @param Fq resonance, range [0.1,1000], logscale
        Reverter(Allocator *alloc, float delay,
                unsigned int srate, int bufsize, float tRef=0.0f, AbsTime *time_=nullptr);
        ~Reverter();
        void filterout(float *smp);
        void setdelay(float value);
        void setphase(float value);
        void setcrossfade(float value);
        void setgain(float value);
        void setsyncMode(SyncMode value);
        void reset();
        void sync(float pos);

        State state;

    private:

        void update_phase(float phase);
        void switchBuffers();
        void nextState();
        void flipState();

        float* input;
        float gain;

        float step(float x);
        float tanhX(const float x);
        float sampleLerp(float *smp, float pos);

        float delay;
        float phase;
        float crossfade;

        float tRef;
        int buffer_offset;
        int buffer_counter;
        float global_offset;
        float reverse_index;
        float phase_offset_old;
        float phase_offset_fade;
        int fading_samples;
        int fade_counter;
        float mav;

        AbsTime *time;

        Allocator &memory;
        unsigned int mem_size;
        int samplerate;
        int buffersize;
        float max_delay;

        unsigned int pos_start;
        float phase_offset, recorded_samples;
        float syncPos, pos_reader, crossfade_offset;
        bool doSync;

        unsigned int pos_writer = 0;
};

}
