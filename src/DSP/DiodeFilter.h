/*
  ZynAddSubFX - a software synthesizer

  Moog Filter.h - Several analog filters (lowpass, highpass...)
  Copyright (C) 2018-2018 Mark McCurry
  Author: Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#pragma once
#include "DiodeLadderFilter.h"
#include "Filter.h"
//~ #include "ImprovedMoog.h"
#include <vector>

namespace zyn {


        
class DiodeFilter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        DiodeFilter(unsigned char Ftype, float Ffreq, float Fq,
                unsigned char Fstages,
                unsigned int srate, int bufsize);
        ~DiodeFilter() override;
        void filterout(float *smp) override;
        void setfreq(float /*frequency*/) override;
        void setfreq_and_q(float frequency, float q_) override;
        void setq(float /*q_*/) override;
        void setgain(float dBgain) override;
        void settype(unsigned char type);
        void setstages(unsigned char stages);
    private:
        unsigned sr;
        float gain;
        int stages;
        unsigned char type;
        DiodeLadderFilter* df;

        
        float step(float x);
        std::vector<float> impulse_response(float alpha, float k);
        

        float fb, f, fc;
        float ff = 0.1; 
        float ff_old = 0.1;
        
};

}
