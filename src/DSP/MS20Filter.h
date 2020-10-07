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
#include "Filter.h"
//~ #include "ImprovedMoog.h"
#include <vector>

namespace zyn {


        
class MS20Filter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        MS20Filter(unsigned char Ftype, float Ffreq, float Fq,
                unsigned char Fstages,
                unsigned int srate, int bufsize);
        ~MS20Filter() override;
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
        //~ ImprovedMoog* im;
        
        float step(float x);
        void make_filter(float f, float q);
        std::vector<float> impulse_response(float alpha, float k);
        

        float tanhdx(const float x, const float d, const float s);
        float tanhd(const float x);
        float smoothABS( float x, const float y);
        float smoothclip(float x, const float a, const float b);
        float tanhsoftlimiter(float x, const float gain, const float offset = 0.0, const float par = 8.0);
        float softlimiter8 (float x, const float drive, const float offset, const float ammount);
        float softlimiter4 (float x, const float drive, const float offset, const float ammount);
        float f_dg(float s);
        float f_g(float s);
        
        // https://github.com/JoepVanlier/JSFX/blob/master/Basics/MS-20.jsfx (MIT License)
        float i, y1, y2, d1, d2, h, hh, k, obs;
        float f, fb;
        float x2;


        
};

}
