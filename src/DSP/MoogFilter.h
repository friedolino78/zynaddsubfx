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


        
class MoogFilter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
                unsigned char Fstages,
                unsigned int srate, int bufsize);
        ~MoogFilter() override;
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
        
// for Cheap non-linear zero-delay filters 
        float tanhdx(const float x, const float d, const float s);
        float tanhd(const float x);
        float smoothABS( float x, const float y);
        float smoothclip(float x, const float a, const float b);
        float tanhsoftlimiter(float x, const float gain, const float offset = 0.0, const float par = 8.0);
        float softlimiter8 (float x, const float drive, const float offset, const float ammount);
        float softlimiter4 (float x, const float drive, const float offset, const float ammount)
        float x2;
        float fb, f;
        float ff;
        
        float b[4] = { 0, 0, 0, 0 };
        float compensation, estimate, c, cm2, cp2, cp3, cp4;
        float xx, y0, y1, y2, y3;
        float t0, t1, t2, t3;
        
        float t2g3;
        float t1g2;
        float t0g1;
        float t1g2t2g3;
        
        float g0, g1, g2, g3;
        float z0, z1, z2, z3;
        float f0, f1, f2, f3;
        float cgfbr, fd2;

        const float a = 2.0f;
        const float s = 0.1f;
        const float d = 1.0f;
        
//ImprovedMoog
        //~ float dV0, dV1, dV2, dV3;
        //~ float V[4];
        //~ float dV[4];
        //~ float tV[4];
        //~ float cutoff, x, g, VT, sr2;

        
};

}
