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

namespace zyn {

class MoogFilter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
                unsigned int srate, int bufsize);
        ~MoogFilter() override;
        void filterout(float *smp) override;
        void setfreq(float /*frequency*/) override;
        void setfreq_and_q(float frequency, float q_) override;
        void setq(float /*q_*/) override;
        void setgain(float dBgain) override;
        void settype(unsigned char type); //

    private:
        unsigned sr;
        float gain;
        unsigned char type;
        
        float step(float x);
        
        // for "Cheap non-linear zero-delay filter" 
        float tanhXdX(const float x);
        float tanhX(const float x);
        float tan_2(const float x);
        float x2, x4;
        float fb;
        float ff;
        float b[4] = {0.0f,0.0f,0.0f,0.0f};
        float compensation, estimate, c, cm2, cp2, cp3, cp4;
        float xx, y0, y1, y2, y3;
        float t0, t1, t2, t3;
        float t1g2t2g3;
        float g0, g1, g2, g3;
        float cmt0, cmt1, cmt2, cmt3;
        float z0, z1, z2, z3;
        float f0, f1, f2, f3;
        float cgfbr, fd2;
        //~ float smp_t1, smp_t0;
};

}
