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
#include <vector>

namespace zyn {



class MoogFilter:public Filter
{
    public:
        //! @param Fq resonance, range [0.1,1000], logscale
        MoogFilter(float Ffreq, float Fq,
                unsigned char non_linear_element /* currently set by "stages" */,
                unsigned int srate, int bufsize);
        ~MoogFilter() override;
        void filterout(float *smp) override;
        void setfreq(float frequency/*frequency*/) override;
        void setfreq_and_q(float frequency, float q_) override;
        void setq(float q_/*q_*/) override;
        void setgain(float dBgain) override;
    private:
        unsigned sr;
        float gain;

        float step(float x);
        std::vector<float> impulse_response(float alpha, float k);
        float tanhd(const float x, const float d, const float s);

        float b[4] = { 0, 0, 0, 0 };
        float compensation, estimate, c2, fb;
        float xx, y0, y1, y2, y3;
        float t0, t1, t2, t3;
        float g0, g1, g2, g3;
        float z0, z1, z2, z3;
        float f0, f1, f2, f3;
        float cgfbr, f, fd2;
        float a = 2.0f;
        float s = 0.1f;
        float d = 1.0f;

        float c, q;
        //float par1, par2;

};

}
