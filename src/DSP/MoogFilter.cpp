#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "MoogFilter.h"

namespace zyn{

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )

MoogFilter::MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
    unsigned char Fstages,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f), stages(Fstages), type(Ftype)
{
    this->make_filter(Ffreq/srate, Fq);
    for (int i = 0; i>4; i++)
    {
        b[i] = 0.0;
    }
}

MoogFilter::~MoogFilter(void)
{

}

void MoogFilter::make_filter(float ff, float q)
{
    fb = cbrtf(q/1000)*4.0; // flattening
    compensation = 1.0f + LIMIT(fb, 1.0, 0);

    ff = LIMIT(ff,0.50f,0.0001f); // prevent aliasing of the resonance and its harmonics
    f = ff * 3.745;
    c = f;
    cm2 = c * 2.0;
    cp2 = c * c;
    cp3 = cp2 * c;
    cp4 = cp2 * cp2;
}

inline float MoogFilter::tanhdx(const float x, const float d = 1.0f, const float s=0.1f)
{    
    return 1.0f - s * (d + 1.0f) * x * x / (d + x * x); // Cest: 1.69
}

inline float MoogFilter::tanhd(const float a)
{    
    float x = 0.98 * a; // reduce drive of the nle to reduce aliasing
    x2 = x*x; 
    return((10 * x2*x + 105*x)/(0.98 * (x2*x2 + 45*x2 + 105))); // Pade approximation of tanh; Cest 0.88
}

float MoogFilter::step(float input)
{

// `Cheap non-linear zero-delay filters',
//// LICENSE TERMS: Copyright 2012 Teemu Voipio
// 
// You can use this however you like for pretty much any purpose,
// as long as you don't claim you wrote it. There is no warranty.
//
// Distribution of substantial portions of this code in source form
// must include this copyright notice and list of conditions.
//

// version from aciddose
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=349859&start=285 (function cascade_4)
// slightly improved for performance (precalculations + tanh approximation) by Friedolino78

    // evaluate the non-linear gains

    t0 = tanhd(b[0] + a);
    t1 = tanhd(b[1] + a);
    t2 = tanhd(b[2] + a);
    t3 = tanhd(b[3] + a);

    // denominators for solutions of individual stages
    g0 = 1.0f / (1.0f + c*t0);
    g1 = 1.0f / (1.0f + c*t1);
    g2 = 1.0f / (1.0f + c*t2);
    g3 = 1.0f / (1.0f + c*t3);

    t2g3 = t2 * g3;
    t1g2 = t1 * g2;
    t0g1 = t0 * g1;
    t1g2t2g3 = t1g2 * t2g3;

    // factored out of the feedback solution
    f3 = c * t2g3;
    f2 = cp2 * t1g2t2g3;
    f1 = cp3 * t0g1 * t1g2t2g3;
    f0 = cp4 * g0 * t0g1 * t1g2t2g3;

    // solve feedback 
    estimate =
        g3 * b[3] +
        f3 * g2 * b[2] +
        f2 * g1 * b[1] +
        f1 * g0 * b[0] +
        f0 * input;

    z0 = c*t0 / (1.0f + c*t0);
    z1 = c*t1 / (1.0f + c*t1);
    z2 = c*t2 / (1.0f + c*t2);
    z3 = c*t3 / (1.0f + c*t3);

    // then solve the remaining outputs (with the non-linear gains here)
    cgfbr = 1.0f / (1.0f + fb * z0*z1*z2*z3);
    xx = input - tanhd(fb * estimate) * cgfbr;
    y0 = t0 * g0 * (b[0] + c * xx);
    y1 = t1 * g1 * (b[1] + c * y0);
    y2 = t2 * g2 * (b[2] + c * y1);
    y3 = t3 * g3 * (b[3] + c * y2);

    // update state
    b[0] += cm2 * (xx - y0);
    b[1] += cm2 * (y0 - y1);
    b[2] += cm2 * (y1 - y2);
    b[3] += cm2 * (y2 - y3);

    // compensate the passband reduction by the negative feedback
    return y3 * compensation;
}

void MoogFilter::filterout(float *smp)
{
    for (int i = 0; i < buffersize; i ++)
        {
            // limit input amplitude to prevent overflow at high input levels + high gain
            // with res, gain and all input levels at max, it may still overflow
            smp[i] = this->step(LIMIT(smp[i]*gain, 0.85, -0.85 )); 
            smp[i] *= outgain;
        }
}

void MoogFilter::setfreq_and_q(float frequency, float q_)
{
    this->make_filter(frequency/sr, q_); // for type cheap
}

void MoogFilter::setfreq(float frequency)
{

}

void MoogFilter::setq(float q_)
{

}

void MoogFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}

void MoogFilter::settype(unsigned char type_)
{
    type = type_;
}

void MoogFilter::setstages(unsigned char stages_)
{

}

};
