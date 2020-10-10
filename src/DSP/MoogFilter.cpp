#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "../Misc/Util.h"
#include "MoogFilter.h"

#define CLAMP(x) tanh_5_4(x) 

namespace zyn{

MoogFilter::MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f), type(Ftype)
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
    fb = cbrtf(q/1000.0)*3.9 + 0.1; // flattening
    compensation = 1.0f + limit(fb, (float)0.0, (float)1.0);

    ff = limit(ff,0.0002f,0.49f);   // limit cutoff to prevent overflow
    c = tan_2(PI * ff);             // pre warp cutoff to map to reality
    
    cm2 = c * 2.0;      // pre calculate some stuff outside the hot zone
    cp2 = c * c;
    cp3 = cp2 * c;
    cp4 = cp2 * cp2;
}

inline float MoogFilter::tan_2(const float x)
{    
    //Pade approximation tan(x) hand tuned to map fCutoff 
    x2 = x*x;
    A = 4.75*(22.3*x - 2*x2*x);
    B = 105 - 45*x2 + x2*x2;
    return (A/B);
}

inline float MoogFilter::tanh_5_4(const float x)
{   
    // Pade approximation of tanh(x) used for input saturation
    x2 = x*x; 
    x4 = x2*x2;
    return x*(945+105*x2+x4)/(945+420*x2+15*x4);
}


inline float MoogFilter::tanhXdX(float x)
{
    // Pade approximation for tanh(x)/x used in filter stages
    x2 = x*x;
    return ((x2 + 105)*x2 + 945) / ((15*x2 + 420)*x2 + 945);
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

    // evaluate the transconductance gM(vIn) = iCtl * ( tanh( vIn ) / vIn )
    t0 = tanhXdX(b[0]);
    t1 = tanhXdX(b[1]);
    t2 = tanhXdX(b[2]);
    t3 = tanhXdX(b[3]);

    // denominators for solutions of individual stages
    cmt0 = c*t0;
    cmt1 = c*t1;
    cmt2 = c*t2;
    cmt3 = c*t3;

    g0 = 1.0f / (1.0f + cmt0);
    g1 = 1.0f / (1.0f + cmt1);
    g2 = 1.0f / (1.0f + cmt2);
    g3 = 1.0f / (1.0f + cmt3);
    
    // pre calc some often used terms
    t2g3 = t2 * g3;
    t1g2t2g3 = t1 * g2 * t2g3;
    t0g1 = t0 * g1;

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

    z0 = cmt0 / (1.0f + cmt0);
    z1 = cmt1 / (1.0f + cmt1);
    z2 = cmt2 / (1.0f + cmt2);
    z3 = cmt3 / (1.0f + cmt3);
    cgfbr = 1.0f / (1.0f + fb * z0*z1*z2*z3);

    // then solve the remaining outputs (with the non-linear gains here)
    xx = input - CLAMP(fb * estimate) * cgfbr;
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
            smp_t0 = tanh_5_4(smp[i]*gain);
            smp[i] = this->step((smp_t0 + smp_t1)/2.0); 
            smp_t1 = smp_t0;
            smp[i] *= outgain;
        }

        assert(smp[buffersize-1]<=10.0);
}

void MoogFilter::setfreq_and_q(float frequency, float q_)
{
    this->make_filter(frequency/sr, q_);
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

};
