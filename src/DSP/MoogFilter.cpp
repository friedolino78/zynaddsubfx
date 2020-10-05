#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "MoogFilter.h"

namespace zyn{

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )
#define NLEstage(x) tanhd(x)
//~ #define NLEstage(x) softlimiter8 (x, 1.0, 0.01, 0.8)
#define NLEfb(x) tanhd(x)
//~ #define NLEfb(x) softlimiter8 (x, 2.0, 0.05, 0.9)

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
    fb = sqrtf(sqrtf(q/1000))*4.0; // flattening
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
    float x = 0.95 * a; // reduce drive of the nle to reduce aliasing
    x2 = x*x; 
    return((10 * x2*x + 105*x)/((x2*x2 + 45*x2 + 105))); // Pade approximation of tanh; Cest 0.88
}

inline float MoogFilter::smoothABS ( float x, const float y) // y controls 'smoothness' usually between 0.002 -> 0.04
{
return (sqrtf((x * x)  + y)) - sqrtf(y);
}

inline float MoogFilter::smoothclip (float x, const float a, const float b) // assuming symmetrical clipping
{
  float  x1 = smoothABS (x, a);
  float  x2 = smoothABS (x, b);
   x = x1 + (a+b);
   x = x - x2;
   x = x * 0.5;
   return x;
}

// // tanh soft limiter
inline float MoogFilter::tanhsoftlimiter (float x, const float gain, const float offset, const float par) // assuming symmetrical clipping
{
    x *= gain;// multiply signal to drive it in the saturation of the function
    x += offset; // add dc offset
    x = x / powf(1.0+powf(fabsf(x), par), 1.0/par);
    x -= offset / powf(1+powf(fabsf(offset), par), 1.0/par);
    return x;
}

// // tanh soft limiter
inline float MoogFilter::softlimiter8 (float x, const float drive, const float offset, const float ammount) // assuming symmetrical clipping
{
    x *= drive;// multiply signal to drive it in the saturation of the function
    x += offset; // add dc offset
    static float x2 = x*x;
    static float x4 = x2*x2;
    static float x8 = x4*x4;
    static float y = x / sqrtf(sqrtf(sqrtf((1.0+x8))));
//    x -= offset / powf(1+powf(fabsf(offset), par), 1/par);
    return (((1-ammount) * x) + (ammount * y)*0.95); // lerp with input
}

// // tanh soft limiter
inline float MoogFilter::softlimiter4 (float x, const float drive, const float offset, const float ammount) // assuming symmetrical clipping
{
    x *= drive;// multiply signal to drive it in the saturation of the function
    x += offset; // add dc offset
    static float x2 = x*x;
    static float x4 = x2*x2;
    static float y = x / sqrtf(sqrtf((1.0+x4)));
//    x -= offset / powf(1+powf(fabsf(offset), par), 1/par);
    return (((1-ammount) * x) + (ammount * y)*0.95); // lerp with input
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

    t0 = NLEstage(b[0] + a);
    t1 = NLEstage(b[1] + a);
    t2 = NLEstage(b[2] + a);
    t3 = NLEstage(b[3] + a);

    // denominators for solutions of individual stages
    g0 = 1.0f / (1.0f + c*t0);
    g1 = 1.0f / (1.0f + c*t1);
    g2 = 1.0f / (1.0f + c*t2);
    g3 = 1.0f / (1.0f + c*t3);

    // pre calculate some values used mutliple times
    t2g3 = t2 * g3;
    t0g1 = t0 * g1;
    t1g2t2g3 = t1 * g2 * t2g3;

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
    xx = input - fb * NLEfb(estimate) * cgfbr;
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
            smp[i] = this->step(LIMIT(smp[i]*gain, 0.8, -0.8 )); 
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
