#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "../Misc/Util.h"
#include "MoogFilter.h"


namespace zyn{

MoogFilter::MoogFilter(unsigned char Ftype, float Ffreq, float Fq,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f), type(Ftype)
{
    setfreq_and_q(Ffreq/srate, Fq);
    settype(type); // q must be set before
    for (int i = 0; i>4; i++)
    {
        b[i] = 0.0f;
    }
}

MoogFilter::~MoogFilter(void)
{

}

inline float MoogFilter::tan_2(const float x)
{    
    //Pade approximation tan(x) hand tuned to map fCutoff 
    x2 = x*x;
    return ((4.75f*(22.3f*x - 2.0f*x2*x))/(105.0f - 45.0f*x2 + x2*x2));
}

inline float MoogFilter::tanhX(const float x)
{   
    // Pade approximation of tanh(x) used for input saturation
    // https://mathr.co.uk/blog/2017-09-06_approximating_hyperbolic_tangent.html
    x2 = x*x;
    return (x*(105+10*x2)/(105+(45+x2)*x2)); // bound to [-1 .. +1]
}


inline float MoogFilter::tanhXdX(float x)
{
    // Pade approximation for tanh(x)/x used in filter stages
    x2 = x*x;
    return ((x2 + 105.0f)*x2 + 945.0f) / ((15.0f*x2 + 420.0f)*x2 + 945.0f);
}

inline float MoogFilter::step(float input)
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
    // based on the version from aciddose
    // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=349859&start=285 (function cascade_4)
    // reading recommendation https://www.native-instruments.com/fileadmin/ni_media/downloads/pdf/VAFilterDesign_2.0.0a.pdf
    //
    // modified for performance (precalculations + tanh approximation) by Friedolino78

    // evaluate the transconductance gM(vIn) = iCtl * ( tanh( vIn ) / vIn )
    t0 = tanhXdX(b[0]);
    t1 = tanhXdX(b[1]);
    t2 = tanhXdX(b[2]);
    t3 = tanhXdX(b[3]);

    // pre calc often used terms
    cmt0 = c*t0;
    cmt1 = c*t1;
    cmt2 = c*t2;
    cmt3 = c*t3;
    
    // denominators for solutions of individual stages
    g0 = 1.0f / (1.0f + cmt0);
    g1 = 1.0f / (1.0f + cmt1);
    g2 = 1.0f / (1.0f + cmt2);
    g3 = 1.0f / (1.0f + cmt3);
    
    // pre calc often used term
    t1g2t2g3 = t1 * g2 * t2 * g3;

    // factored out of the feedback solution
    f3 = c * t2 * g3;
    f2 = cp2 * t1g2t2g3;
    f1 = cp3 * t0 * g1 * t1g2t2g3;
    f0 = cp4 * g0 * t0 * g1 * t1g2t2g3;

    // estimate feedback for t0
    estimate =
        g3 * b[3] +
        f3 * g2 * b[2] +
        f2 * g1 * b[1] +
        f1 * g0 * b[0] +
        f0 * input;

    // feedback gain coefficient
    z0 = cmt0 / (1.0f + cmt0);
    z1 = cmt1 / (1.0f + cmt1);
    z2 = cmt2 / (1.0f + cmt2);
    z3 = cmt3 / (1.0f + cmt3);
    cgfbr = 1.0f / (1.0f + fb * z0*z1*z2*z3);


    // calculate input for the fist stage
    xx = input - tanhX(fb * estimate) * cgfbr;
    // calculate output of all stages
    y0 = t0 * g0 * (b[0] + c * xx);
    y1 = t1 * g1 * (b[1] + c * y0);
    y2 = t2 * g2 * (b[2] + c * y1);
    y3 = t3 * g3 * (b[3] + c * y2);

    // update state
    b[0] += cm2 * (xx - y0);
    b[1] += cm2 * (y0 - y1);
    b[2] += cm2 * (y1 - y2);
    b[3] += cm2 * (y2 - y3);

    // calculate multimode filter output
    return (a0 * xx
          + a1 * y0 
          + a2 * y1 
          + a3 * y2 
          + a4 * y3);
}

void MoogFilter::filterout(float *smp)
{
    for (int i = 0; i < buffersize; i ++)
    {
        smp[i] = step(tanhX(smp[i]*gain));
        smp[i] *= outgain;
    }
}

void MoogFilter::setfreq_and_q(float frequency, float q_)
{
    setfreq(frequency/sr);
    setq(q_);
}

void MoogFilter::setfreq(float ff)
{
    // limit cutoff to prevent overflow
    ff = limit(ff,0.0002f,0.48f);
    // pre warp cutoff to map to reality
    c = tan_2(PI * ff);
    // pre calculate some stuff outside the hot zone
    cm2 = c * 2.0;
    cp2 = c * c;
    cp3 = cp2 * c;
    cp4 = cp2 * cp2;
}

void MoogFilter::setq(float q)
{
    // flattening the Q input
    fb = cbrtf(q/1000.0f)*4.0f + 0.1f; 
    // compensation factor for passband reduction by the negative feedback
    compensation = 1.0f + limit(fb, 0.0f, 1.0f); 
}

void MoogFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}

void MoogFilter::settype(unsigned char type_)
{
    type = type_;
    // set coefficients for each filter type
    // theory from "THE ART OF VA FILTER DESIGN"
    // by Vadim Zavalishin
    switch (type) 
    {
        case 1:
            a0 = 0.0; a1 = 0.0; a2 = 4.0; a3 =-8.0; a4 = 4.0;
            break;
        case 0:
            a0 = 1.0; a1 =-4.0; a2 = 6.0; a3 =-4.0; a4 = 1.0;
            break;
        case 2:
        default:
            a0 = 0.0; a1 = 0.0; a2 = 0.0; a3 = 0.0; a4 = compensation;
            break;
    }
}

};
