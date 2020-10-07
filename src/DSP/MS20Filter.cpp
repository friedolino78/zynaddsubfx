#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "MS20Filter.h"

namespace zyn{

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )
#define NLEstage(x) tanhd(x)
//~ #define NLEstage(x) softlimiter8 (x, 1.0, 0.01, 0.8)
#define NLEfb(x) tanhd(x)
//~ #define NLEfb(x) softlimiter8 (x, 2.0, 0.05, 0.9)

#define sign(x) x<0 ? -1 : 1 

MS20Filter::MS20Filter(unsigned char Ftype, float Ffreq, float Fq,
    unsigned char Fstages,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f), stages(Fstages), type(Ftype)
{
    d1 = 0.0; d2 = 0.0;
    y1 = 0.0; y2 = 0.0;
}

MS20Filter::~MS20Filter(void)
{

}

void MS20Filter::make_filter(float ff, float q)
{
      f =  2.5 * ff * PI;
      h  = f / 1.0; // relative frequency with oversampling
      hh = 0.5 * h;
      fb = sqrtf(sqrtf(q/1000))*1.2; // flattening
      k  = 2*fb;// - 0.2*reso*freq;
      
}

inline float MS20Filter::tanhdx(const float x, const float d = 1.0f, const float s=0.1f)
{    
    return 1.0f - s * (d + 1.0f) * x * x / (d + x * x); // Cest: 1.69
}

inline float MS20Filter::tanhd(const float a)
{    
    float x = 0.95 * a; // reduce drive of the nle to reduce aliasing
    x2 = x*x; 
    return((10 * x2*x + 105*x)/((x2*x2 + 45*x2 + 105))); // Pade approximation of tanh; Cest 0.88
}

inline float MS20Filter::smoothABS ( float x, const float y) // y controls 'smoothness' usually between 0.002 -> 0.04
{
return (sqrtf((x * x)  + y)) - sqrtf(y);
}

inline float MS20Filter::smoothclip (float x, const float a, const float b) // assuming symmetrical clipping
{
  float  x1 = smoothABS (x, a);
  float  x2 = smoothABS (x, b);
   x = x1 + (a+b);
   x = x - x2;
   x = x * 0.5;
   return x;
}

// // tanh soft limiter
inline float MS20Filter::tanhsoftlimiter (float x, const float gain, const float offset, const float par) // assuming symmetrical clipping
{
    x *= gain;// multiply signal to drive it in the saturation of the function
    x += offset; // add dc offset
    x = x / powf(1.0+powf(fabsf(x), par), 1.0/par);
    x -= offset / powf(1+powf(fabsf(offset), par), 1.0/par);
    return x;
}

// // tanh soft limiter
inline float MS20Filter::softlimiter8 (float x, const float drive, const float offset, const float ammount) // assuming symmetrical clipping
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
inline float MS20Filter::softlimiter4 (float x, const float drive, const float offset, const float ammount) // assuming symmetrical clipping
{
    x *= drive;// multiply signal to drive it in the saturation of the function
    x += offset; // add dc offset
    static float x2 = x*x;
    static float x4 = x2*x2;
    static float y = x / sqrtf(sqrtf((1.0+x4)));
//    x -= offset / powf(1+powf(fabsf(offset), par), 1/par);
    return (((1-ammount) * x) + (ammount * y)*0.95); // lerp with input
}


inline float MS20Filter::f_dg(float s){
    return abs(s) > 1 ? 0.25 : 1.0;
}

inline float MS20Filter::f_g(float s){
    return fabs(s) > 1.0 ? s - 0.75*sign(s)*(fabs(s)-1.0) : s;
}

float MS20Filter::step(float x)
{
    // https://github.com/JoepVanlier/JSFX/blob/master/Basics/MS-20.jsfx (MIT License)
    float gd2k, ky2, gky2, dgky2,
          f1, f2, a, b, c, d, norm, sig1, thsig1, thsig1sq, sig2, thsig2, thsig2sq, tanhterm1, tanhterm2, hhthsig1sqm1, hhthsig2sqm1;

    gd2k = f_g(d2*k);

    tanhterm1 = tanhd(-d1 +  x - gd2k);
    tanhterm2 = tanhd(d1 - d2 + gd2k);

    for(int i = 0; i<4; i++)
    {
        ky2 = k*y2;
        //~ gky2 = f_g(ky2);
        gky2 = softlimiter8(ky2, 1.0, 0.0, 0.9);
        dgky2 = f_dg(ky2);

        sig1 = x - y1 - gky2;
        thsig1 = tanh(sig1);
        thsig1sq = thsig1 * thsig1;
      
        sig2 = y1 - y2 + gky2;
        thsig2 = tanh(sig2);
        thsig2sq = thsig2 * thsig2;
        hhthsig1sqm1 = hh*(thsig1sq - 1);
        hhthsig2sqm1 = hh*(thsig2sq - 1);
      
        f1 = y1 - d1 - hh*(tanhterm1 + thsig1);
        f2 = y2 - d2 - hh*(tanhterm2 + thsig2);
        a = -hhthsig1sqm1 + 1;
        b = -k*hhthsig1sqm1*dgky2;
        c = hhthsig2sqm1;
        d = (k*dgky2 - 1)*hhthsig2sqm1 + 1;
        
        norm = 1 / ( a*d - b*c );
        y1 = y1 - ( d*f1 - b*f2 ) * norm;
        y2 = y2 - ( a*f2 - c*f1 ) * norm;
    }
      
    d1 = y1;
    d2 = y2;
    return y2;
}

void MS20Filter::filterout(float *smp)
{
    for (int i = 0; i < buffersize; i ++)
        {
            // limit input amplitude to prevent overflow at high input levels + high gain
            // with res, gain and all input levels at max, it may still overflow
            smp[i] = this->step(LIMIT(smp[i]*gain, 0.85, -0.85 )); 
            smp[i] *= outgain;
        }
}

void MS20Filter::setfreq_and_q(float frequency, float q_)
{
    this->make_filter(frequency/sr, q_); // for type cheap
}

void MS20Filter::setfreq(float frequency)
{

}

void MS20Filter::setq(float q_)
{

}

void MS20Filter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}

void MS20Filter::settype(unsigned char type_)
{
    type = type_;
}

void MS20Filter::setstages(unsigned char stages_)
{

}

};
