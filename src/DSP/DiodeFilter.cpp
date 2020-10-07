#include <cassert>
#include <cstdio>
#include <cmath>
#include <stdio.h>

#include "DiodeFilter.h"

namespace zyn{

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )

DiodeFilter::DiodeFilter(unsigned char Ftype, float Ffreq, float Fq,
    unsigned char Fstages,
    unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f), stages(Fstages), type(Ftype)
{
    df = new DiodeLadderFilter;
    setfreq_and_q(Ffreq, Fq);
}

DiodeFilter::~DiodeFilter(void)
{
    delete df;

}

void DiodeFilter::filterout(float *smp)
{
    for (int i = 0; i < buffersize; i ++)
        {
            fc = (ff != ff_old) ? ((i*ff) + ((buffersize-i)* ff_old)) / buffersize : ff;
            smp[i] = df->tick(LIMIT(smp[i]*gain, 0.8, -0.8 ), fc*2); 
            smp[i] *= outgain;
        }
        ff_old = ff;
        
}

void DiodeFilter::setfreq_and_q(float frequency, float q_)
{
    setfreq(frequency);
    setq(q_);
}

void DiodeFilter::setfreq(float frequency)
{
    ff = LIMIT(frequency,0.4999f,0.0001f);
}

void DiodeFilter::setq(float q_)
{
    df->set_q(sqrtf(sqrtf(q_/1000))*1.0);
}

void DiodeFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}

void DiodeFilter::settype(unsigned char type_)
{
    type = type_;
}

void DiodeFilter::setstages(unsigned char stages_)
{

}

};
