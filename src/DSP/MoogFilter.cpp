#include <cstdio>
#include <cmath>
#include "MoogFilter.h"

namespace zyn{

// `Cheap non-linear zero-delay filters',
// version from aciddose
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=349859&start=285 (function cascade_4)

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )
#define CLAMP(x, e) (x / powf(1.0+powf(fabs(x), e), 1.0/e)) //linear between -1 and 1 with round edges

inline float MoogFilter::tanhd(const float x, const float d = 1.0f, const float s=0.1f)
	{
		return 1.0f - s * (d + 1.0f) * x*x / (d + x*x);
	}


float MoogFilter::step(float input)
{
   // per-sample computation
   t0 = tanhd(b[0] + a, 1.0, par1);
   t1 = tanhd(b[1] + a, 1.0, par1);
   t2 = tanhd(b[2] + a, 1.0, par1);
   t3 = tanhd(b[3] + a, 1.0, par1);

   g0 = 1.0f / (1.0f + c*t0);
   g1 = 1.0f / (1.0f + c*t1);
   g2 = 1.0f / (1.0f + c*t2);
   g3 = 1.0f / (1.0f + c*t3);

   z0 = c*t0 / (1.0f + c*t0);
   z1 = c*t1 / (1.0f + c*t1);
   z2 = c*t2 / (1.0f + c*t2);
   z3 = c*t3 / (1.0f + c*t3);

   f3 = c       * t2*g3;
   f2 = c*c     * t1*g2 * t2*g3;
   f1 = c*c*c   * t0*g1 * t1*g2 * t2*g3;
   f0 = c*c*c*c *    g0 * t0*g1 * t1*g2 * t2*g3;


   estimate =
       g3 * b[3] +
       f3 * g2 * b[2] +
       f2 * g1 * b[1] +
       f1 * g0 * b[0] +
       f0 * input;

    // feedback gain coefficient, absolutely critical to get this correct
    cgfbr = 1.0f / (1.0f + fb * z0*z1*z2*z3);

    // clamp can be a hard clip, a diode + highpass is better
    xx = input - CLAMP(fb * estimate, par3) * cgfbr;
    y0 = t0 * g0 * (b[0] + c * xx);
    y1 = t1 * g1 * (b[1] + c * y0);
    y2 = t2 * g2 * (b[2] + c * y1);
    y3 = t3 * g3 * (b[3] + c * y2);

    b[0] += c2 * (xx - y0);
    b[1] += c2 * (y0 - y1);
    b[2] += c2 * (y1 - y2);
    b[3] += c2 * (y2 - y3);

    // you must limit the compensation if feedback is clamped
    return y3 * compensation;
}

std::vector<float> MoogFilter::impulse_response(float alpha, float k)
{
    this->setfreq_and_q(alpha, k);
    int ir_len = 10000;
    std::vector<float> output;
    output.push_back(step(1.0f));
    for(int i=0; i<ir_len-1; ++i)
        output.push_back(this->step(0.0f));
    return output;
}

MoogFilter::MoogFilter(float Ffreq, float Fq,
        unsigned char non_linear_element,
        unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f)
{
    this->setfreq_and_q(Ffreq/srate, Fq);
}

MoogFilter::~MoogFilter(void)
{

}

void MoogFilter::filterout(float *smp)
{
    for(int i=0; i<buffersize; ++i)
        smp[i] = this->step(gain*smp[i]);
}

void MoogFilter::setfreq(float frequency/*frequency*/)
{
	frequency /= ( 1 - (1.49*par1) ); // compensate nonlinerarity impact
    c = LIMIT(frequency/sr,0.41f,0.0001f); // makes overflow if higher
	c = tan(3.16f * c); // tuned to AnalogFilter
	c2 = 2.0 * c; // precalc of often used values
}

void MoogFilter::setq(float q_/*q_*/)
{
	fb = cbrtf(q_/1000)*12.0 - 0.45; // flattening
	compensation = 1.3f + LIMIT(fb, 1.0, 0);
}

void MoogFilter::setfreq_and_q(float frequency, float q_)
{
    this->setfreq(frequency);
		this->setq(q_);
}

void MoogFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}
};
