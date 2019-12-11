#include <cstdio>

#include <cmath>

#include "../Misc/Matrix.h"
#include "MoogFilter.h"

namespace zyn{

// `Cheap non-linear zero-delay filters',
// version from aciddose
// https://www.kvraudio.com/forum/viewtopic.php?f=33&t=349859&start=285 (function cascade_4)

#define LIMIT(x,lu, ll) (x>lu ? lu : (x<ll ? ll : x) )
#define CLAMP(x, f) (x/f / powf(1.0+powf(fabs(x/f), 12.0), 1.0/12.0)*f)

inline float MoogFilter::tanhd(const float x, const float d = 1.0f, const float s=0.1f)
	{
		return 1.0f - s * (d + 1.0f) * x*x / (d + x*x);
		//return tanh(x) /x;
	}


float MoogFilter::step(float input)
{
   // per-sample computation
   t0 = tanhd(b[0] + a, d, s);
   t1 = tanhd(b[1] + a, d, s);
   t2 = tanhd(b[2] + a, d, s);
   t3 = tanhd(b[3] + a, d, s);

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

    //b[0] = LIMIT(b[0], 1.0f, -1.0f);

   estimate =
       g3 * b[3] +
       f3 * g2 * b[2] +
       f2 * g1 * b[1] +
       f1 * g0 * b[0] +
       f0 * input;

    // feedback gain coefficient, absolutely critical to get this correct
    // i believe in the original this is computed incorrectly?
    cgfbr = 1.0f / (1.0f + fb * z0*z1*z2*z3);

    // clamp can be a hard clip, a diode + highpass is better
    // if you implement a highpass do not forget to include it in the computation of the gain coefficients!
    xx = input - CLAMP(fb * estimate, 1.0f) * cgfbr;
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

void MoogFilter::make_filter(float ff, float q)
{
				ff = LIMIT(ff,0.41f,0.0001f); // makes overflow if higher
				f = tan(3.64f * ff); // tunes to AnalogFilter
        c = f;
        c2 = 2.0 * f;
        fb = cbrtf(q/1000)*10.0 - 0.45; // flattening
        compensation = 1.0f + LIMIT(fb, 1.0, 0);
        printf("fb = %f\n", fb);
				printf("ff = %f\n", ff);
}

std::vector<float> MoogFilter::impulse_response(float alpha, float k)
{
    this->make_filter(alpha, k);

    int ir_len = 10000;

    std::vector<float> output;
    output.push_back(step(1.0f));
    for(int i=0; i<ir_len-1; ++i)
        output.push_back(this->step(0.0f));


    make_filter(2*alpha, k);
    for(int i=0; i<ir_len-1; ++i)
        output.push_back(this->step(0.0f));

    return output;
}

MoogFilter::MoogFilter(float Ffreq, float Fq,
        unsigned char non_linear_element,
        unsigned int srate, int bufsize)
    :Filter(srate, bufsize), sr(srate), gain(1.0f)
{
    (void) non_linear_element; // TODO

    this->make_filter(Ffreq/srate, Fq);

}

MoogFilter::~MoogFilter(void)
{

}
void MoogFilter::filterout(float *smp)
{
    for(int i=0; i<buffersize; ++i)
        smp[i] = this->step(gain*smp[i]);
}
void MoogFilter::setfreq(float /*frequency*/)
{
    //Dummy
}

void MoogFilter::setfreq_and_q(float frequency, float q_)
{
    this->make_filter(frequency/sr, q_);

}

void MoogFilter::setq(float /*q_*/)
{
}

void MoogFilter::setgain(float dBgain)
{
    gain = dB2rap(dBgain);
}
};
