
// This code is released under the MIT license (see below).
//
// The MIT License
// 
// Copyright (c) 2012 Dominique Wurtz (www.blaukraut.info)
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef __DIODE_LADDER_FILTER_HPP__
#define __DIODE_LADDER_FILTER_HPP__

#include <cmath>
#include <algorithm>

// Emulation of Diode ladder lowpass filter as found in Roland TB303 or EMS VCS3
// Version 0.1 (04/03/2012)
namespace zyn{
	
class DiodeLadderFilter
{
public:

	DiodeLadderFilter()
	{
		std::fill(z, z + 4, 0);
		set_q(0);
	}

	void reset()
	{
		if (k < 17) std::fill(z, z + 4, 0);
	}

	// q: resonance in the range [0..1]
	void set_q(const double q)
	{
		assert(q >= 0 && q <= 1.);
		k = 20 * q;
		A = 1 + 0.5*k; // resonance gain compensation 
	}

	// Process one sample.
	//
	// x: input signal
	// fc: normalized cutoff frequency in the range [0..1] => 0 HZ .. Nyquist
	double tick(const double x, const double fc)
	{
		assert(fc > 0 && fc < 1);
		const double wc = 2 * tan(0.5* PI_HALF * fc); // PI is Nyquist frequency 
		//~ wc = 2 * tan(0.5*wc); // dewarping, not required with 2x oversampling
		const double wc2 = wc*wc;
		const double wc3 = wc2*wc;
		const double wc4 = wc3*wc;
		const double b = 1 / (1+8*wc+20*wc2+16*wc3+2*wc4);
		const double g = 2*wc4 * b;

		// current state
		const double s = (z[0]*wc3 + z[1]*(wc2+2*wc3) + z[2]*(wc+4*wc2+2*wc3) + z[3]*(1+6*wc+9*wc2+2*wc3)) * b;
		
		// solve feedback loop (linear)
		double y4 = (g*x + s) / (1 + g*k);

		// input clipping
		const double y0 = fast_tanh(x - k*y4);

		// Compute all integrator outputs (y1, y2, y3, y4).
		// Unlike in the well-known Moog transistor ladder, this gets quite nasty due the
		// inherent coupling between filter stages.
		const double y1 = (y0*(2*wc+12*wc2+20*wc3+8*wc4) + z[0]*(1+6*wc+10*wc2+4*wc3) +
			z[1]*(2*wc+8*wc2+6*wc3) + z[2]*(2*wc2+4*wc3) + z[3]*2*wc3)*b;
		const double y2 = (y0*(2*wc2+8*wc3+6*wc4) + z[0]*(wc+4*wc2+3*wc3) +
			z[1]*(1+6*wc+11*wc2+6*wc3) + z[2]*(wc+4*wc2+4*wc3) + z[3]*(wc2+2*wc3))*b;
		const double y3 = (y0*(2*wc3+4*wc4) + z[0]*(wc2+2*wc3) +
			z[1]*(wc+4*wc2+4*wc3) + z[2]*(1+6*wc+10*wc2+4*wc3) + z[3]*(wc+4*wc2+2*wc3))*b;
		y4 = g*y0 + s;

		// update filter state
		z[0] += 4*wc*(y0 - y1 + y2);
		z[1] += 2*wc*(y1 - 2*y2 + y3);
		z[2] += 2*wc*(y2 - 2*y3 + y4);
		z[3] += 2*wc*(y3 - 2*y4);

		return A*y4;
	}
	
private:
	double k, A;
	double z[4];
	double PI_HALF = M_PI/2.0;

	static inline double fast_tanh(const double x)
	{
		return x / (1 + abs(x));
	}
};
}
#endif // __DIODE_LADDER_FILTER_HPP__
