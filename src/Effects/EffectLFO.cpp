/*
  ZynAddSubFX - a software synthesizer

  EffectLFO.cpp - Stereo LFO used by some effects
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Author: Nasca Octavian Paul

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include "EffectLFO.h"
#include "../Misc/Util.h"

#include <cmath>
#include "globals.h"

namespace zyn {

MonoEffectLFO::MonoEffectLFO(float srate_f, float bufsize_f)
  : phase(0.0f),
    ampl1(RND),
    ampl2(RND),
    lfornd(0.0f),
    samplerate_f(srate_f),
    buffersize_f(bufsize_f) {
updateparams(40, 0, 0);
}

MonoEffectLFO::~MonoEffectLFO() {}

// Update the changed parameters
void MonoEffectLFO::updateparams(unsigned char Pfreq, unsigned char Prandomness, unsigned char PLFOtype_) {
    float lfofreq = (powf(2.0f, Pfreq / 127.0f * 10.0f) - 1.0f) * 0.03f;
    PLFOtype = PLFOtype_;
    inc = fabsf(lfofreq) * buffersize_f / samplerate_f;
    if (inc > 0.49999999f) {
      inc = 0.499999999f; // Limit the Frequency
    }

    lfornd = Prandomness / 127.0f;
    lfornd = (lfornd > 1.0f) ? 1.0f : lfornd;

    if (PLFOtype > 1) {
      PLFOtype = 1; // this has to be updated if more lfo's are added
    }
}

// Compute the shape of the LFO
float MonoEffectLFO::getlfoshape(float x) {
    x = fmodf(x, 1.0f);
    float out;
    switch (PLFOtype) {
      case 1: // EffectLFO_TRIANGLE
        if ((x > 0.0f) && (x < 0.25f)) {
          out = 4.0f * x;
        } else if ((x > 0.25f) && (x < 0.75f)) {
          out = 2.0f - 4.0f * x;
        } else {
          out = 4.0f * x - 4.0f;
        }
        break;
      // when adding more, ensure ::updateparams() gets updated
      default:
        out = cosf(x * 2.0f * M_PI); // EffectLFO_SINE
    }
    return out;
}

// LFO output
float MonoEffectLFO::effectlfoout(float phaseOffset, float stereoOffset) {
    float out = getlfoshape(phase + phaseOffset+stereoOffset);
    if ((PLFOtype == 0) || (PLFOtype == 1)) {
      out *= (ampl1 + phase * (ampl2 - ampl1));
    }
    out = (out + 1.0f) * 0.5f;

    // update phase for master lfo
    if (phaseOffset == 0.0f) {
      phase += inc;
      if (phase > 1.0f) {
        phase -= 1.0f;
        ampl1 = ampl2;
        ampl2 = (1.0f - lfornd) + lfornd * RND;
      }
    }

    return out;
}

EffectLFO::EffectLFO(float srate_f, float bufsize_f)
    :Pfreq(40),
      Prandomness(0),
      PLFOtype(0),
      Pstereo(64),
      left(srate_f, bufsize_f),
      right(srate_f, bufsize_f)
{
    updateparams();
}

EffectLFO::~EffectLFO() {
    }

//Update the changed parameters
void EffectLFO::updateparams(void)
{
    left.updateparams(Pfreq, Prandomness, PLFOtype);
    right.updateparams(Pfreq, Prandomness, PLFOtype);
    stereoOffset = ((Pstereo - 64.0f) / 127.0f) + 1.0f;
}

//LFO output
void EffectLFO::effectlfoout(float *outl, float *outr, float phaseOffset)
{
    *outl = left.effectlfoout(phaseOffset);
    *outr = right.effectlfoout(phaseOffset, stereoOffset);

}

}
