/*
  ZynAddSubFX - a software synthesizer

  Hysteresis.cpp - Hysteresis effect
  Copyright (C) 2002-2005 Nasca Octavian Paul
  Copyright (C) 2009-2010 Mark McCurry
  Author: Nasca Octavian Paul
          Mark McCurry

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#include <cmath>
#include <rtosc/ports.h>
#include <rtosc/port-sugar.h>
#include "../Misc/Allocator.h"
#include "Hysteresis.h"


namespace zyn {

#define rObject Hysteresis
#define rBegin [](const char *msg, rtosc::RtData &d) {
#define rEnd }

rtosc::Ports Hysteresis::ports = {
    {"preset::i", rOptions(Hysteresis 1, Hysteresis 2)
                  rProp(parameter)
                  rDoc("Instrument Presets"), 0,
                  rBegin;
                  rObject *o = (rObject*)d.obj;
                  if(rtosc_narguments(msg))
                      o->setpreset(rtosc_argument(msg, 0).i);
                  else
                      d.reply(d.loc, "i", o->Ppreset);
                  rEnd},
    rEffParVol(rDefault(64)),
    rEffParPan(),
    rEffPar(Pdrive,   2, rShort("drive"),
            "Drive of Hysteresis"),
    rEffPar(Premanence,   3, rShort("rem"),
            "Remanence of Hysteresis"),
    rEffPar(Pcoercivity,  4, rShort("coerc"),
            "Coercivity of Hysteresis"),

};
#undef rBegin
#undef rEnd
#undef rObject

Hysteresis::Hysteresis(EffectParams pars)
    :Effect(pars),
      Pvolume(50),
      remanence(0.5f),
      coercivity(0.5f)
{
    //~ hyst_l = memory.alloc<JilesAtherton>(float(Pcoercivity)*8.0f, float(Premanence)/128.0f);
    //~ hyst_r = memory.alloc<JilesAtherton>(float(Pcoercivity)*8.0f, float(Premanence)/128.0f);
}

Hysteresis::~Hysteresis()
{
    //~ memory.devalloc(hyst_l);
    //~ memory.devalloc(hyst_r);
}

//Initialize the delays
void Hysteresis::init(void)
{
    //~ hyst_l->init();
    //~ hyst_r->init();
}
#define ALPHA 0.1f
#define BETA 0.2f
void Hysteresis::out(const Stereo<float *> &input)
{
    for(int i = 0; i < buffersize; ++i) {
        
        // Vorverarbeitung des Eingangssignals unter Berücksichtigung von Remanenz und Koerzitivität
        float processed_input_l = input.l[i];
        if (fabs(input.l[i]) > coercivity) {
            if(input.l[i] > 0)
                processed_input_l = input.l[i] + remanence * tanh(drive * input.l[i] - coercivity);
            else
                processed_input_l = input.l[i] + remanence * tanh(drive * input.l[i] + coercivity);
        }
        float processed_input_r = input.r[i];
        if (fabs(input.r[i]) > coercivity) {
            if(input.r[i] > 0)
                processed_input_r = input.r[i] + remanence * tanh(input.r[i] - coercivity);
            else
                processed_input_r = input.r[i] + remanence * tanh(input.r[i] + coercivity);
        }


       // Aktualisierung der Zustände mit dem vorverarbeiteten Eingang
       state_l = ALPHA * state_l + BETA * tanh(processed_input_l);
       state_r = ALPHA * state_r + BETA * tanh(processed_input_r);
       efxoutl[i] = state_l;
       efxoutr[i] = state_r;
    }
}


//Parameter control
void Hysteresis::setvolume(unsigned char _Pvolume)
{
    Pvolume = _Pvolume;

    if(insertion == 0) {
        if (Pvolume == 0) {
            outvolume = 0.0f;
        } else {
            outvolume = powf(0.01f, (1.0f - Pvolume / 127.0f)) * 4.0f;
        }
        volume    = 1.0f;
    }
    else
        volume = outvolume = Pvolume / 127.0f;
}

void Hysteresis::setdrive(unsigned char Pdrive)
{
    drive   = Pdrive / 64.0f;
}

void Hysteresis::setremanence(unsigned char Premanence)
{
    remanence   = Premanence / 127.0f;
    //~ hyst_l->setMr(float(Premanence+3)/32.0f);
    //~ hyst_r->setMr(float(Premanence+3)/32.0f);
}

void Hysteresis::setcoercivity(unsigned char Pcoercivity)
{
    coercivity   = Pcoercivity / 1270.0f;
    //~ hyst_l->setHc(float((Pcoercivity+2)/8.0f));
    //~ hyst_r->setHc(float((Pcoercivity+2)/8.0f));
}

// Add the setter functions for alpha and beta
//~ void Hysteresis::setalpha(unsigned char Palpha) {
    //~ hyst_l->setAlpha(float(Palpha)/128.0f);
    //~ hyst_r->setAlpha(float(Palpha)/128.0f);
//~ }

//~ void Hysteresis::setbeta(unsigned char Pbeta) {
    //~ hyst_l->setBeta(float(Pbeta)/128.0f);
    //~ hyst_r->setBeta(float(Pbeta)/128.0f);
//~ }

unsigned char Hysteresis::getpresetpar(unsigned char npreset, unsigned int npar)
{
#define	PRESET_SIZE 3
#define	NUM_PRESETS 1
    static const unsigned char presets[NUM_PRESETS][PRESET_SIZE] = {
        {67, 64, 35 }, //Hysteresis 1

    };
    if(npreset < NUM_PRESETS && npar < PRESET_SIZE) {
        if(npar == 0 && insertion != 0) {
            /* lower the volume if this is insertion effect */
            return presets[npreset][npar] / 2;
        }
        return presets[npreset][npar];
    }
    return 0;
}

void Hysteresis::setpreset(unsigned char npreset)
{
    if(npreset >= NUM_PRESETS)
        npreset = NUM_PRESETS - 1;
    for(int n = 0; n != 128; n++)
        changepar(n, getpresetpar(npreset, n));
    Ppreset = npreset;
}

void Hysteresis::changepar(int npar, unsigned char value)
{
    switch(npar) {
        case 0:
            setvolume(value);
            break;
        case 1:
            setpanning(value);
            break;
        case 2:
            setdrive(value);
            break;
        case 3:
            setremanence(value);
            break;
        case 4:
            setcoercivity(value);
            break;

    }
}

unsigned char Hysteresis::getpar(int npar) const
{
    switch(npar) {
        case 0:  return Pvolume;
        case 1:  return Ppanning;
        case 2:  return int(drive*64.0f);
        case 3:  return int(remanence*127.0f);
        case 4:  return int(coercivity*1270.0f);
        default: return 0; // in case of bogus parameter number
    }
}

}