/*
  ZynAddSubFX - a software synthesizer

  SEQParams.h - Parameters for LFO
  Copyright (C) 2021 Michael Kirchner
  Author: Michael Kirchner

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.
*/

#ifndef SEQ_PARAMS_H
#define SEQ_PARAMS_H

#include <Misc/Time.h>
#include <rtosc/ports.h>
#include "Presets.h"

#define MAX_CUTOFF 40.0f

namespace zyn {

class XMLwrapper;

class SEQParams:public Presets
{
    public:
        SEQParams(const AbsTime* time_ = nullptr);
        SEQParams(consumer_location_t loc,
                  const AbsTime* time_ = nullptr);
        SEQParams(float freq_,
                  float cutoff_,
                  unsigned char steps_,
                  float delay_,
                  bool continous,
                  consumer_location_t loc,
                  const AbsTime* time_ = nullptr);
        ~SEQParams() override;

        void add2XML(XMLwrapper& xml) override;
        void defaults();
        /**Loads the SEQ from the xml*/
        void getfromXML(XMLwrapper& xml);
        void paste(SEQParams &);


        float freq=120.0f;
        float cutoff=0.0f; /**<cutoff frequency of LP filter (0.0f=off) */
        float delay=0.0f; /**<delay (0=off)*/
        float intensity=0.0f;
        unsigned char continous=0;
        unsigned char steps=8;
        int           numerator=0;  /**<numerator for integer ratio between system tempo and LFO freq (0=off)*/
        int           denominator=4;/**<denominator for integer ratio between system tempo and LFO freq (0=off)*/
        float         sequence[NUM_SEQ_STEPS]={};

        //! what kind is the SEQ (0 - frequency, 1 - amplitude, 2 - filter)
        consumer_location_type_t fel;
        int loc; //!< consumer location

        const AbsTime *time;
        int64_t last_update_timestamp=0;

        static const rtosc::Ports &ports;
    private:
        //! common functionality of ctors
        void setup();

        /* Default parameters */
        float         Dfreq=120.0f;
        unsigned char Dcutoff=0.0f;
        float         Ddelay=0.0f;
        unsigned char Dcontinous=0;
        unsigned char Dsteps=0;
        float         Dsequence[NUM_SEQ_STEPS]={};
};

}

#endif
