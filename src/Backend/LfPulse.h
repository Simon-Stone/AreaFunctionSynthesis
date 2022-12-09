// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __LF_PULSE_H__
#define __LF_PULSE_H__

#include "Signal.h"

// ****************************************************************************
/// This class represents the model of glottal flow introduced by Liljencrants
/// and Fant.
// ****************************************************************************

class LfPulse
{
  // **************************************************************************
  // Pulse parameters.
  // **************************************************************************

public:
  double F0;      // in Hz
  double AMP;     // in cm^3/s
  double OQ;      // [0, 1] Open quotient
  double SQ;      // [1, 2] Speed quotient
  double TL;      // [0, 0.2] Spectral tilt

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  LfPulse();
  void resetParams();
  void getPulse(Signal& s, int numSamples, bool getDerivative);

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  double getEpsilon(double ta, double te);
  double getAlpha(double tp, double te, double ta, double epsilon);
  double getB(double AMP, double tp, double alpha);
};

#endif
