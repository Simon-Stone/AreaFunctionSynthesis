// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __IIRFILTER_H__
#define __IIRFILTER_H__

#include <complex>
#include "Signal.h"

const int MAX_IIR_ORDER = 32;
const int IIR_BUFFER_MASK   = 63; 
const int IIR_BUFFER_LENGTH = 64;


typedef std::complex<double> ComplexValue;

// ****************************************************************************
/// This class represents a recursive infinit impulse response filter.
/// The transfer function has the form:
///
/// H(z) = [a0 + a1*z^(-1) + a2*z^(-2) + ...] / [1 - b1*z^(-1) - b2*z^(-2) - ...].
///
/// The recursion formula is thus:
///
/// y[n] = a0*x[n] + a1*x[n-1] + a2*x[n-2] + ... + b1*y[n-1] + b2*y[n-2] + ...
///
// ****************************************************************************

class IirFilter
{
  public:
    IirFilter();
    void resetBuffers(double initialValue = 0.0);

    double getOutputSample(double nextInputSample);
    ComplexValue getFrequencyResponse(double freqRatio);
    void getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength);
    void getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength, int SR, double F0);

    void setGain(double gain);
    void setCoefficients(const double *A, const double *B, const int newOrder);
    bool combineWithFilter(const IirFilter *f, bool cascade);

    // Creation of some simple filters.

    void createUnityFilter();
    void createSinglePoleLowpass(double cutoffFreqRatio);
    void createSinglePoleHighpass(double cutoffFreqRatio);
    void createSecondOrderLowpass(double freqRatio, double Q);
    void createChebyshev(double cutoffFreqRatio, bool isHighpass, int numPoles);

public:
    double a[MAX_IIR_ORDER+1];
    double b[MAX_IIR_ORDER+1];
    int order;

private:
    int pos;
    double inputBuffer[IIR_BUFFER_LENGTH];
    double outputBuffer[IIR_BUFFER_LENGTH];

    void clearCoefficients();
};


#endif
