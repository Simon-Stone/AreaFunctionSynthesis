// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __TL_MODEL_H__
#define __TL_MODEL_H__

#include "Matrix2x2.h"
#include "Signal.h"
#include "Dsp.h"
#include "Tube.h"
#include "Constants.h"


// ****************************************************************************
/// This class encapsulates the simulation of vocal tract acoustics in the 
/// frequency domain on the basis of the transmission line analogy with lumped
/// elements.
// ****************************************************************************

class TlModel
{
  // ************************************************************************
  // Public data.
  // ************************************************************************

public:
  enum RadiationType 
  { 
    NO_RADIATION, 
    PISTONINSPHERE_RADIATION, 
    PISTONINWALL_RADIATION, 
    PARALLEL_RADIATION,
    NUM_RADIATION_OPTIONS
  };

  enum SpectrumType  
  { 
    INPUT_IMPEDANCE, OUTPUT_IMPEDANCE, FLOW_SOURCE_TF, PRESSURE_SOURCE_TF, RADIATION
  };

  /// Options for the acoustic simulation.

  struct Options
  {
    RadiationType radiation;
    bool boundaryLayer;
    bool heatConduction;
    bool softWalls;
    bool hagenResistance;
    bool innerLengthCorrections;
    bool lumpedElements;
    bool paranasalSinuses;
    bool piriformFossa;
    bool staticPressureDrops;
  };

  Options options;
  Tube tube;

  // ************************************************************************
  // Public functions.
  // ************************************************************************

public:
  TlModel();
  
  void getImpulseResponseWindow(Signal *window, int length);
  void getImpulseResponse(Signal *impulseResponse, int lengthExponent);
  void getSpectrum(SpectrumType type, ComplexSignal *spectrum, int spectrumLength, int section);

  int getMostConstrictedSection();
  double getMeanFlow(double lungPressure_dPa);
  void setLungPressure(double lungPressure_dPa);
  void getFormants(double *formantFreq, double *formantBW, int &numFormants, 
    const int MAX_FORMANTS, bool &frictionNoise, bool &isClosure, bool &isNasal);

  static double getCircumference(double area);

  // ************************************************************************
  // Private data.
  // ************************************************************************

private:
  /// Maximal number of sampling points in the spectrum
  static const int    MAX_NUM_FREQ = 1024;
  static const double MIN_AREA_CM2;
  static const double MIN_FREQ_RAD;

  // Options and tube geometry used in the previous calculation
  Options prevOptions;
  Tube prevTube;

  /// The product of the tube section matrices within a branch.
  Matrix2x2 matrixProduct[Tube::NUM_SECTIONS][MAX_NUM_FREQ];

  bool resetCalculations;   ///< Must the calculations be reset, because some parameter has changed
  double f0;                ///< Current frequency resolution
  int numFreq;
  double lungPressure_dPa;   ///< Currently set lung pressure in Pa

  double discreteOmega[MAX_NUM_FREQ];
  Complex mouthRadiationImpedance[MAX_NUM_FREQ];
  Complex noseRadiationImpedance[MAX_NUM_FREQ];
  Complex lungTerminationImpedance[MAX_NUM_FREQ];
  Complex radiationCharacteristic[MAX_NUM_FREQ];


  // ************************************************************************
  // Private functions.
  // ************************************************************************

private:
  void prepareCalculations();

  Complex getRadiationCharacteristic(double omega);
  Complex getRadiationImpedance(double omega, double radiationArea_cm2);
  void getLumpedSectionImpedances(double omega, Tube::Section *ts, Complex &Za, Complex &Zb);
  Matrix2x2 getSectionMatrix(double omega, int section);
  Complex getJunctionImpedance(double omega, double A1_cm2, double A2_cm2);

  Complex getInputImpedance(int freqIndex, int section);
  Complex getOutputImpedance(int freqIndex, int section);
  Complex getPressureSourceTF(int freqIndex, int section);
  Complex getFlowSourceTF(int freqIndex, int section);
};

#endif
