// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __TRIANGULAR_GLOTTIS_H__
#define __TRIANGULAR_GLOTTIS_H__

#include "Glottis.h"
#include "IirFilter.h"


// ****************************************************************************
// This class defines the triangular-glottis two-mass-model (my invention).
// ****************************************************************************

class TriangularGlottis : public Glottis
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:

  enum ControlParamIndex 
  { 
    // Frequency and lung pressure always must be the first two parameters
    FREQUENCY, 
    PRESSURE,
    // Rest displacement of the vocal processes from the midline
    REST_DISP_1,
    REST_DISP_2,
    ARY_AREA, 
    ASPIRATION_STRENGTH,
    NUM_CONTROL_PARAMS  
  };

  enum StaticParamIndex
  {
    REST_LENGTH,
    REST_THICKNESS_1,
    REST_THICKNESS_2,
    MASS_1,
    MASS_2,
    DAMPING_RATIO_1,
    DAMPING_RATIO_2,
    SPRING_K_1,
    SPRING_K_2,
    CONTACT_SPRING_K_1,
    CONTACT_SPRING_K_2,
    COUPLING_SPRING_K,
    INLET_LENGTH,
    OUTLET_LENGTH,
    // The "natural F0" is the F0, when the tension parameter Q = 1
    NATURAL_F0,
    F0_DIV_Q,
    NUM_STATIC_PARAMS
  };

  enum DerivedParamIndex
  {
    RELATIVE_DISP_1,
    RELATIVE_DISP_2,
    ABSOLUTE_DISP_1,
    ABSOLUTE_DISP_2,
    CURRENT_LENGTH,
    CURRENT_THICKNESS_1,
    CURRENT_THICKNESS_2,
    CURRENT_AREA_1,
    CURRENT_AREA_2,
    CURRENT_TENSION,            // The tension parameter Q
    CONTACT_AREA,
    NUM_DERIVED_PARAMS
  };

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  TriangularGlottis();

  // Functions that overwrite the virtual functions in the base class.

  string getName();
  void resetMotion();
  void incTime(const double timeIncrement_s, const double pressure_dPa[]);
  void calcGeometry();
  void getTubeData(double *length_cm, double *area_cm2);
  int getApertureParamIndex();
  
  virtual double getAspirationStrength_dB();

  // Additional functions
  
  double getTensionParameter(double f0);
  void getLengthAndThickness(const double Q, double &length_cm, double thickness[]);
  void getOpenCloseDimensions(double openLength[], double contactLength[], 
    double meanOpenWidth[], double meanContactZ[]);
  double getContactArea(double backX[], double frontX[], double length, double thickness);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  static const int BUFFER_LENGTH = 4;
  static const int BUFFER_MASK = 3;

  double relativeDisplacementBuffer[2][BUFFER_LENGTH];

  IirFilter supraglottalPressureFilter;

  /// Absolute position in samples.
  int pos;
};

#endif
