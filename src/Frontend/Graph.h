// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __GRAPH_H__
#define __GRAPH_H__

#include <wx/wx.h>

// ****************************************************************************
// Internally, all numerical values for the quantities are in CGS units!
// ****************************************************************************

struct PhysicalQuantity
{
  char name[64];
  char symbol[64];
  char mksUnit[64];
  char cgsUnit[64];
  double mksToCgsFactor;
};

enum PhysicalQuantityIndex 
{ 
  PQ_LENGTH,
  PQ_MASS,
  PQ_TIME,
  PQ_AREA,
  PQ_FORCE,
  PQ_MOMENTUM,
  PQ_FREQUENCY,
  PQ_ANGLE,
  PQ_ANGULAR_VELOCITY,
  PQ_PRESSURE,
  PQ_MASSFLOW,
  PQ_DENSITY,
  PQ_VOLUME_VELOCITY,
  PQ_TEMPERATURE,
  PQ_RATIO,
  PQ_VELOCITY,
  NUM_PHYSICAL_QUANTITIES
};

const PhysicalQuantity physicalQuantity[NUM_PHYSICAL_QUANTITIES] = 
{
  "length",           "l",      "m",      "cm",         100.0,
  "mass",             "m",      "kg",     "g",          1000.0,
  "time",             "t",      "s",      "s",          1.0,
  "area",             "A",      "m^2",    "cm^2",       10000.0,
  "force",            "F",      "N",      "dyne",       100000.0,
  "momentum",         "(mv)",   "N-s",    "dyne-s",     100000.0,
  "frequency",        "f",      "Hz",     "Hz",         1.0,
  "angle",            "phi",    "rad",    "rad",        1.0,
  "angular velocity", "omega",  "rad/s",  "rad/s",      1.0,
  "pressure",         "P",      "Pa",     "dPa",        10.0,
  "mass flow",        "dm/dt",  "kg/s",   "g/s",        1000.0,
  "density",          "rho",    "kg/m^3", "g/cm^3",     0.001,
  "volume velocity",  "dV/dt",  "m^3/s",  "cm^3/s",     1000000.0,
  "temperature",      "T",      "K",      "-",          1.0,
  "ratio",            "",       "",       "",           1.0,
  "velocity",         "v",      "m/s",    "cm/s",       100.0
};

// ****************************************************************************
// ****************************************************************************

class Graph
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  struct LinearDomain
  {
    PhysicalQuantityIndex quantity;
    // All in CGS units!
    double reference;
    double scaleDivision;
  
    double negativeLimitMin;
    double negativeLimitMax;
    double negativeLimit;

    double positiveLimitMin;
    double positiveLimitMax;
    double positiveLimit;

    int    numZoomSteps;
    int    postDecimalPositions;        // Number of post-decimal positions
    bool   useCgsUnit;                  // Scale inscriptions in cgs-units?
    bool   useRelativeInscription;      // Inscription relative to reference value
    bool   showGrayLines;
  } 
  abscissa, linearOrdinate;
  
  struct LogDomain
  {
    // All in CGS units!
    double reference;
    double scaleDivision;
    double lowerLevelMin;
    double lowerLevelMax;
    double lowerLevel;

    double upperLevelMin;
    double upperLevelMax;
    double upperLevel;

    bool   showGrayLines;
    double zoomStep;                    // in dB      
  } 
  logOrdinate;

  // Options.

  bool isLinearOrdinate;
  bool abscissaAtBottom;                    // bottom or top
  bool ordinateAtLeftSide;                  // left or right

  // **************************************************************************
  // Public functions.
  // All parameters of the functions are in CGS units!
  // **************************************************************************

public:
  Graph();

  void init(wxWindow *ctrl, int leftMargin, int rightMargin, int topMargin, int bottomMargin);
  void getDimensions(int& x, int& y, int& w, int& h);
  void getMargins(int& left, int& right, int& top, int& bottom);

  void initAbscissa(PhysicalQuantityIndex quantity, double reference, double scaleDivision,
                    double negativeLimitMin, double negativeLimitMax, double negativeLimit,
                    double positiveLimitMin, double positiveLimitMax, double positiveLimit,
                    int numZoomSteps, int postDecimalPositions, 
                    bool useCgsUnit, bool useRelativeInscription, bool showGrayLines);
  
  void initLinearOrdinate(PhysicalQuantityIndex quantity, double reference, double scaleDivision,
                    double negativeLimitMin, double negativeLimitMax, double negativeLimit,
                    double positiveLimitMin, double positiveLimitMax, double positiveLimit,
                    int numZoomSteps, int postDecimalPositions, 
                    bool useCgsUnit, bool useRelativeInscription, bool showGrayLines);

  void initLogOrdinate(double reference, double scaleDivision,
                      double lowerLevelMin, double lowerLevelMax, double lowerLevel,
                      double upperLevelMin, double upperLevelMax, double upperLevel,
                      bool showGrayLines, double zoomStep);

  void paintAbscissa(wxDC &dc);
  void paintOrdinate(wxDC &dc);

  void zoomInAbscissa(bool negativeLimit, bool positiveLimit);
  void zoomOutAbscissa(bool negativeLimit, bool positiveLimit);
  void zoomInOrdinate(bool negativeLimit, bool positiveLimit);
  void zoomOutOrdinate(bool negativeLimit, bool positiveLimit);

  // Transformationen phys. Gr??e <-> Pixelzeile/-spalte **********

  int getXPos(double absXValue);
  int getYPos(double absYValue);
  double getAbsXValue(int xPos);
  double getAbsYValue(int yPos);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  static const double EPSILON;
  
  // Margins of the graph to the edge of the picture
  int leftMargin, rightMargin, bottomMargin, topMargin;
  wxWindow *control;

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  void getZoomFactors(LinearDomain *domain, double& positiveZoomFactor, double& negativeZoomFactor);
};

#endif