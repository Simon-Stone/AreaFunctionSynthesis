// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// The new solver (Cholesky factorization) was added by Johann Marwitz.
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __TDS_MODEL__
#define __TDS_MODEL__

#include <cmath>
#include "Dsp.h"
#include "Geometry.h"
#include "IirFilter.h"
//#include "GlottisTitze.h"
#include "Tube.h"
#include "Constants.h"
#include "TimeFunction.h"


// ****************************************************************************
/// Class for the simulation of vocal tract acoustics in the time domain on the
/// basis of a branched tube model of the vocal tract.
/// All simulation is calculated in CGS-units.
// ****************************************************************************

class TdsModel
{
public:

  // ************************************************************************
  // Constants.
  // ************************************************************************

  /// Two additional output currents at the mouth and nose opening
  static const int NUM_BRANCH_CURRENTS = Tube::NUM_SECTIONS + 4;

  // Max. number of non-zero places per row in the matrix (for Gauss-Seidel)
  static const int MAX_CONCERNED_MATRIX_COLUMNS = 16;
  static const int MAX_CONSTRICTIONS = 4;

  static const int NUM_NOISE_BUFFER_SAMPLES = 8;
  static const int NOISE_BUFFER_MASK = 7;

  // Max. number of non-zero places per row in the matrix when it is saved in symmetric envelope structure (for cholesky factorization)
  static const int MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE = 57;
  // Max. number of non-zero places per	column in the matrix when it is saved in symmetric envelope structure (for cholesky factorization)
  static const int MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE = 10;

  // Max. number of non-zero places per row in the symmetric saved matrix
  static const int MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC = 3;
  // Max. number of non-zero places per	column the symmetric saved matrix
  static const int MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC = 4;

  // For the temporal discretization
  static const double THETA;
  static const double THETA1;

  static const double MIN_AREA_CM2;
  static const double NOISE_CUTOFF_FREQ;


  // ************************************************************************
  /// Some options
  // ************************************************************************

  enum SolverType
  {
    SOR_GAUSS_SEIDEL,
    CHOLESKY_FACTORIZATION,
    NUM_SOLVER_TYPES
  };

  enum GlottisLossOptions
  {
    STANDARD_ENTRANCE_LOSS,
    ENTRANCE_LOSS_VAN_DEN_BERG,
    VARIABLE_ENTRANCE_LOSS,
    NUM_GLOTTIS_LOSS_OPTIONS
  };

  struct Options
  {
    bool turbulenceLosses;        ///< Consider fluid dynamic losses due to turbulence
    bool softWalls;               ///< Consider losses due to soft walls
    bool generateNoiseSources;    ///< Generate noise sources during the simulation
    bool radiationFromSkin;       ///< Allow sound radiation from the skin
    bool piriformFossa;           ///< Include the piriform fossa
    bool innerLengthCorrections;  ///< Additional inductivities between adjacent sections
    bool transvelarCoupling;      ///< Sound transmission through the velum tissue?
    int glottisLossOption;        ///< One of the GlottisLossOptions
    double flowSeparationAreaRatio; ///< At which glottal A2/A1 ratio does flow separation occur?
    SolverType solverType;
  };

  // ************************************************************************
  /// Information for one individual lumped noise source
  // ************************************************************************

  struct NoiseSource
  {
    bool   isFirstOrder;
    double cutoffFrequency;
    double dampingFactor;
    double targetAmp;
    double currentAmp;
    // Input and output sample buffers for the spectral shaping filter
    double inputBuffer[NUM_NOISE_BUFFER_SAMPLES];
    double outputBuffer[NUM_NOISE_BUFFER_SAMPLES];
    double sample;        ///< The current sampling point of the noise source
  };

  // ************************************************************************
  /// Information about a supraglottal constriction
  // ************************************************************************

  struct Constriction
  {
    int firstSection;
    int lastSection;
    int narrowestSection;
    double obstaclePos;
    Tube::Articulator articulator;
    double laterality;
  };


  // ************************************************************************
  /// Structure for one individual branch current in the electrical
  /// network. The identity of a branch current is defined by the
  /// indices of the tube sections from where it comes and where
  /// it goes.
  // ************************************************************************

  struct BranchCurrent
  {
    int sourceSection;
    int targetSection;

    double magnitude;
    double magnitudeRate;
    double noiseMagnitude;    ///< The flow low-pass filtered at 1000 Hz
  };

  // ************************************************************************
  /// An individual short homogeneous tube section.
  // ************************************************************************

  struct TubeSection
  {
    bool isDynamic;        ///< Can the network components R, C, L, ... change ?

    double pos;
    double area;
    double areaRate;
    double length;
    double volume;
    double volumeRate;        ///< Volume change per unit time
    Tube::Articulator articulator;
    double laterality;

    NoiseSource monopoleSource; ///< Is created in the center of the tube section
    NoiseSource dipoleSource;   ///< Is created at the entrance of the tube section

    double pressure;
    double pressureRate;

    /// \name Indices of the inflowing and outflowing currents
    /// @{
    int currentIn;
    int currentOut[2];
    /// @}

    /// \name Wall properties
    /// @{
    double Mw;    ///< Mass per unit-area
    double Bw;    ///< Resistance per unit-area
    double Kw;    ///< Stiffness per unit-area
    /// @}

    /// \name Wall currents
    /// @{
    double wallCurrent;       ///< Current flow "into" the walls
    double wallCurrentRate;   ///< 1st derivative of the wall-flow
    double wallCurrentRate2;  ///< 2nd derivative of the wall-flow
    /// @}

    double L;        ///< Inductivity
    double C;        ///< Capacity
    double R[2];     ///< Ohm's resistance left and right
    double S;        ///< Pressure source at the inlet of the section (A constant in the pressure-difference eq.)

    // For the wall vibration
    double alpha;
    double beta;

    // Temporary values
    double D;
    double E;
  };

  // ************************************************************************
  // Public variables.
  // ************************************************************************

  double flowSourceAmp;
  int    flowSourceSection;   ///< Where is the periodic volume velocity source ?
  double pressureSourceAmp;
  int    pressureSourceSection;  ///< Where is the static pressure source ?

  NoiseSource lipsDipoleSource;   ///< Last noise source at the mouth opening

  TubeSection tubeSection[Tube::NUM_SECTIONS];
  BranchCurrent branchCurrent[NUM_BRANCH_CURRENTS];

  // Help variables to effectively solve the system of eqs. with Gauss-Seidel
  int numFilledRowValues[NUM_BRANCH_CURRENTS];
  int filledRowIndex[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_COLUMNS];

  // Help variables to effectively solve the system of eqs. with cholesky factorization
  int numFilledRowValuesSymmetricEnvelope[NUM_BRANCH_CURRENTS];
  int filledRowIndexSymmetricEnvelope[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE];
  int numFilledColumnValuesSymmetricEnvelope[NUM_BRANCH_CURRENTS];
  int filledColumnIndexSymmetricEnvelope[NUM_BRANCH_CURRENTS][MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE];

  bool doNetworkInitialization;
  double timeStep;
  /// Aspiration strength from -40 dB to 0 dB.
  double aspirationStrength_dB;

  double matrix[NUM_BRANCH_CURRENTS][NUM_BRANCH_CURRENTS];
  double factorizationMatrix[NUM_BRANCH_CURRENTS][NUM_BRANCH_CURRENTS];
  double solutionVector[NUM_BRANCH_CURRENTS];
  double flowVector[NUM_BRANCH_CURRENTS];

  int SORIterations;      ///< For external evaluation of the iterations needed to solve the system of eqs.

  IirFilter glottalToneFilter;
  IirFilter transglottalPressureFilter;
  double glottalBernoulliFactor;

  IirFilter transvelarCouplingFilter1;
  IirFilter transvelarCouplingFilter2;

  Constriction constriction[MAX_CONSTRICTIONS];
  int numConstrictions;

  Options options;


  // ************************************************************************
  // Public methods of the class
  // ************************************************************************

public:
  TdsModel();
  ~TdsModel();

  void initModel();
  void resetMotion();

  // **************************************************************
  /// \name These functions should be called for each time step
  // **************************************************************
  /// @{
  void setTube(Tube *tube, bool filtering = false);
  void getTube(Tube *tube);
  void setFlowSource(double flow_cm3_s, int section);
  void setPressureSource(double pressure_dPa, int section = Tube::FIRST_TRACHEA_SECTION);
  double proceedTimeStep(const string &matrixFileName = "");
  /// @}

  // **************************************************************

  void solveEquationsSor(const string &matrixFileName = "");
  void solveEquationsCholesky();
  int getSampleIndex() { return position; }
  void getSectionFlow(int sectionIndex, double &inflow, double &outflow);
  double getSectionPressure(int sectionIndex);
  void checkGlottalEntranceLossCoeffFlucher2011();

  // ************************************************************************
  // Private data.
  // ************************************************************************

private:
  int position;         ///< Internal counter for the sampling position
  double teethPosition; ///< Position of the teeth (from the glottis)
  
  // ************************************************************************
  // Private functions.
  // ************************************************************************

private:
  void prepareTimeStep();

  double getGlottalEntranceLossCoeffFlucher2011();
  double getGlottalEntranceLossCoeffFlucher2011(double pressure_dPa, double d_cm);

  double getJunctionInductance(double A1_cm2, double A2_cm2);

  void calcMatrix();
  void updateVariables();
  
  void calcNoiseSources();
  void calcNoiseSample(NoiseSource *s, double ampThreshold);

  double getCurrentIn(const int section);
  double getCurrentOut(const int section);
  double getCurrentIn(const TubeSection *ts);
  double getCurrentOut(const TubeSection *ts);
};


#endif
