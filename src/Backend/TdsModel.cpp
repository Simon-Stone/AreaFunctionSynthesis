// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "TdsModel.h"
#include "TridiagonalEquationSet.h"
#include <fstream>
#include <iomanip>
#include <cstdlib>


const double TdsModel::THETA  = 0.515;          // ideal is 0.515
const double TdsModel::THETA1 = 1.0 - TdsModel::THETA;

// For MIN_AREA_CM2, 0.1 mm^2 seems to be a very good value.
const double TdsModel::MIN_AREA_CM2 = 0.1E-2; // = 0.1 mm^2;

// Cutoff-freq. of the low-pass filter for the flow that induces friction
const double TdsModel::NOISE_CUTOFF_FREQ = 500.0;     


// ****************************************************************************
/// Constructor.
// ****************************************************************************

TdsModel::TdsModel()
{
  // ****************************************************************
  // Acoustic options
  // ****************************************************************

  options.turbulenceLosses = true;
  options.softWalls = true;
  options.generateNoiseSources = true;
  options.radiationFromSkin = true;
  options.piriformFossa = false;
  options.innerLengthCorrections = true;
  options.transvelarCoupling = false;
  options.glottisLossOption = STANDARD_ENTRANCE_LOSS;
  options.flowSeparationAreaRatio = 1.0;
  options.solverType = CHOLESKY_FACTORIZATION; // CHOLESKY_FACTORIZATION | SOR_GAUSS_SEIDEL

  // ****************************************************************

  Tube tube;
  setTube(&tube);

  initModel();
  resetMotion();
}

// ****************************************************************************
/// Destructor.
// ****************************************************************************

TdsModel::~TdsModel()
{
}


// ****************************************************************************
// ****************************************************************************

void TdsModel::initModel()
{
  TubeSection *ts;
  BranchCurrent *bc;
  int i, j, k;

  SORIterations = 0;
  timeStep = 1.0 / (double)TDS_SAMPLING_RATE;

  // The network components of the static tube segments must be
  // initialized once at the beginning.
  doNetworkInitialization = true;

  // ****************************************************************
  // Tube section for a (static) pressure or volume velocity source.
  // ****************************************************************

  flowSourceSection      = -1;
  pressureSourceSection  = -1;

  // ****************************************************************
  // Assign the source and target sections for the branch currents.
  // ****************************************************************

  // Every tube section has one in-flowing current

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    bc = &branchCurrent[i];
    bc->sourceSection = i - 1;
    bc->targetSection = i;
  }

  // Branch currents with special source sections.

  branchCurrent[Tube::FIRST_TRACHEA_SECTION].sourceSection = -1;
  branchCurrent[Tube::FIRST_NOSE_SECTION].sourceSection = Tube::LAST_PHARYNX_SECTION;
  branchCurrent[Tube::FIRST_FOSSA_SECTION].sourceSection = 
    Tube::FIRST_PHARYNX_SECTION + Tube::FOSSA_COUPLING_SECTION;
  
  for (i=0; i < Tube::NUM_SINUS_SECTIONS; i++)
  {
    branchCurrent[Tube::FIRST_SINUS_SECTION + i].sourceSection = 
      Tube::FIRST_NOSE_SECTION + Tube::SINUS_COUPLING_SECTION[i];
  }

  // The two currents from the mouth opening in to free space.

  bc = &branchCurrent[Tube::NUM_SECTIONS];
  bc->sourceSection = Tube::LAST_MOUTH_SECTION;
  bc->targetSection = -1;

  bc = &branchCurrent[Tube::NUM_SECTIONS + 1];
  bc->sourceSection = Tube::LAST_MOUTH_SECTION;
  bc->targetSection = -1;

  // The two currents from the nostrils into free space.

  bc = &branchCurrent[Tube::NUM_SECTIONS + 2];
  bc->sourceSection = Tube::LAST_NOSE_SECTION;
  bc->targetSection = -1;

  bc = &branchCurrent[Tube::NUM_SECTIONS + 3];
  bc->sourceSection = Tube::LAST_NOSE_SECTION;
  bc->targetSection = -1;


  // ****************************************************************
  // Assign each tube section the currents that flow into and out of
  // the section.
  // ****************************************************************

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];

    ts->currentIn     = -1;
    ts->currentOut[0] = -1;
    ts->currentOut[1] = -1;

    ts->isDynamic = true;
  }

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    bc = &branchCurrent[i];

    // Where does it come from ?
    if (bc->sourceSection != -1)
    {
      ts = &tubeSection[bc->sourceSection];
      if (ts->currentOut[0] == -1) 
        { ts->currentOut[0] = i; } 
      else 
        { ts->currentOut[1] = i; }
    }

    // Where does it flow to ?
    if (bc->targetSection != -1)
    {
      ts = &tubeSection[bc->targetSection];
      ts->currentIn = i;
    }
  }


  // ****************************************************************
  // Init the matrix and the help structures for the Cholesky 
  // factorization.
  // ****************************************************************

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    numFilledRowValuesSymmetricEnvelope[i] = 0;
    for (k = 0; k < MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE; k++)
    {
      filledRowIndexSymmetricEnvelope[i][k] = 0;
    }
  }

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    numFilledColumnValuesSymmetricEnvelope[i] = 0;
    for (k = 0; k < MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE; k++)
    {
      filledColumnIndexSymmetricEnvelope[i][k] = 0;
    }
  }

  // Fill the whole matrix with zeros.

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    for (j = 0; j < NUM_BRANCH_CURRENTS; j++)
    {
      matrix[i][j] = 0.0;
    }
  }

  // Fill the whole factorizationMatrix with zeros.

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    for (j = 0; j < NUM_BRANCH_CURRENTS; j++)
    {
      factorizationMatrix[i][j] = 0.0;
    }
  }

  // Set a 1 to all non-zero places in the matrix.

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    if (branchCurrent[i].sourceSection != -1)
    {
      ts = &tubeSection[branchCurrent[i].sourceSection];
      if (ts->currentIn != -1) { matrix[i][ts->currentIn] = 1.0; }
      if (ts->currentOut[0] != -1) { matrix[i][ts->currentOut[0]] = 1.0; }
      if (ts->currentOut[1] != -1) { matrix[i][ts->currentOut[1]] = 1.0; }
    }

    if (branchCurrent[i].targetSection != -1)
    {
      ts = &tubeSection[branchCurrent[i].targetSection];
      if (ts->currentIn != -1) { matrix[i][ts->currentIn] = 1.0; }
      if (ts->currentOut[0] != -1) { matrix[i][ts->currentOut[0]] = 1.0; }
      if (ts->currentOut[1] != -1) { matrix[i][ts->currentOut[1]] = 1.0; }
    }

    // Keep in mind the columns elements with a number different from
    // zero in one row in the lower left triangular matrix and the 
    // columns elements between the first non zeros in one row and 
    // the main diagonal (symmetric envelope format).
    // The indices on the main diagonal are not kept in mind explicitly.

    numFilledRowValuesSymmetricEnvelope[i] = 0;
    int betweenNonZeros = 0;

    for (j = 0; j < i; j++)
    {
      if ((matrix[i][j] != 0.0) || betweenNonZeros != 0)
      {
        betweenNonZeros = 1;
        matrix[i][j] = 1;  // for the calculation of the envelpe column indices

        if (numFilledRowValuesSymmetricEnvelope[i] < MAX_CONCERNED_MATRIX_COLUMNS_SYMMETRIC_ENVELOPE)
        {
          filledRowIndexSymmetricEnvelope[i][numFilledRowValuesSymmetricEnvelope[i]] = j;
          numFilledRowValuesSymmetricEnvelope[i]++;
        }
        else
        {
          printf("Error: Attention: The max. number of used rows and columns has been "
            "exceeded in prepareCholeskyFactorization().\n");
        }
      }
    }
  }

  // Keep in mind the rows elements with a number different from zero
  // in one column in the lower left triangular matrix.
  // (symmetric envelope format)
  // The indices on the main diagonal are not kept in mind explicitly.

  for (j = 0; j < NUM_BRANCH_CURRENTS; j++)
  {
    numFilledColumnValuesSymmetricEnvelope[j] = 0;

    for (i = j + 1; i < NUM_BRANCH_CURRENTS; i++)
    {
      if (matrix[i][j] != 0.0)
      {

        if (numFilledColumnValuesSymmetricEnvelope[j] < MAX_CONCERNED_MATRIX_ROWS_SYMMETRIC_ENVELOPE)
        {
          filledColumnIndexSymmetricEnvelope[j][numFilledColumnValuesSymmetricEnvelope[j]] = i;
          numFilledColumnValuesSymmetricEnvelope[j]++;
        }
        else
        {
          printf("Error: Attention: The max. number of used rows and columns has been "
            "exceeded in prepareCholeskyFactorization().\n");
        }
      }
    }
  }

  // ****************************************************************
  // Init the matrix and the help structures for the Gauss-Seidel 
  // method.
  // ****************************************************************

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    numFilledRowValues[i] = 0;
    for (k=0; k < MAX_CONCERNED_MATRIX_COLUMNS; k++)
    {
      filledRowIndex[i][k] = 0;
    }
  }

  // Fill the whole matrix with zeros.

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    for (j=0; j < NUM_BRANCH_CURRENTS; j++) 
    { 
      matrix[i][j] = 0.0; 
    }
  }

  // Set a 1 to all non-zero places in the matrix.

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    if (branchCurrent[i].sourceSection != -1)
    {
      ts = &tubeSection[ branchCurrent[i].sourceSection ];
      if (ts->currentIn     != -1) { matrix[i][ ts->currentIn ] = 1.0; }
      if (ts->currentOut[0] != -1) { matrix[i][ ts->currentOut[0] ] = 1.0; }
      if (ts->currentOut[1] != -1) { matrix[i][ ts->currentOut[1] ] = 1.0; }
    }

    if (branchCurrent[i].targetSection != -1)
    {
      ts = &tubeSection[ branchCurrent[i].targetSection ];
      if (ts->currentIn     != -1) { matrix[i][ ts->currentIn ] = 1.0; }
      if (ts->currentOut[0] != -1) { matrix[i][ ts->currentOut[0] ] = 1.0; }
      if (ts->currentOut[1] != -1) { matrix[i][ ts->currentOut[1] ] = 1.0; }
    }

    // Keep in mind the columns with a number different from zero. The indices
    // on the main diagonal are not kept in mind explicitly.
    
    numFilledRowValues[i] = 0;
    
    for (j=0; j < NUM_BRANCH_CURRENTS; j++)
    {
      if ((matrix[i][j] != 0.0) && (i != j))
      {
        if (numFilledRowValues[i] < MAX_CONCERNED_MATRIX_COLUMNS)
        {
          filledRowIndex[i][ numFilledRowValues[i] ] = j;
          numFilledRowValues[i]++;
        }
        else
        {
          printf("Error: Attention: The max. number of used rows and columns has been "
            "exceeded in prepareGaussSeidel().\n");
        }
      }
    }
  }

}


// ****************************************************************************
/// Sets all flows, pressures, etc. into the equilibrium state.
// ****************************************************************************

void TdsModel::resetMotion()
{
  const double Q = 0.5;
  TubeSection *ts;
  BranchCurrent *bc;
  int i, k;

  // ****************************************************************
  // The tube sections.
  // ****************************************************************

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];

    ts->areaRate         = 0.0;

    ts->pressure         = 0.0;
    ts->pressureRate     = 0.0;
    ts->wallCurrent      = 0.0;
    ts->wallCurrentRate  = 0.0;
    ts->wallCurrentRate2 = 0.0;
    
    // Intermediate values ********************************
    ts->L                = 0.0;
    ts->C                = 0.0;
    ts->R[0] = ts->R[1]  = 0.0;
    ts->S                = 0.0;
    ts->alpha            = 0.0;
    ts->beta             = 0.0;
    ts->D                = 0.0;
    ts->E                = 0.0;

    // Set the noise sources to 0 ***********************************

    ts->monopoleSource.targetAmp = 0.0;
    ts->monopoleSource.currentAmp = 0.0;
    ts->monopoleSource.isFirstOrder = true;
    ts->monopoleSource.cutoffFrequency = 3000.0;
    ts->monopoleSource.dampingFactor = 0.5;
    ts->monopoleSource.sample = 0.0;
    for (k = 0; k < NUM_NOISE_BUFFER_SAMPLES; k++)
    {
      ts->monopoleSource.inputBuffer[k] = 0.0;
      ts->monopoleSource.outputBuffer[k] = 0.0;
    }

    ts->dipoleSource.targetAmp = 0.0;
    ts->dipoleSource.currentAmp = 0.0;
    ts->dipoleSource.isFirstOrder = true;
    ts->dipoleSource.cutoffFrequency = 3000.0;
    ts->dipoleSource.dampingFactor = 0.5;
    ts->dipoleSource.sample = 0.0;
    for (k = 0; k < NUM_NOISE_BUFFER_SAMPLES; k++)
    {
      ts->dipoleSource.inputBuffer[k] = 0.0;
      ts->dipoleSource.outputBuffer[k] = 0.0;
    }
  }

  // Extra dipole source at the lips ********************************

  lipsDipoleSource.targetAmp = 0.0;
  lipsDipoleSource.currentAmp = 0.0;
  lipsDipoleSource.isFirstOrder = true;
  lipsDipoleSource.cutoffFrequency = 3000.0;
  lipsDipoleSource.dampingFactor = 0.5;
  lipsDipoleSource.sample = 0.0;
  for (k = 0; k < NUM_NOISE_BUFFER_SAMPLES; k++)
  {
    lipsDipoleSource.inputBuffer[k] = 0.0;
    lipsDipoleSource.outputBuffer[k] = 0.0;
  }

  // ****************************************************************
  // The branch currents.
  // ****************************************************************

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    bc = &branchCurrent[i];

    bc->magnitude = 0.0;
    bc->magnitudeRate = 0.0;
    bc->noiseMagnitude = 0.0;
  }

  flowSourceAmp = 0.0;
  pressureSourceAmp = 0.0;
  position = 0;             // The internal position counter

  numConstrictions = 0;

  aspirationStrength_dB = Tube::DEFAULT_ASPIRATION_STRENGTH_DB;

  doNetworkInitialization = true;

  glottalBernoulliFactor = 0.0;

  // The volume velocity vector is 0. *******************************
  
  for (i=0; i < NUM_BRANCH_CURRENTS; i++) 
  { 
    flowVector[i] = 0.0; 
  }

  transglottalPressureFilter.createChebyshev(50.0 / (double)TDS_SAMPLING_RATE, false, 4);
  transglottalPressureFilter.resetBuffers();


  // ****************************************************************
  // Reset the two filters H1(s) and H2(s) for transvelar coupling.
  // They work to create a volume velocity U_nose in the nasal cavity from the
  // difference of the pressures P_below and P_above below and above the velum.
  // The filter coefficients were derived from a matched z-transform of the
  // analog circuit model for the velum according to Dang et al. (2016, JASA)
  // for a sampling rate of 44100 Hz.
  // U_nose(s) = P_below(s)*H1(s) + P_above(s)*H2(s).
  // ****************************************************************

  const int COUPLING_FILTER_ORDER = 4;

  // The high precision of the coefficients is strictly necessary.
  // A precision of only 4 post decimal positions would make the filter "explode".
  
  // Filter for cgs units !!!!
  const double a1[COUPLING_FILTER_ORDER + 1] =
  {
    5.027640021717718e-007,
    -7.995535578908732e-007,
    2.967895557191014e-007,
    0.0,
    0.0
  };

  const double b1[COUPLING_FILTER_ORDER + 1] =
  {
    0.0,
    3.986308869708467,
    -5.959669638387298,
    3.960408461107104,
    -0.987047716233603
  };

  const double a2[COUPLING_FILTER_ORDER + 1] =
  {
    6.589309727087047e-004,
    -0.001972281980771,
    0.001968000742164,
    -6.546497341015677e-004,
    0.0
  };

  const double b2[COUPLING_FILTER_ORDER + 1] =
  {
    0.0,
    3.986308869708467,
    -5.959669638387298,
    3.960408461107104,
    -0.987047716233603
  };
  
  transvelarCouplingFilter1.setCoefficients(a1, b1, COUPLING_FILTER_ORDER);
  transvelarCouplingFilter2.setCoefficients(a2, b2, COUPLING_FILTER_ORDER);
  transvelarCouplingFilter1.resetBuffers();
  transvelarCouplingFilter2.resetBuffers();

  // ****************************************************************
  // Reset the filter to transform the pressure just above the 
  // vocal folds to a volume velocity of the outer larynx skin for
  // the glottal tone (voice bar).
  // Use the same filter here as for the transvelar coupling.
  // ****************************************************************

//  glottalToneFilter.createSecondOrderLowpass(200.0 / (double)TDS_SAMPLING_RATE, 1 / sqrt(2.0));

  glottalToneFilter.setCoefficients(a1, b1, COUPLING_FILTER_ORDER);
  glottalToneFilter.resetBuffers();
}


// ****************************************************************************
// ****************************************************************************

void TdsModel::setTube(Tube *tube, bool filtering)
{
  int i;
  TubeSection *target = NULL;
  Tube::Section *source = NULL;

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    source = tube->section[i];
    target = &tubeSection[i];

    target->pos    = source->pos_cm;
    target->area   = source->area_cm2;
    target->length = source->length_cm;
    target->volume = source->volume_cm3;

    target->Mw     = source->wallMass_cgs;
    target->Bw     = source->wallResistance_cgs;
    target->Kw     = source->wallStiffness_cgs;

    target->volumeRate = 0.0;
    target->areaRate   = 0.0;

    target->articulator = source->articulator;
    target->laterality = source->laterality;
  }

  // Set the teeth position
  teethPosition = tube->teethPosition_cm;

  // Set the aspiration strength
  aspirationStrength_dB = tube->aspirationStrength_dB;
}


// ****************************************************************************
// ****************************************************************************

void TdsModel::getTube(Tube *tube)
{
  int i;
  TubeSection *source = NULL;
  Tube::Section *target = NULL;

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    source = &tubeSection[i];
    target = tube->section[i];

    target->pos_cm     = source->pos;
    target->area_cm2   = source->area;
    target->length_cm  = source->length;
    target->volume_cm3 = source->volume;

    target->wallMass_cgs       = source->Mw;
    target->wallResistance_cgs = source->Bw;
    target->wallStiffness_cgs  = source->Kw;

    target->articulator = source->articulator;
    target->laterality = source->laterality;
  }

  // Set the teeth position
  tube->teethPosition_cm = teethPosition;

  // Set the aspiration strength
  tube->aspirationStrength_dB = aspirationStrength_dB;
}


// ****************************************************************************
/// Sets the flow source into the given tube section and to the given strength.
/// \param flow_cm3_s Strength of the source in cm3/s.
/// \param section Section of the source. Set section to -1 to disable the 
/// flow source.
// ****************************************************************************

void TdsModel::setFlowSource(double flow_cm3_s, int section)
{
  flowSourceSection = section;
  flowSourceAmp     = flow_cm3_s;
}


// ****************************************************************************
/// Sets the pressure source into the given tube section and to the given 
/// strength.
/// \param pressure_dPa Strength of the source in deci-Pascal.
/// \param section Section of the source. Set section to -1 to disable the 
/// pressure source.
// ****************************************************************************

void TdsModel::setPressureSource(double pressure_dPa, int section)
{
  pressureSourceSection = section;
  pressureSourceAmp = pressure_dPa;
}


// ****************************************************************************
/// Calculate a new time step in the simulation.
/// Returns the sum flow radiated from the mouth, the nostrils and the vocal
/// tract walls.
/// \param matrixFileName The name of a text file the coefficient matrix shall
/// be written to.
// ****************************************************************************

double TdsModel::proceedTimeStep(const string &matrixFileName)
{
  TubeSection *ts = NULL;

  // Calculation of resistors and other values.
  prepareTimeStep();
  
  // Calculation of the matrix coefficients.
  calcMatrix();

  if (options.solverType == CHOLESKY_FACTORIZATION)
  {
    // Solve the system of eqs. with cholesky factorization
    solveEquationsCholesky();
  }
  else
  {
    // Default: Solve the system of eqs. with the SOR method.
    solveEquationsSor(matrixFileName);
  }

  // Recalculate all currents, pressures and their derivatives.
  updateVariables();

  // ****************************************************************
  // Get the radiated flow.
  // ****************************************************************

  double radiatedFlow = 0.0;

  ts = &tubeSection[Tube::LAST_MOUTH_SECTION];
  if (ts->currentOut[0] != -1) { radiatedFlow+= branchCurrent[ ts->currentOut[0] ].magnitude; }
  if (ts->currentOut[1] != -1) { radiatedFlow+= branchCurrent[ ts->currentOut[1] ].magnitude; }
  
  ts = &tubeSection[Tube::LAST_NOSE_SECTION];
  if (ts->currentOut[0] != -1) { radiatedFlow+= branchCurrent[ ts->currentOut[0] ].magnitude; }
  if (ts->currentOut[1] != -1) { radiatedFlow+= branchCurrent[ ts->currentOut[1] ].magnitude; }

  // ****************************************************************
  // Consider the sound radiation from the skin near the glottis.
  // ****************************************************************

  if (options.radiationFromSkin)
  {
    ts = &tubeSection[Tube::FIRST_PHARYNX_SECTION];
    radiatedFlow += glottalToneFilter.getOutputSample(ts->pressure);
  }
  
  // Increase the internal position counter to the next sample.
  position++;     

  return radiatedFlow;
}


// ****************************************************************************
/// Calculate the network components.
// ****************************************************************************

void TdsModel::prepareTimeStep()
{
  TubeSection *ts = NULL;
  int i;
  double d;
  double Lw, Rw, Cw;
  double surface;
  double a, b;
  double u;
  double circ;

  // Calculate the values D and E and the components L, R, W, C for
  // all tube sections.

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];
    if (ts->area < MIN_AREA_CM2) 
    { 
      ts->area = MIN_AREA_CM2; 
    }

    ts->S = 0.0;      // No pressure source at the inlet of the section
    circ = 2.0*sqrt(ts->area*M_PI);

    // **************************************************************
    // Recalculate the components of the dynamic tube sections.
    // **************************************************************

    if ((ts->isDynamic) || (doNetworkInitialization))
    {
      // The Helmholtz-Resonators are special.

      if ((i >= Tube::FIRST_SINUS_SECTION) && (i <= Tube::LAST_SINUS_SECTION))
      {
        ts->L    = AMBIENT_DENSITY_CGS*(ts->length / ts->area);
        ts->C    = ts->volume / (AMBIENT_DENSITY_CGS*SOUND_VELOCITY_CGS*SOUND_VELOCITY_CGS);
        ts->R[0] = (8.0*AIR_VISCOSITY_CGS * M_PI * ts->length) / (ts->area * ts->area);
        ts->R[1] = ts->R[0];
      }
      else

      // Normal tube segments.
      {
        // Assume a circular cross-section.
        a = b = sqrt(ts->area / M_PI);      // Half-axes of the circle

        // If the radius of the circle becomes smaller than a threshold,
        // then make the cross-section elliptical with the threshold
        // being the length of the bigger half-axis.
        // This will strongly increase the viscous resistance, esp.
        // just before supraglottal or glottal closures and make it
        // roughly equivalent to that of a rectangular slit.

        double MIN_RADIUS_CM = 0.8;    // = 8 mm (corresponds to area of 2.0 cm^2)

        // For non-glottal sections we need a bigger minimum radius in order
        // to increase the damping and avoid clicks during the creation
        // of closures like in "at" or "et".

        if ((i != Tube::LOWER_GLOTTIS_SECTION) && (i != Tube::UPPER_GLOTTIS_SECTION))
        {
          // More than 1.6 can make some vowels like /y/ sound muffled.
          MIN_RADIUS_CM = 1.6;
        }

        if (a < MIN_RADIUS_CM)
        {
          a = MIN_RADIUS_CM;
          b = ts->area / (M_PI*a);
        }

        ts->L    = (AMBIENT_DENSITY_CGS * 0.5 * ts->length) / ts->area;
        ts->C    = ts->volume / (AMBIENT_DENSITY_CGS * SOUND_VELOCITY_CGS * SOUND_VELOCITY_CGS);
        ts->R[0] = ((2.0 * AIR_VISCOSITY_CGS * ts->length) * (a * a + b * b)) / (M_PI*a*a*a*b*b*b);
        ts->R[1] = ts->R[0];
      }
    }

    // **************************************************************
    // The alpha and beta values for the incorporation of wall 
    // vibration.
    // **************************************************************

    ts->alpha = 0.0;
    ts->beta = 0.0;

  	if ((options.softWalls) && (i != Tube::LOWER_GLOTTIS_SECTION) && (i != Tube::UPPER_GLOTTIS_SECTION))
  	{
      // What is the area of the wall ?
      if ((i >= Tube::FIRST_SINUS_SECTION) && (i <= Tube::LAST_SINUS_SECTION))
      {
        surface = 4.0*M_PI*pow((3.0*ts->volume)/(4.0*M_PI), 2.0/3.0); // Surface of a sphere
      }
      else
      {
        surface = circ*ts->length;    // Surface of a cylinder barrel
      }

      if (surface < MIN_AREA_CM2)
      {
        surface = MIN_AREA_CM2;
      }

      Rw = ts->Bw / surface;
      Lw = ts->Mw / surface;
      Cw = surface / ts->Kw;

      ts->alpha = 1.0 / (Lw / (timeStep*timeStep*THETA*THETA) + Rw / (timeStep*THETA) + 1.0/Cw);
      ts->beta = ts->alpha*(
        ts->wallCurrent*(Lw/(timeStep*timeStep*THETA*THETA) + Rw/(timeStep*THETA)) +
        ts->wallCurrentRate*(Lw*(THETA1/THETA + 1.0)/(timeStep*THETA) + Rw*(THETA1/THETA)) +
        ts->wallCurrentRate2*Lw*(THETA1/THETA)
        );
  	}

  }


  // ****************************************************************
  // End of the tube section loop
  // ****************************************************************

  doNetworkInitialization = false;  

  // ****************************************************************
  // If a wide tube section follows a narrow tube section,
  // the complete loss of the kinetic pressure is assumed.
  // From a wide to a narrow section, the Bernoulli equation is
  // applied.
  // ****************************************************************

  if (options.turbulenceLosses)
  {
    for (i = Tube::FIRST_PHARYNX_SECTION + 1; i <= Tube::LAST_MOUTH_SECTION; i++)
    {
      ts = &tubeSection[i - 1];

      // There must be a simple current from section i-1 to section i 
      // (no branching off).
      // This is to avoid click artifacts due to the coupling of the 
      // sinus pririformis and the nasal cavity.

      if ((ts->currentOut[0] != -1) && (ts->currentOut[1] == -1))
      {
        u = getCurrentOut(ts);

        // Add Bernoulli resistance when the flow is from a wide into
        // a narrow section.
        if (((tubeSection[i].area < tubeSection[i - 1].area) && (u > 0)) ||
          ((tubeSection[i].area > tubeSection[i - 1].area) && (u < 0)))
        {
          ts = &tubeSection[i - 1];
          ts->R[1] -= u*0.5*AMBIENT_DENSITY_CGS / (ts->area*ts->area);
          ts = &tubeSection[i];
          ts->R[0] += u*0.5*AMBIENT_DENSITY_CGS / (ts->area*ts->area);
        }
      }
    }
  }

  // ****************************************************************
  // Possibly decouple the sinus piriformis.
  // ****************************************************************

  if (options.piriformFossa == false)
  {
    // Set entrance of piriform fossae to the (very high) flow resistance
    // of the smallest possible area.
    ts = &tubeSection[Tube::FIRST_FOSSA_SECTION];
    ts->R[0] = 8.0 * AIR_VISCOSITY_CGS * ts->length * M_PI / (MIN_AREA_CM2 * MIN_AREA_CM2);
  }
  
  // ****************************************************************
  // Init the nonlinear resistances at the glottis.
  // ****************************************************************

  double sourceArea;
  double targetArea;
  
  // Entrance loss coefficient of the glottis: 1.0 corresponds to 
  // perfect Benoulli flow.
  double k_ent = 1.0;

  if (options.glottisLossOption == ENTRANCE_LOSS_VAN_DEN_BERG)
  {
    k_ent = 1.375;
  }
  else
  if (options.glottisLossOption == VARIABLE_ENTRANCE_LOSS)
  {
    k_ent = getGlottalEntranceLossCoeffFlucher2011();
  }

  sourceArea = tubeSection[Tube::LAST_TRACHEA_SECTION].area;
  targetArea = tubeSection[Tube::LOWER_GLOTTIS_SECTION].area;
  u = getCurrentIn(Tube::LOWER_GLOTTIS_SECTION);
  
  if (u > 0)
  {
    tubeSection[Tube::LOWER_GLOTTIS_SECTION].R[0] +=
      k_ent * 0.5 * AMBIENT_DENSITY_CGS * fabs(u) *
      (1.0 / (targetArea*targetArea) - 1.0 / (sourceArea*sourceArea));
  }


  // ****************************************************************
  // Transition between the glottal sections:
  // Assume Bernoulli flow, as long as 
  // A_upper < flowSeparationAreaRatio * A_lower 
  // (flowSeparationAreaRatio = 1.0 ... 1.3).
  // ****************************************************************

  sourceArea = tubeSection[Tube::LOWER_GLOTTIS_SECTION].area;
  targetArea = tubeSection[Tube::UPPER_GLOTTIS_SECTION].area;

  double bernoulliTargetFactor = 0.0;
  if (targetArea < options.flowSeparationAreaRatio * sourceArea)
  {
    bernoulliTargetFactor = 1.0;
  }

  const double FILTER_COEFF = 0.8;  // A factor of 0.8 corresponds to a time constant of 0.1 ms
  glottalBernoulliFactor = FILTER_COEFF * glottalBernoulliFactor + (1.0 - FILTER_COEFF) * bernoulliTargetFactor;

  u = getCurrentOut(Tube::LOWER_GLOTTIS_SECTION);

  if (u > 0)
  {
    tubeSection[Tube::LOWER_GLOTTIS_SECTION].R[1] +=
      glottalBernoulliFactor * fabs(u) * 0.5 * AMBIENT_DENSITY_CGS * 
      (1.0 / (targetArea*targetArea) - 1.0 / (sourceArea*sourceArea));
  }


  // At the glottal exit -> total loss of kinetic pressure.
  // Nothing to do here.

  // ****************************************************************
  // Calculate the position and strength of the noise sources.
  // ****************************************************************
 
  if (options.generateNoiseSources) 
  { 
    calcNoiseSources(); 
  }

  // ****************************************************************
  // Calc. the extra flow into the nasal cavity due to transvelar
  // coupling.
  // ****************************************************************

  double transvelarCouplingFlow = 0.0;

  if (options.transvelarCoupling)
  {
    // The extra flow is calculated from the filtered pressures 
    // below and above the velum.
    double p1 = tubeSection[Tube::FIRST_MOUTH_SECTION + 2].pressure;
    double p2 = tubeSection[Tube::FIRST_NOSE_SECTION + 2].pressure;

    transvelarCouplingFlow = transvelarCouplingFilter1.getOutputSample(p1) + transvelarCouplingFilter2.getOutputSample(p2);
  }

  // ****************************************************************
  // Calculate D and E for each tube section.
  // ****************************************************************

  double sourceAmp;

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];

    sourceAmp = ts->monopoleSource.sample;

    if (i == Tube::FIRST_NOSE_SECTION + 2) 
    { 
      sourceAmp += transvelarCouplingFlow; 
    }

    if (i == flowSourceSection) 
    { 
      sourceAmp += flowSourceAmp; 
    }

    d = timeStep*THETA / (ts->C + ts->alpha);
    ts->E = d;
    ts->D = ts->pressure + timeStep*THETA1*ts->pressureRate - 
            d*(ts->beta - sourceAmp);
  }

}


// ****************************************************************************
/// Calculates the glottis entrance loss coefficient according to Flucher
/// (JASA, 2011) based on the minimal diameter of the glottis and P_sub
/// according to their approximation formula in ref. 16.
// ****************************************************************************

double TdsModel::getGlottalEntranceLossCoeffFlucher2011()
{
  if (pressureSourceSection != Tube::FIRST_TRACHEA_SECTION)
  {
    return 1.0;
  }

  double transglottalPressure_dPa = 
    getSectionPressure(Tube::FIRST_TRACHEA_SECTION + Tube::NUM_TRACHEA_SECTIONS - 1) - 
    getSectionPressure(Tube::FIRST_PHARYNX_SECTION);

  // Filter with a Cheby. low-pass filter with 50 Hz cutoff-freq.
  double filteredTransglottalPressure_dPa = transglottalPressureFilter.getOutputSample(transglottalPressure_dPa);

  double GLOTTIS_LENGTH_CM = 1.25;
  double d_cm = tubeSection[Tube::LOWER_GLOTTIS_SECTION].area / GLOTTIS_LENGTH_CM;

  double k_ent = getGlottalEntranceLossCoeffFlucher2011(filteredTransglottalPressure_dPa, d_cm);

  return k_ent;
}


// ****************************************************************************
/// Calculates the glottis entrance loss coefficient according to Flucher
/// (JASA, 2011) based on the given minimal diameter d_cm (in centimeters) and the
/// given transglottal pressure pressure_dPa (in deci-pascal).
// ****************************************************************************

double TdsModel::getGlottalEntranceLossCoeffFlucher2011(double pressure_dPa, double d_cm)
{
  static const double MIN_D_CM = 0.001;  // For numeric reasons
  static const double MIN_PRESSURE_CMH2O = 0.001;  // For numeric reasons

  // Basis is the approximation formula in ref. 16. of Fulcher et al. (2011).
  double pressure_cmH2O = pressure_dPa / 979.7;
  double D = d_cm;     // D is the diameter in cm.

  if (pressure_cmH2O < MIN_PRESSURE_CMH2O)
  {
    pressure_cmH2O = MIN_PRESSURE_CMH2O;
  }

  if (D < MIN_D_CM)
  {
    D = MIN_D_CM;
  }
  
  double logD = log10(D);
  double c1 = 0.7953;
  double c2 = 1.4741;
  double c3 = 0.6529;
  double a = pow(10, c1*logD*logD + c2*logD + c3);

  double d1 = -0.7427;
  double d2 = -1.6209;
  double d3 = -0.875;
  double b = d1*logD*logD + d2*logD + d3;

  double k_ent = a * pow(pressure_cmH2O, b);

  // Safety check!
  if (k_ent < 0.6)
  {
    k_ent = 0.6;
  }

  if (k_ent > 18.0)
  {
    k_ent = 18.0;
  }

  return k_ent;
}


// ****************************************************************************
/// Test function for the calculation and output (console) of the values in the 
/// result table in the Fulcher paper.
// ****************************************************************************

void TdsModel::checkGlottalEntranceLossCoeffFlucher2011()
{
  // Reproduce Table I from Fulcher et al. (2011).
  const double CM_H2O_TO_DPA = 97.97;   // To deci-Pascal
  double d_cm = 0.0;

  printf("Values for k_ent:\n");

  d_cm = 0.005;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.0075;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.01;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.02;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.04;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.08;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.16;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  d_cm = 0.32;
  printf("d=%f cm: %4.3f  %4.3f  %4.3f  %4.3f  %4.3f\n",
    d_cm,
    getGlottalEntranceLossCoeffFlucher2011(3 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(5 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(10 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(15 * CM_H2O_TO_DPA, d_cm),
    getGlottalEntranceLossCoeffFlucher2011(25 * CM_H2O_TO_DPA, d_cm));

  printf("\n");
}


// ****************************************************************************
/// Calculate the positions and amplitudes of the noise sources.
// ****************************************************************************

void TdsModel::calcNoiseSources()
{
  const double MAX_CONSTRICTION_AREA = 1.0;   // 1.0 cm^2
  const double MAX_DELTA_AREA = 0.2;          // 0.2 cm^2
  const double MAX_TEETH_DISTANCE = 2.0;      // = 2 cm

  Constriction *cons = NULL;
  int i, k;
  TubeSection *ts = NULL;
  NoiseSource *s = NULL;

  // ****************************************************************
  // Reset the target amplitude of all noise sources.
  // ****************************************************************

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    tubeSection[i].monopoleSource.targetAmp = 0.0;
    tubeSection[i].dipoleSource.targetAmp = 0.0;
  }
  lipsDipoleSource.targetAmp = 0.0;

  // ****************************************************************
  // The noise source right above the glottis is always there.
  // ****************************************************************

  numConstrictions = 1;
  cons = &constriction[0];

  cons->firstSection = Tube::LOWER_GLOTTIS_SECTION;
  cons->lastSection = Tube::UPPER_GLOTTIS_SECTION;
  // Always take the upper glottis section as the "narrowest" section
  // to avoid flipping of the two sections all the time.
  cons->narrowestSection = Tube::UPPER_GLOTTIS_SECTION;
  cons->obstaclePos = 1.5;   // 1.5 cm above the glottis.
  cons->articulator = Tube::VOCAL_FOLDS;
  cons->laterality = 0.0;


  // ****************************************************************
  // Determine the most constricted parts made with the tongue.
  // ****************************************************************

  double minTongueAreaForTeethSource = 1000000.0;   // Extremely high area.
  double minTongueArea = 1000000.0;   // Extremely high area.
  int minTongueSection = -1;

  for (i = Tube::FIRST_PHARYNX_SECTION; i <= Tube::LAST_MOUTH_SECTION; i++)
  {
    if ((tubeSection[i].articulator == Tube::TONGUE) && 
        (tubeSection[i].area < minTongueArea))
    {
      minTongueArea = tubeSection[i].area;
      minTongueSection = i;
    }
  }

  if (minTongueArea < MAX_CONSTRICTION_AREA)
  {
    cons = &constriction[numConstrictions];
    numConstrictions++;

    cons->articulator = Tube::TONGUE;
    cons->laterality = 0.0;
    cons->narrowestSection = minTongueSection;
    cons->firstSection = minTongueSection;
    cons->lastSection = minTongueSection;

    // Find the first and last section of this constricted region.
    
    double maxArea = minTongueArea + MAX_DELTA_AREA; 

    while ((tubeSection[ cons->firstSection ].area < maxArea) &&
           (tubeSection[ cons->firstSection ].articulator == Tube::TONGUE) &&
           (cons->firstSection > Tube::FIRST_PHARYNX_SECTION)) { cons->firstSection--; }
  
    while ((tubeSection[ cons->lastSection ].area < maxArea) &&
           (tubeSection[ cons->lastSection ].articulator == Tube::TONGUE) &&
           (cons->lastSection < Tube::LAST_MOUTH_SECTION)) { cons->lastSection++; }

    cons->firstSection++;
    cons->lastSection--;

    // Determine the laterality for this constriction (maximum within the channel).

    for (i = cons->firstSection; i <= cons->lastSection; i++)
    {
      if (tubeSection[i].laterality > cons->laterality)
      {
        cons->laterality = tubeSection[i].laterality;
      }
    }

    // Position and distance to an obstacle.

    double jetPos = tubeSection[ cons->lastSection ].pos + tubeSection[ cons->lastSection ].length;

    // The case for /s,S/
    if (teethPosition - jetPos < MAX_TEETH_DISTANCE)
    {
      cons->obstaclePos = teethPosition;
      minTongueAreaForTeethSource = tubeSection[ cons->narrowestSection ].area;
    }
    else
    // The case for /x,ch/
    {
      // The obstacle is in the middle of the section immediately
      // downstream from the jet section. This way, we get 2 pressure
      // sources at each end of the section to approximate the
      // source distribution.
      cons->obstaclePos = tubeSection[ cons->lastSection + 1 ].pos + 
        0.5*tubeSection[ cons->lastSection + 1 ].length;
    }
  }


  // ****************************************************************
  // Is there potentially a second constriction in the tongue region?
  // ****************************************************************

  if ((numConstrictions > 0) && (constriction[numConstrictions-1].articulator == Tube::TONGUE))
  {
    Constriction *prevCons = &constriction[numConstrictions-1];
    minTongueArea = 1000000.0;   // Extremely high area.
    minTongueSection = -1;

    for (i = Tube::FIRST_PHARYNX_SECTION; i <= Tube::LAST_MOUTH_SECTION; i++)
    {
      if ((tubeSection[i].articulator == Tube::TONGUE) && 
          (tubeSection[i].area < minTongueArea) && 
          ((i < prevCons->firstSection) || (i > prevCons->lastSection)))
      {
        minTongueArea = tubeSection[i].area;
        minTongueSection = i;
      }
    }

    if (minTongueArea < MAX_CONSTRICTION_AREA)
    {
      cons = &constriction[numConstrictions];
      numConstrictions++;

      cons->articulator = Tube::TONGUE;
      cons->laterality = 0.0;
      cons->narrowestSection = minTongueSection;
      cons->firstSection = minTongueSection;
      cons->lastSection = minTongueSection;

      // Find the first and last section of this constricted region.
    
      double maxArea = minTongueArea + MAX_DELTA_AREA; 

      while ((tubeSection[ cons->firstSection ].area < maxArea) &&
             (tubeSection[ cons->firstSection ].articulator == Tube::TONGUE) &&
             (cons->firstSection > Tube::FIRST_PHARYNX_SECTION)) { cons->firstSection--; }
  
      while ((tubeSection[ cons->lastSection ].area < maxArea) &&
             (tubeSection[ cons->lastSection ].articulator == Tube::TONGUE) &&
             (cons->lastSection < Tube::LAST_MOUTH_SECTION)) { cons->lastSection++; }

      cons->firstSection++;
      cons->lastSection--;

      if ((cons->firstSection > prevCons->lastSection + 1) || (cons->lastSection < prevCons->firstSection - 1))
      {
        // Determine the laterality for this constriction (maximum within the channel).

        for (i = cons->firstSection; i <= cons->lastSection; i++)
        {
          if (tubeSection[i].laterality > cons->laterality)
          {
            cons->laterality = tubeSection[i].laterality;
          }
        }

        // Position and distance to an obstacle.

        double jetPos = tubeSection[ cons->lastSection ].pos + tubeSection[ cons->lastSection ].length;

        // The case for /s,S/
        if (teethPosition - jetPos < MAX_TEETH_DISTANCE)
        {
          cons->obstaclePos = teethPosition;
          minTongueAreaForTeethSource = tubeSection[ cons->narrowestSection ].area;
        }
        else
        // The case for /x,ch/
        {
          // The obstacle is in the middle of the section immediately
          // downstream from the jet section. This way, we get 2 pressure
          // sources at each end of the section to approximate the
          // source distribution.
          cons->obstaclePos = tubeSection[ cons->lastSection + 1 ].pos + 
            0.5*tubeSection[ cons->lastSection + 1 ].length;
        }
      }
      else
      {
        // Remove this constriction again.
        numConstrictions--;
      }

    }
  }


  // ****************************************************************
  // Determine the most constricted part in the region of the lower 
  // lip.
  // ****************************************************************

  double maxArea = 0.0;
  double minLipArea = 1000000.0;    // Very big area.
  int minLipSection = -1;

  for (i = Tube::FIRST_PHARYNX_SECTION; i <= Tube::LAST_MOUTH_SECTION; i++)
  {
    if ((tubeSection[i].articulator == Tube::LOWER_LIP) && 
        (tubeSection[i].area < minLipArea))
    {
      minLipArea = tubeSection[i].area;
      minLipSection = i;
    }
  }

  // The constriction area has to be smaller than a potential tongue
  // constriction with a source at the upper incisors.

  if ((minLipArea < MAX_CONSTRICTION_AREA) && (minLipArea < minTongueAreaForTeethSource))
  {
    cons = &constriction[numConstrictions];
    numConstrictions++;

    cons->articulator = Tube::LOWER_LIP;
    cons->laterality = 0.0;
    cons->narrowestSection = minLipSection;
    cons->firstSection = minLipSection;
    cons->lastSection = minLipSection;

    // Find the first and last section of this constricted region.
    
    maxArea = minLipArea + MAX_DELTA_AREA; 

    while ((tubeSection[ cons->firstSection ].area < maxArea) &&
           (tubeSection[ cons->firstSection ].articulator == Tube::LOWER_LIP) &&
           (cons->firstSection > Tube::FIRST_PHARYNX_SECTION)) { cons->firstSection--; }
  
    while ((tubeSection[ cons->lastSection ].area < maxArea) &&
           (tubeSection[ cons->lastSection ].articulator == Tube::LOWER_LIP) &&
           (cons->lastSection < Tube::LAST_MOUTH_SECTION)) { cons->lastSection++; }

    cons->firstSection++;
    cons->lastSection--;

    // The obstacle pos. is right at the end of the constriction.
    cons->obstaclePos = tubeSection[cons->lastSection + 1].pos;
  }

    
  // ****************************************************************
  // Determine the parameters of the dipole noise source(s) for all
  // individual constrictions.
  // Always use two dipole sources: One at the upstream end of the
  // tube section with the constriction, and one at the downstream 
  // end. The amplitude of both sources is scaled in relation to the
  // distances of both tube ends to the constriction location.
  // ****************************************************************

  for (k=0; k < numConstrictions; k++)
  {
    cons = &constriction[k];

    // Determine the tube section with the obstacle.

    int obstacleSection = -1;

    for (i = Tube::FIRST_PHARYNX_SECTION; (i <= Tube::LAST_MOUTH_SECTION) && (obstacleSection == -1); i++)
    {
      if ((tubeSection[i].pos <= cons->obstaclePos) &&
          (tubeSection[i].pos + tubeSection[i].length >= cons->obstaclePos))
      {
        obstacleSection = i;
      }
    }

    if (obstacleSection != -1)
    {
      // Set the parameters of the two noise sources.

      NoiseSource *upstreamSource = NULL;
      NoiseSource *downstreamSource = NULL;

      if (obstacleSection < Tube::LAST_MOUTH_SECTION) 
      { 
        upstreamSource = &tubeSection[obstacleSection].dipoleSource; 
        downstreamSource = &tubeSection[obstacleSection + 1].dipoleSource; 
      }
      else
      { 
        upstreamSource = &tubeSection[obstacleSection].dipoleSource; 
        downstreamSource = &lipsDipoleSource; 
      }

      // Factors between 0 and 1 for the contributions of the two sources.
      double downstreamFactor = (cons->obstaclePos - tubeSection[obstacleSection].pos) / 
        tubeSection[obstacleSection].length;
      double upstreamFactor = 1.0 - downstreamFactor;
     
      // Determine the particle velocity and the area of the constriction.

      // This must be greater than the "normal" minimum area
      // to avoid click artifacts in the noise sources.
      const double MIN_AREA_CM2 = 0.1;    // 10 mm^2

      ts = &tubeSection[cons->narrowestSection];
      double A = ts->area;

      if (A < MIN_AREA_CM2)
      {
        A = MIN_AREA_CM2;
      }

      double flow = 0.0;
      if (ts->currentOut[0] != -1) 
      { 
        flow+= branchCurrent[ts->currentOut[0]].noiseMagnitude;
      }
      if (ts->currentOut[1] != -1) 
      { 
        flow+= branchCurrent[ts->currentOut[1]].noiseMagnitude; 
      }

      // Generate noise sources only for outgoing flow.
      // Otherwise we might get click artifacts.
      if (flow < 0.0)
      {
        flow = 0.0;
      }

      double v = flow / A;
      double cutoffFreq = 6000.0;
      double gain = 0.0;

      if (cons->articulator == Tube::LOWER_LIP)
      {
        // Make /f/ less noisy and give it an essentially flat spectrum 
        // (see Badin 1989).
        gain = 2.0e-7;
        cutoffFreq = 6000.0;
      }
      else
      if (cons->articulator == Tube::VOCAL_FOLDS)
      {
        // The gain K must vary between 0.05 for good breathy phonation
        // (not too much aspiration noise) and 3.0 for good aspiration
        // after voiceless aspirated plosives.
        // Hence there is a range of variation of 
        // X = 20*log_10(3.0/0.05) = 35 dB.
        gain = 0.5e-7 * pow(10.0, aspirationStrength_dB / 20.0);

        // According to Badin et al. (1994), the aspiration noise source
        // spectrum is essentially flat -> therefore set a very high
        // cutoff frequency.
        cutoffFreq = 6000.0;
      }
      else
      {
        // "Normal" fricatives: /s,sch,ch,x/
        
        // If the source is at the teeth, it is louder compared to a wall source.
        if (fabs(cons->obstaclePos - teethPosition) < 0.0001)
        {
          double d = sqrt(4.0*A / M_PI);
          cutoffFreq = 0.15*v/d;
          gain = 10.0e-7;
        }
        else
        {
          // The case for a wall source for /x/
          gain = 5.0e-7;
          double d = sqrt(4.0*A / M_PI);
          cutoffFreq = 0.15*v/d;
        }
      }

      double fullAmp = gain*fabs(v)*v*v*sqrt(A);   // Stevens' book

      // Strongly reduce the amplitude if this is a lateral constriction.

      const double LATERALITY_THRESHOLD = 0.1;
      if (cons->laterality > LATERALITY_THRESHOLD)
      {
        fullAmp = 0.0;
      }

      // Safety check:
      // According to Steven's book (p. 105), f_c is typically rather low:
      // e.g. 600 Hz for u = 160 cm^3/s and 1600 Hz for u = 416 cm^3/s
      if (cutoffFreq < 50.0)
      {
        cutoffFreq = 50.0;
      }
      if (cutoffFreq > 2000.0)
      {
        cutoffFreq = 2000.0;
      }

      upstreamSource->targetAmp = upstreamFactor * fullAmp;   
      upstreamSource->isFirstOrder = true;
      upstreamSource->cutoffFrequency = cutoffFreq;

      downstreamSource->targetAmp = downstreamFactor * fullAmp;   
      downstreamSource->isFirstOrder = true;
      downstreamSource->cutoffFrequency = cutoffFreq;
    }     // of (obstacleSection != -1)

  }

  // ****************************************************************
  // Calculate the new noise samples at the positions of the sources.
  // ****************************************************************

  const double MIN_MONOPOLE_AMP = 0.001;          // cm^3/s
  const double MIN_DIPOLE_AMP = 0.001;          // deci-Pascal

  // Run through all tube sections.

  for (i = Tube::FIRST_PHARYNX_SECTION; i <= Tube::LAST_MOUTH_SECTION; i++)
  {
    ts = &tubeSection[i];
    calcNoiseSample(&ts->monopoleSource, MIN_MONOPOLE_AMP);
    calcNoiseSample(&ts->dipoleSource, MIN_DIPOLE_AMP);
  }
  calcNoiseSample(&lipsDipoleSource, MIN_DIPOLE_AMP);
}


// ****************************************************************************
/// Calculates a new noise sample for the given noise source and the given
/// amplitude threshold.
// ****************************************************************************

void TdsModel::calcNoiseSample(NoiseSource *s, double ampThreshold)
{
  // ****************************************************************
  // Force a smooth onset and offset of the noise source amplitudes 
  // using a 1st oder low-pass filter.
  // ****************************************************************

  const double FILTER_CUTOFF_FREQ = 40.0;
  const double F = 1.0 - exp(-2.0*M_PI*FILTER_CUTOFF_FREQ*timeStep);

  double oldAmp = s->currentAmp;
  s->currentAmp += F*(s->targetAmp - s->currentAmp);   // 1st order approach

  // ****************************************************************
  // The source has just turned off.
  // ****************************************************************

  if ((oldAmp >= ampThreshold) && (s->currentAmp < ampThreshold))
  {
    // Clear the noise sample buffers.
    int k;
    for (k = 0; k < NUM_NOISE_BUFFER_SAMPLES; k++)
    {
      s->inputBuffer[k] = 0.0;
      s->outputBuffer[k] = 0.0;
    }
  }

  // ****************************************************************
  // The source is not active -> return now.
  // ****************************************************************

  if (s->currentAmp < ampThreshold)
  {
    s->sample = 0.0;
    return;
  }

  // ****************************************************************
  // The source is active.
  // ****************************************************************
  
  // Setup the IIR-filter to get the filter coefficients.
  IirFilter filter;

  if (s->isFirstOrder)
  {
    filter.createSinglePoleLowpass(s->cutoffFrequency*timeStep);
  }
  else
  {
    double D = s->dampingFactor;
    if (D < 0.000001) { D = 0.000001; }
    filter.createSecondOrderLowpass(s->cutoffFrequency*timeStep, 1.0 / (2.0*D));
  }

  // Insert a new random number into the input buffer and
  // do the recursive filtering.

  // inputSample is a random number with the standard deviation sqrt(12)
  double inputSample = 
    rand() + rand() + rand() + rand() + rand() + rand() +
    rand() + rand() + rand() + rand() + rand() + rand();

  inputSample /= (double)RAND_MAX;
  inputSample -= 6.0;
  inputSample /= sqrt(12.0);    // Divide by the standard deviation of the above sum

  s->inputBuffer[position & NOISE_BUFFER_MASK] = inputSample;
  double sum = filter.a[0] * inputSample;
  int k;
  for (k = 1; k <= filter.order; k++)
  {
    sum += filter.a[k] * s->inputBuffer[(position - k) & NOISE_BUFFER_MASK];
    sum += filter.b[k] * s->outputBuffer[(position - k) & NOISE_BUFFER_MASK];
  }
  s->outputBuffer[position & NOISE_BUFFER_MASK] = sum;
  s->sample = sum * s->currentAmp;      // Resulting magnitude
}


// ****************************************************************************
/// Returns the current flow in the center of the given section.
// ****************************************************************************

void TdsModel::getSectionFlow(int sectionIndex, double &inflow, double &outflow)
{
  inflow = getCurrentIn(sectionIndex);
  outflow = getCurrentOut(sectionIndex);
}

// ****************************************************************************
/// Returns the current pressure in the center of the given section.
// ****************************************************************************

double TdsModel::getSectionPressure(int sectionIndex)
{
  if ((sectionIndex < 0) ||(sectionIndex >= Tube::NUM_SECTIONS))
  {
    return 0.0;
  }

  TubeSection *ts = &tubeSection[sectionIndex];
  return ts->pressure;
}


// ****************************************************************************
/// Calculates the junction inductance between two adjacent tube sections with
/// the given areas according to SONDHI (1983).
/// This additional inductance is a kind of "inner length correction" that
/// affects the formant frequencies when there are abrupt changes in the area
/// function.
// ****************************************************************************

double TdsModel::getJunctionInductance(double A1_cm2, double A2_cm2)
{
  double a, b;    // Radii of the bigger and smaller tube section

  // Safety checks.

  if (A1_cm2 < MIN_AREA_CM2)
  {
    A1_cm2 = MIN_AREA_CM2;
  }
  
  if (A2_cm2 < MIN_AREA_CM2)
  {
    A2_cm2 = MIN_AREA_CM2;
  }

  // ****************************************************************

  if (A1_cm2 > A2_cm2)
  {
    a = sqrt(A1_cm2 / M_PI);
    b = sqrt(A2_cm2 / M_PI);
  }
  else
  {
    a = sqrt(A2_cm2 / M_PI);
    b = sqrt(A1_cm2 / M_PI);
  }

  double H = 1.0 - b/a;
  double L = 8.0*AMBIENT_DENSITY_CGS*H / (3.0*M_PI*M_PI*b);

  return L;
}


// ****************************************************************************
/// Calculate the coefficients of the matrix for the system of eqs.
// ****************************************************************************

void TdsModel::calcMatrix()
{
  int i;

  // Clear the solution vector. *************************************

  for (i=0; i < NUM_BRANCH_CURRENTS; i++) 
  { 
    solutionVector[i] = 0.0; 
  }

  // ****************************************************************
  // Fill one row of the matrix for each branch current.
  // ****************************************************************

  double F, G, H;
  BranchCurrent *bc       = NULL;
  TubeSection   *sourceTs = NULL;
  TubeSection   *targetTs = NULL;

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    bc = &branchCurrent[i];
    if (bc->sourceSection != -1) 
    { 
      sourceTs = &tubeSection[bc->sourceSection]; 
    } 
    else 
    { 
      sourceTs = NULL; 
    }
    
    if (bc->targetSection != -1) 
    { 
      targetTs = &tubeSection[bc->targetSection]; 
    } 
    else 
    { 
      targetTs = NULL; 
    }

    // **************************************************************
    // Both the source and the target section are INvalid.
    // **************************************************************

    if ((!sourceTs) && (!targetTs))
    {
      printf("Error in calcMatrix(): Both the source and target section of the "
        "branch current are invalid!\n");
      return;
    }

    // **************************************************************
    // The branch current flows into free space.
    // **************************************************************

    if (targetTs == NULL)
    {
      // ************************************************************
      // There is no radiation impedance.
      // ************************************************************

      if ((sourceTs->currentOut[0] == -1) || (sourceTs->currentOut[1] == -1))
      {
        printf("Error in calcMatrix(): There are no 2 parallel currents for "
          "the radiation impedance.\n");
        return;
      }

      int resistanceCurrent  = sourceTs->currentOut[0];
      int inductivityCurrent = sourceTs->currentOut[1];

      double uR = branchCurrent[resistanceCurrent].magnitude;
      double uL = branchCurrent[inductivityCurrent].magnitude;
      double uR_rate = branchCurrent[resistanceCurrent].magnitudeRate;
      double uL_rate = branchCurrent[inductivityCurrent].magnitudeRate;

      double L_A = sourceTs->L;
      double R_A = sourceTs->R[1];
      double S   = -lipsDipoleSource.sample;    // Spannungsquelle am Ende des Rohrabschnitts

      double radiationArea = sourceTs->area;

      // ************************************************************
      // Current through the radiation resistor.
      // ************************************************************

      if (i == resistanceCurrent)
      {
        double R_rad = (128 * AMBIENT_DENSITY_CGS * SOUND_VELOCITY_CGS) / 
          (9.0 * M_PI * M_PI * radiationArea);

        F = L_A / (timeStep*THETA) + R_A + R_rad;
        G = L_A / (timeStep*THETA) + R_A;
        H = -(L_A / (timeStep*THETA)) * (uR+uL) - (L_A*(THETA1/THETA))*(uR_rate + uL_rate) + S;
      }
      else

      // ************************************************************
      // Current through the radiation inductivity.
      // ************************************************************

      if (i == inductivityCurrent)
      {
        double L_rad = (8.0*AMBIENT_DENSITY_CGS) / (3.0*M_PI*sqrt(radiationArea*M_PI));
        double L_AB = L_A + L_rad;

        F = L_A/(timeStep*THETA) + R_A;
        G = L_AB/(timeStep*THETA) + R_A;
        H = -(1.0/(timeStep*THETA))*(L_A*uR + L_AB*uL) - (THETA1/THETA)*(L_A*uR_rate + L_AB*uL_rate) + S;
      }
      else
      {
        printf("Error in calcMatrix(): The branch current into the free field has "
          "not a valid type.\n");
        return;
      }

      // Inflowing currents
      if (sourceTs->currentIn != -1) { matrix[i][sourceTs->currentIn] = sourceTs->E; }
        
      // Branch currents through the inductivity and the resistance
      matrix[i][resistanceCurrent]  = -sourceTs->E - F;
      matrix[i][inductivityCurrent] = -sourceTs->E - G;

      solutionVector[i] = H - sourceTs->D;
    }

    else

    // **************************************************************
    // The target tube section is valid in any case.
    // **************************************************************
    {
      // Inductivities, resistances, ...

      double L_B = targetTs->L;
      double R_B = targetTs->R[0];

      double L_A = 0.0;
      double R_A = 0.0;
      if (sourceTs != NULL)
      {
        L_A = sourceTs->L;
        R_A = sourceTs->R[1];
      }

      double L_AB = L_A + L_B;
      double R_AB = R_A + R_B;

      // Is there a birfurcation between the sections A and B ?

      int branchingOffCurrent = -1;

      if (sourceTs != NULL)
      {
        if (sourceTs->currentOut[0] == i) 
          { branchingOffCurrent = sourceTs->currentOut[1]; }
        else
          { branchingOffCurrent = sourceTs->currentOut[0]; }
      }

      // Strength of an additional pressure source between both sections
      
      double S = targetTs->S;
      S-= targetTs->dipoleSource.sample;

      if (bc->targetSection == pressureSourceSection) 
      { 
        S-= pressureSourceAmp; 
      }

      // ************************************************************
      // There is a parallel outflowing current.
      // ************************************************************

      if (branchingOffCurrent != -1)
      {
        double uB      = bc->magnitude;
        double uB_rate = bc->magnitudeRate;
        double uD;
        double uD_rate;

        uD      = branchCurrent[branchingOffCurrent].magnitude;
        uD_rate = branchCurrent[branchingOffCurrent].magnitudeRate;

        F = L_AB/(timeStep*THETA) + R_AB;
        G = L_A/(timeStep*THETA)  + R_A;
        H = - (1.0/(timeStep*THETA))*(L_AB*uB + L_A*uD)
            - (THETA1/THETA)*(L_AB*uB_rate + L_A*uD_rate) + S;

        matrix[i][branchingOffCurrent] = -sourceTs->E - G;    // the parallel current.

        // In A inflowing currents
        if (sourceTs != NULL)
        {
          if (sourceTs->currentIn != -1) { matrix[i][sourceTs->currentIn] = sourceTs->E; }
        }

        // This current
        matrix[i][i] = -targetTs->E - sourceTs->E - F;

        // From B outflowing currents
        if (targetTs->currentOut[0] != -1) { matrix[i][targetTs->currentOut[0]] = targetTs->E; }
        if (targetTs->currentOut[1] != -1) { matrix[i][targetTs->currentOut[1]] = targetTs->E; }

        // Solution value
        solutionVector[i] = H + targetTs->D - sourceTs->D;
      }
      else

      // ************************************************************
      // Simple current from section A to B.
      // ************************************************************

      {
        double u      = bc->magnitude;
        double u_rate = bc->magnitudeRate;

        // Apply the "inner tube length correction" to the junction 
        // between these two sections in terms of an additional
        // inductivity (SONDHI, 1983).

        if ((options.innerLengthCorrections) && (bc->sourceSection >= Tube::FIRST_PHARYNX_SECTION) &&
          (bc->targetSection <= Tube::LAST_MOUTH_SECTION))
        {
          L_AB+= getJunctionInductance(tubeSection[bc->sourceSection].area, tubeSection[bc->targetSection].area);
        }

        G = L_AB / (timeStep*THETA) + R_AB;
        H = - u_rate*L_AB*(THETA1/THETA) - (L_AB*u) / (timeStep*THETA) + S;

        // In A inflowing currents
        if (sourceTs != NULL)
        {
          if (sourceTs->currentIn != -1) { matrix[i][sourceTs->currentIn] = sourceTs->E; }
        }

        // This current
        matrix[i][i] = -targetTs->E - G;
        if (sourceTs != NULL) { matrix[i][i]-= sourceTs->E; }

        // From B outflowing currents
        if (targetTs->currentOut[0] != -1) { matrix[i][targetTs->currentOut[0]] = targetTs->E; }
        if (targetTs->currentOut[1] != -1) { matrix[i][targetTs->currentOut[1]] = targetTs->E; }

        // Solution value
        solutionVector[i] = H + targetTs->D;
        if (sourceTs != NULL) { solutionVector[i]-= sourceTs->D; }
      }

    }
  }   // All branch currents

}


// ****************************************************************************
/// Recalculate all currents, pressures and their temporal derivatives.
// ****************************************************************************

void TdsModel::updateVariables()
{
  int i;
  TubeSection *ts = NULL;
  BranchCurrent *bc = NULL;

  double oldPressure;
  double oldCurrent;
  double oldCurrentRate;
  double netFlow;

  // ****************************************************************
  // The new currents and their derivatives
  // ****************************************************************

  double noiseFilterCoeff = exp(-2.0*M_PI*NOISE_CUTOFF_FREQ*timeStep);

  for (i=0; i < NUM_BRANCH_CURRENTS; i++)
  {
    bc = &branchCurrent[i];
    oldCurrent = bc->magnitude;
    oldCurrentRate = bc->magnitudeRate;

    bc->magnitude = flowVector[i];
    bc->magnitudeRate = (bc->magnitude - oldCurrent)/(timeStep*THETA) - (THETA1/THETA)*bc->magnitudeRate;
    bc->noiseMagnitude = (1.0-noiseFilterCoeff)*bc->magnitude + noiseFilterCoeff*bc->noiseMagnitude;
  }

  // ****************************************************************
  // The new pressures and their derivatives
  // ****************************************************************

  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];

    netFlow = getCurrentIn(ts) - getCurrentOut(ts);

    oldPressure = ts->pressure;
    ts->pressure = ts->D + ts->E*netFlow;
    ts->pressureRate = (ts->pressure - oldPressure)/(timeStep*THETA) - ts->pressureRate*(THETA1/THETA);

    // The current "into the wall".

    oldCurrent = ts->wallCurrent;
    oldCurrentRate = ts->wallCurrentRate;
    
    ts->wallCurrent = ts->pressureRate*ts->alpha + ts->beta;
    ts->wallCurrentRate  = (ts->wallCurrent - oldCurrent)/(timeStep*THETA) - oldCurrentRate*(THETA1/THETA);
    ts->wallCurrentRate2 = (ts->wallCurrentRate - oldCurrentRate)/(timeStep*THETA) - ts->wallCurrentRate2*(THETA1/THETA);
  }

}

// ****************************************************************************
/// Solve the linear system of equations using the SOR method (Gauss-Seidel 
/// method with over-relaxation).
// ****************************************************************************

void TdsModel::solveEquationsSor(const string &matrixFileName)
{
  const int MAX_ITERATIONS = 100; // 100 is the original;
  const double EPSILON = 0.1;    // 0.01 = Original

  // OMEGA=1 corresponds to the Gauss-Seidel-method; The Optimum for OMEGA
  // with 1 <= OMEGA <= 2 must be determined experimetally
  const double OMEGA = 1.25;   // Relaxation parameter

  int i, j, k;
  int iteration;
  double sum;
  TubeSection *ts = NULL;
  double d;
  double residualNorm;      // The magnitude of the residuum squared
  bool isActive[NUM_BRANCH_CURRENTS];

  // ****************************************************************
  // All currents that flow into or out of a closed section are marked
  // as inactive and not calculated.
  // ****************************************************************

  for (i=0; i < NUM_BRANCH_CURRENTS; i++) 
  { 
    isActive[i] = true; 
  }
  
  for (i=0; i < Tube::NUM_SECTIONS; i++)
  {
    ts = &tubeSection[i];
    
    if (ts->area <= 1.01*MIN_AREA_CM2)
    {
      if (ts->currentIn != -1)     { isActive[ts->currentIn] = false; }
      if (ts->currentOut[0] != -1) { isActive[ts->currentOut[0]] = false; }
      if (ts->currentOut[1] != -1) { isActive[ts->currentOut[1]] = false; }
    }
  }
  
  // The initial solution vector.

  for (i=0; i < NUM_BRANCH_CURRENTS; i++) 
  { 
    flowVector[i] = 0.0; 
  }

  // Perform a couple of iterations *********************************
  iteration = 0;

  do
  {
    // Run through all rows of the matrix ***************************
    residualNorm = 0.0;

    for (i=NUM_BRANCH_CURRENTS-1; i >= 0; --i)
    {
      // Shall this current be calculated at all? 
      if (isActive[i])
      {
        // Calculation of the row sum
        sum = matrix[i][i]*flowVector[i];     // Element on the main diagonal

        for (k=numFilledRowValues[i]-1; k >= 0; --k)
        {
          j = filledRowIndex[i][k];
          sum+= flowVector[j]*matrix[i][j];
        }

        d = solutionVector[i] - sum;    // d is the difference between the Shall- and Is-result
        residualNorm+= d*d;
        flowVector[i]+= OMEGA*d/matrix[i][i];
      }
    }

    iteration++;
  }
  while ((iteration < MAX_ITERATIONS) && (residualNorm > EPSILON*EPSILON));


  // For an external evaluation only.
 
  SORIterations = iteration;

  // **************************************************************************
  // If matrixFileName != NULL, write the coefficient matrix, the solution vector,
  // and the RHS into a text file.
  // **************************************************************************

  if (matrixFileName != "")
  {
    ofstream file;
    file.open(matrixFileName);
    if (file)
    {
      printf("Writing matrix to file %s.\n", matrixFileName.c_str());

      int precision = std::numeric_limits<double>::max_digits10 + 1;
      file << std::setprecision(precision);

      file << "Num. iterations: " << iteration << endl;
      file << "Matrix (" << NUM_BRANCH_CURRENTS << " columns); vector of unknowns(1 column); constant terms(right - hand side, 1 column)" << endl;

      int x, y;
      for (y = 0; y < NUM_BRANCH_CURRENTS; y++)
      {
        for (x = 0; x < NUM_BRANCH_CURRENTS; x++)
        {
          file << matrix[y][x] << " ";
        }
        file << flowVector[y] << " ";
        file << solutionVector[y] << " ";
        file << endl;
      }
    }
    else
    {
      printf("Error: Failed to open matrix file %s for writing!\n", matrixFileName.c_str());
    }
  }
}


// ****************************************************************************
/// Solve the linear system of equations using the Cholesky factorization
// ****************************************************************************

void TdsModel::solveEquationsCholesky()
{
  int k, i, j, u, v;

  // ****************************************************************
  // Negate the whole matrix and the right-hand side vector.
  // ****************************************************************

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    factorizationMatrix[i][i] = -matrix[i][i];
    // calculate row i without main diagonal elements
    for (v = 0; v < numFilledRowValuesSymmetricEnvelope[i]; v++)
    {
      j = filledRowIndexSymmetricEnvelope[i][v];
      factorizationMatrix[i][j] = -matrix[i][j];
    }
  }

  for (i = 0; i < NUM_BRANCH_CURRENTS; i++)
  {
    solutionVector[i] = -solutionVector[i];
  }

  // ****************************************************************
  // Cholesky factorization
  // ****************************************************************

  for (k = 0; k < NUM_BRANCH_CURRENTS; k++)
  {
    for (v = 0; v < numFilledRowValuesSymmetricEnvelope[k]; v++)
    {
      j = filledRowIndexSymmetricEnvelope[k][v];
      factorizationMatrix[k][k] -= factorizationMatrix[k][j] * factorizationMatrix[k][j];
    }

    if (factorizationMatrix[k][k] < 0) printf("Error: Cholesky factorization: Matrix is not positive definite!\n");
    factorizationMatrix[k][k] = sqrt(factorizationMatrix[k][k]);

    // calculate column k without main diagonal elements
    for (u = 0; u < numFilledColumnValuesSymmetricEnvelope[k]; u++)
    {
      i = filledColumnIndexSymmetricEnvelope[k][u];

      // calculate row k without main diagonal elements
      for (v = 0; v < numFilledRowValuesSymmetricEnvelope[k]; v++)
      {
        j = filledRowIndexSymmetricEnvelope[k][v];
        factorizationMatrix[i][k] -= factorizationMatrix[i][j] * factorizationMatrix[k][j];
      }

      factorizationMatrix[i][k] = factorizationMatrix[i][k] / factorizationMatrix[k][k];
    }

  }

  // ****************************************************************
  // forward substitution
  // ****************************************************************

  for (k = 0; k < NUM_BRANCH_CURRENTS; k++)
  {
    for (v = 0; v < numFilledRowValuesSymmetricEnvelope[k]; v++)
    {
      j = filledRowIndexSymmetricEnvelope[k][v];
      solutionVector[k] -= factorizationMatrix[k][j] * solutionVector[j];
    }
    solutionVector[k] /= factorizationMatrix[k][k];
  }

  // ****************************************************************
  // backward substitution
  // ****************************************************************

  for (k = NUM_BRANCH_CURRENTS - 1; k >= 0; --k)
  {
    for (u = 0; u < numFilledColumnValuesSymmetricEnvelope[k]; u++)
    {
      i = filledColumnIndexSymmetricEnvelope[k][u];
      solutionVector[k] -= factorizationMatrix[i][k] * flowVector[i];
    }
    flowVector[k] = solutionVector[k] / factorizationMatrix[k][k];
  }
}


// ****************************************************************************
/// Returns the volume velocity into the given tube section.
/// \param section The tube section
// ****************************************************************************

double TdsModel::getCurrentIn(const int section)
{
  double flow = 0.0;
  if ((section >= 0) && (section < Tube::NUM_SECTIONS))
  {
    TubeSection *ts = &tubeSection[section];
    if (ts->currentIn != -1) 
    { 
      flow+= branchCurrent[ts->currentIn].magnitude; 
    }
  }
  return flow;
}

// ****************************************************************************
/// Returns the volume velocity out of the given tube section.
/// \param section The tube section
// ****************************************************************************
    
double TdsModel::getCurrentOut(const int section)
{
  double flow = 0.0;
  if ((section >= 0) && (section < Tube::NUM_SECTIONS))
  {
    TubeSection *ts = &tubeSection[section];
    if (ts->currentOut[0] != -1) { flow+= branchCurrent[ts->currentOut[0]].magnitude; }
    if (ts->currentOut[1] != -1) { flow+= branchCurrent[ts->currentOut[1]].magnitude; }
  }
  return flow;
}

// ****************************************************************************
/// Returns the volume velocity into the given tube section.
/// \param ts Pointer to the tube section
// ****************************************************************************

double TdsModel::getCurrentIn(const TubeSection *ts)
{
  double flow = 0.0;
  if (ts != NULL)
  {
    if (ts->currentIn != -1) 
    { 
      flow+= branchCurrent[ts->currentIn].magnitude; 
    }
  }
  return flow;
}

// ****************************************************************************
/// Returns the volume velocity out of the given tube section.
/// \param ts Pointer to the tube section
// ****************************************************************************
    
double TdsModel::getCurrentOut(const TubeSection *ts)
{
  double flow = 0.0;
  if (ts != NULL)
  {
    if (ts->currentOut[0] != -1) { flow+= branchCurrent[ts->currentOut[0]].magnitude; }
    if (ts->currentOut[1] != -1) { flow+= branchCurrent[ts->currentOut[1]].magnitude; }
  }
  return flow;
}


// ****************************************************************************
