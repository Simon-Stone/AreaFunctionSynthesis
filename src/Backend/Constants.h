#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

// ****************************************************************************
// Physical constants according to Flanagan (1965).
// ****************************************************************************

const double STATIC_PRESSURE_CGS = 1.013e6;   // deci-Pa = ubar
const double AMBIENT_DENSITY_CGS = 1.14e-3;   // g/cm^3
const double ADIABATIC_CONSTANT = 1.4;      
const double SOUND_VELOCITY_CGS = 3.5e4;      // cm/s
const double AIR_VISCOSITY_CGS = 1.86e-4;     // dyne-s/cm^2
const double SPECIFIC_HEAT_CGS = 0.24;        // cal/g-K
const double HEAT_CONDUCTION_CGS = 0.055e-3;  // cal/cm-s-K

const double CRITICAL_REYNOLDS_NUMBER = 1800.0;

// ****************************************************************************
// Constants for the synthesis.
// ****************************************************************************

const int SAMPLING_RATE = 22050;
// By which factor is the sampling rate in the time-domain simulation higher 
// than in the final signal:
const int TDS_SR_FACTOR = 1;
const int TDS_SAMPLING_RATE = TDS_SR_FACTOR*SAMPLING_RATE;
// The synthetic speech signal is finally restricted to this upper frequency:
const double SYNTHETIC_SPEECH_BANDWIDTH_HZ = 8000.0;

#endif

