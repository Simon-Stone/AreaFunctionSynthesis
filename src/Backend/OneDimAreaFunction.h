#ifndef __ONE_DIM_AREA_FUNCTION_H__
#define __ONE_DIM_AREA_FUNCTION_H__

#include <wx/textfile.h>
#include "Tube.h"

// ************************************************************************************
/// This is the three-parameter 1D area function according to Fant & Bavegard (1997).
// ************************************************************************************

class OneDimAreaFunction
{
	// **************************************************************************
	// Public data.
	// **************************************************************************

public:
	// ****************************************************************
	// Time-variant area function parameters.
	// ****************************************************************
	struct Param
	{
		double x;         ///< Parameter value set by the user or defined as target
		double limitedX;  ///< Parameter value limited by biomechanical constraints
		double min;
		double max;
		double neutral;
		wxString abbr;      ///< Abbreviation, e.g. "XLAR_CM"
		wxString name;      ///< Long name, e.g. " x position of larynx tube end in cm"
	};

	// ****************************************************************

	enum ParamIndex
	{
		LLAR_CM, ALAR_CM2, POWLAR,
		XP_CM, AP_CM2, POWP,
		XC_CM, AC_CM2, POWC,
		XA_CM, AA_CM2, POWA,
		XIN_CM, AIN_CM2, 
		LVT_CM, ALIP_CM2,		
		NUM_AF_PARAMS
	};

	// ****************************************************************
	// A defined area function, for example for the 
	// vowels [i], [a], [u].
	// ****************************************************************
	struct Shape
	{
		wxString name;
		Param param[NUM_AF_PARAMS];    ///< Parameters
	};

	// Number of samples of the vocal tract area function
	static const int NUM_SECTIONS = Tube::NUM_PHARYNX_MOUTH_SECTIONS;
		// Maximum number of samples in high resolution area function used for drawing
	static const int NUM_AREA_SAMPLES = 100;
	// High-resolution area function for drawing
	double areaFunction_cm2[NUM_AREA_SAMPLES];
	// Step width for high-res aera function
	double deltaX_cm;
  
	// Is reference area function from VocalTractLab loaded
	bool hasRefAreaFunction;
	// Length of reference area function in samples
	int refLength;
	// X coordinate of each area sample in reference function
	double* refX_cm;
	// Reference area values
	double* refAreaFunction_cm2;
	// Flags indicating if corresponding parameter value is fixed or constrained
	bool isFixed[NUM_AF_PARAMS];
	
	// **************************************************************************
	// Public functions.
	// **************************************************************************

public:
	OneDimAreaFunction();

	// Calculates the area at the given x position
	double calculateArea(double x_cm);
	// Calculates the high-res area function
	void calculateHiResAreaFunction();
	// Calculates the discrete tube function
	void calculateOneDimTubeFunction(Tube *tubeFunction, bool alignTubesWithRef = false);
	// Calcualates a dependent parameter 
	void calculateParameter(Param *param, int paramIdx);
	// Calculate the vocal tract area tube function according to the reference AF.
	void calculateRefTubeFunction(Tube *tubeFunction);
	

	void reset();

	// Returns the current paramters of the vocal tract area function
	void getParameters(Param *currentParams);
	void setParameters(Param *newParams);
	
	bool setRefAreaFunction(int length, double *x_cm, double *area_cm2);
	bool getRefAreaFunction(int *length, double **x_cm, double **area_cm2);
	

	// **************************************************************************
	// Private data.
	// **************************************************************************

private:
	// parameter values of the area function model
	Param afParams[NUM_AF_PARAMS];
	
	

	// **************************************************************************
	// Private functions.
	// **************************************************************************

private:
	
};

#endif
