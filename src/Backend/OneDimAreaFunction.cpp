#include "OneDimAreaFunction.h"
#include <float.h>
#include <iostream>
#include <math.h>
#include <wx/msgdlg.h>
#include <wx/tokenzr.h>
// ****************************************************************************
/// Constructor. Initializes the area function.
// ****************************************************************************
OneDimAreaFunction::OneDimAreaFunction()
{
	hasRefAreaFunction = false;
	for (int i = 0; i < NUM_AF_PARAMS; i++)
	{
		isFixed[i] = false;
	}
	reset();
}

// *******************************************************************************
/// Calculate the area at a given position in the vocal tract.
// *******************************************************************************
double OneDimAreaFunction::calculateArea(double pos_cm)
{
	double area_cm2 = 0.0;

	if (pos_cm <= afParams[LLAR_CM].x)
	{
		area_cm2 = afParams[ALAR_CM2].x;
	}
	else if (pos_cm <= afParams[XP_CM].x)
	{
		area_cm2 = (afParams[AP_CM2].x + afParams[ALAR_CM2].x) / 2 + (afParams[AP_CM2].x - afParams[ALAR_CM2].x) / 2 * cos(M_PI*pow((afParams[XP_CM].x - pos_cm) / (afParams[XP_CM].x - afParams[LLAR_CM].x), afParams[POWLAR].x));
	}
	else if (pos_cm <= afParams[XC_CM].x)
	{
		area_cm2 = (afParams[AC_CM2].x + afParams[AP_CM2].x) / 2 + (afParams[AC_CM2].x - afParams[AP_CM2].x) / 2 * cos(M_PI*pow((afParams[XC_CM].x - pos_cm) / (afParams[XC_CM].x - afParams[XP_CM].x), afParams[POWP].x));
	}
	else if (pos_cm <= afParams[XA_CM].x)
	{
		area_cm2 = (afParams[AA_CM2].x + afParams[AC_CM2].x) / 2 + (afParams[AA_CM2].x - afParams[AC_CM2].x) / 2 * cos(M_PI*pow((afParams[XA_CM].x - pos_cm) / (afParams[XA_CM].x - afParams[XC_CM].x), afParams[POWC].x));
	}
	else if (pos_cm <= afParams[XIN_CM].x)
	{
		area_cm2 = (afParams[AIN_CM2].x + afParams[AA_CM2].x) / 2 + (afParams[AIN_CM2].x - afParams[AA_CM2].x) / 2 * cos(M_PI*pow((afParams[XIN_CM].x - pos_cm) / (afParams[XIN_CM].x - afParams[XA_CM].x), afParams[POWA].x));
	}
	else
	{
		area_cm2 = afParams[ALIP_CM2].x;
	}
	
	if (area_cm2 < 0.0)
	{
		area_cm2 = 0.0;
	}

	return area_cm2;
}

// *******************************************************************************
/// Calculate the "continuous" vocal tract area function.
// *******************************************************************************
void OneDimAreaFunction::calculateHiResAreaFunction()
{
	deltaX_cm = afParams[LVT_CM].x / NUM_AREA_SAMPLES;
	for (int i = 0; i < NUM_AREA_SAMPLES; i++)
	{
		areaFunction_cm2[i] = calculateArea(i*deltaX_cm);
	}
}

// *******************************************************************************
/// Calculate the vocal tract area tube function based on the model AF.
// *******************************************************************************
void OneDimAreaFunction::calculateOneDimTubeFunction(Tube *tubeFunction, bool alignTubesWithRef)
{
	// *******************************************************************************
	// The tube function takes the minimum value of the continuous function in the 
	// following interval.
	// *******************************************************************************
	double secWidth_cm = afParams[LVT_CM].x / NUM_SECTIONS;
	double minArea_cm2 = DBL_MAX;
	double currentArea_cm2 = 0.0;
	double stepSize = secWidth_cm*0.01;
		
	double tubeArea_cm2[NUM_SECTIONS] = { 0.0 };
	double length_cm[NUM_SECTIONS] = { 0.0 };
	Tube::Articulator articulator[NUM_SECTIONS] = { Tube::Articulator::OTHER_ARTICULATOR };
	double laterality[NUM_SECTIONS] = { 0.0 };

	int j = 0;
	double x_cm = 0.0;

	for (int i = 0; i < NUM_SECTIONS; i++)
	{
		minArea_cm2 = DBL_MAX;		
		while (x_cm < (i + 1)*secWidth_cm)
		{
			
			currentArea_cm2 = calculateArea(x_cm);
			
			if (currentArea_cm2 < minArea_cm2)
			{
				minArea_cm2 = currentArea_cm2;
			}			

			if (alignTubesWithRef && refX_cm[j] < (i+1)* secWidth_cm)
			{
				x_cm = refX_cm[j];
				j++;
			}
			else
			{
				x_cm += stepSize;
			}
		}
		length_cm[i] = secWidth_cm;
		tubeArea_cm2[i] = minArea_cm2;
		if (x_cm <= afParams[XP_CM].x)
		{
			articulator[i] = Tube::Articulator::OTHER_ARTICULATOR;
		}
		else if (x_cm <= afParams[XIN_CM].x && x_cm + secWidth_cm < afParams[XIN_CM].x)
		{
			articulator[i] = Tube::Articulator::TONGUE;
		}
		else if (x_cm <= afParams[XIN_CM].x)
		{
			articulator[i] = Tube::Articulator::LOWER_INCISORS;
		}
		else if (x_cm > afParams[XIN_CM].x)
		{
			articulator[i] = Tube::Articulator::LOWER_LIP;
		}
	}
	
	tubeFunction->setPharynxMouthGeometry(length_cm, tubeArea_cm2, articulator, &laterality[0], afParams[XIN_CM].x);
}

// *******************************************************************************
/// Calculate a dependent or fixed parameter.
// *******************************************************************************
void OneDimAreaFunction::calculateParameter(Param *param, int paramIdx)
{
	switch (paramIdx)
	{
	case LLAR_CM:
		param[LLAR_CM].x = 2.0;
		break;
	case ALAR_CM2:
		param[ALAR_CM2].x = 1.5125;
		break;
	case POWLAR:
		param[POWLAR].x = 1.0;
		break;
	case XP_CM:
		param[XP_CM].x = 3.0948;
		break;
	case AP_CM2:
		break;
	case POWP:
		break;
	case XC_CM:
		break;
	case AC_CM2:
		break;
	case POWC:
		break;
	case XA_CM:
		param[XA_CM].x = -2.3466
			- 0.0606 * param[XC_CM].x
			- 2.0517 * param[POWC].x
			- 0.1592 * param[AA_CM2].x
			+ 1.1607 * param[XIN_CM].x
			+ 0.1425 * param[XC_CM].x * param[POWC].x;
		break;
	case AA_CM2:
		break;
	case POWA:
		break;
	case XIN_CM:
		break;
	case AIN_CM2:
		break;
	case LVT_CM:
		break;
	case ALIP_CM2:
		break;
	default:
		break;
	}	
	// Make sure x coordinates of the control points are still monotonous
	if (param[XP_CM].x < param[LLAR_CM].x)
	{
		param[XP_CM].x = param[LLAR_CM].x;
	}
	if (param[XC_CM].x < param[XP_CM].x)
	{
		param[XC_CM].x = param[XP_CM].x;
	}
	if (param[XA_CM].x < param[XC_CM].x)
	{
		param[XA_CM].x = param[XC_CM].x;
	}
	if (param[XIN_CM].x < param[XA_CM].x)
	{
		param[XIN_CM].x = param[XA_CM].x;
	}
	if (param[LVT_CM].x < param[XIN_CM].x)
	{
		param[LVT_CM].x = param[XIN_CM].x;
	}

	// Make sure no exponent is negative
	if (param[POWLAR].x < 0)
	{
		param[POWLAR].x = 0;
	}
	if (param[POWP].x < 0)
	{
		param[POWP].x = 0;
	}
	if (param[POWC].x < 0)
	{
		param[POWC].x = 0;
	}
	if (param[POWA].x < 0)
	{
		param[POWA].x = 0;
	}
}

// *******************************************************************************
/// Calculate the vocal tract area tube function based on the reference AF.
// *******************************************************************************
void OneDimAreaFunction::calculateRefTubeFunction(Tube *tubeFunction)
{
	// *******************************************************************************
	// The tube function takes the minimum value of the continuous function in the 
	// following interval.
	// *******************************************************************************
	double secWidth_cm = refX_cm[refLength - 1] / NUM_SECTIONS;
	double minArea_cm2 = DBL_MAX;
	double currentArea_cm2 = 0.0;

	double tubeArea_cm2[NUM_SECTIONS] = { 0.0 };
	double length_cm[NUM_SECTIONS] = { 0.0 };
	double xMax_cm;
	int lastJ = 0;

	Tube::Articulator articulator[NUM_SECTIONS] = { Tube::Articulator::OTHER_ARTICULATOR };

	for (int i = 0; i < NUM_SECTIONS; i++)
	{
		xMax_cm = (i+1)*secWidth_cm;
		minArea_cm2 = DBL_MAX;
		for (int j = lastJ; j < refLength; j++)
		{
			if (refX_cm[j] < xMax_cm)
			{
				currentArea_cm2 = refAreaFunction_cm2[j];
				if (currentArea_cm2 < minArea_cm2)
				{
					minArea_cm2 = currentArea_cm2;
				}
			}
			else
			{
				lastJ = j;
				break;
			}			
		}
		length_cm[i] = secWidth_cm;
		tubeArea_cm2[i] = minArea_cm2;
	}
	double laterality[Tube::NUM_PHARYNX_MOUTH_SECTIONS] = { 0.0 };
	tubeFunction->setPharynxMouthGeometry(length_cm, tubeArea_cm2, articulator, &laterality[0], afParams[XIN_CM].x);
}


void OneDimAreaFunction::reset()
{
	/* Init larynx tube */
	afParams[LLAR_CM].name = "Length of the larynx tube in cm";
	afParams[LLAR_CM].abbr = "Llar_cm";
	afParams[LLAR_CM].x = 2.0;
	afParams[ALAR_CM2].name = "Area of the larynx tube in cm2";
	afParams[ALAR_CM2].abbr = "Alar_cm2";
	afParams[ALAR_CM2].x = 1.0;
	afParams[POWLAR].name = "Exponent of the segment between larynx and posterior control point";
	afParams[POWLAR].abbr = "powLar";
	afParams[POWLAR].x = 1.0;

	/* Init variables for schwa area function */
	afParams[XP_CM].name = "X-coordinate of posterior control point in cm";
	afParams[XP_CM].abbr = "xp_cm";
	afParams[XP_CM].x = 3.02;
	afParams[AP_CM2].name = "Area at posterior control point in cm2";
	afParams[AP_CM2].abbr = "Ap_cm2";
	afParams[AP_CM2].x = 5.609;
	afParams[POWP].name = "Exponent of segment between posterior and constriction control point";
	afParams[POWP].abbr = "powP";
	afParams[POWP].x = 1.0;

	afParams[XC_CM].name = "X-coordinate of constriction control point in cm";
	afParams[XC_CM].abbr = "xc_cm";
	afParams[XC_CM].x = 5.92;
	afParams[AC_CM2].name = "Area at constriction control point in cm^2";
	afParams[AC_CM2].abbr = "Ac_cm2";
	afParams[AC_CM2].x = 2.879;
	afParams[POWC].name = "Exponent of segment between constriction and anterior control point";
	afParams[POWC].abbr = "powC";
	afParams[POWC].x = 1.0;

	afParams[XA_CM].name = "X-coordinate of anterior control point in cm";
	afParams[XA_CM].abbr = "xa_cm";
	afParams[XA_CM].x = 8.48;
	afParams[AA_CM2].name = "Area at anterior control point in cm^2";
	afParams[AA_CM2].abbr = "Aa_cm2";
	afParams[AA_CM2].x = 4.238;
	afParams[POWA].name = "Exponent of segment between anterior and incisor control point";
	afParams[POWA].abbr = "powA";
	afParams[POWA].x = 1.0;

	afParams[XIN_CM].name = "X-coordinate of incisor position in cm";
	afParams[XIN_CM].abbr = "xin_cm";
	afParams[XIN_CM].x = 15.31;
	afParams[AIN_CM2].name = "Area at incisor position in cm2";
	afParams[AIN_CM2].abbr = "Ain_cm2";
	afParams[AIN_CM2].x = 0.701;

	afParams[LVT_CM].name = "Length of vocal tract and end of lip tube in cm";
	afParams[LVT_CM].abbr = "Lvt_cm";
	afParams[LVT_CM].x = 16.44;
	afParams[ALIP_CM2].name = "Area of lip opening in cm2";
	afParams[ALIP_CM2].abbr = "Alip_cm2";
	afParams[ALIP_CM2].x = 1.65;

	calculateHiResAreaFunction();
}

// *******************************************************************************
/// Return the current parameters of the area function.
// *******************************************************************************
void OneDimAreaFunction::getParameters(Param *currentParams)
{
	for (int i = 0; i < NUM_AF_PARAMS; i++)
	{
		currentParams[i] = afParams[i];	
	}	
}

// *******************************************************************************
/// Set the parameters of the area function.
// *******************************************************************************
void OneDimAreaFunction::setParameters(Param *newParams)
{
	for (int i = 0; i < NUM_AF_PARAMS; i++)
	{
		afParams[i].x = newParams[i].x;
	}
	for (int i = 0; i < NUM_AF_PARAMS; i++)
	{
		if (isFixed[i])
		{
			calculateParameter(afParams, i);
		}
	}
	calculateHiResAreaFunction();
}


// *******************************************************************************
/// Set the reference area function.
// *******************************************************************************
bool OneDimAreaFunction::setRefAreaFunction(int length, double *x_cm, double *area_cm2)
{
	if (length <= 0){ return false; }
	refLength = length;
	refX_cm = new double[refLength];
	refAreaFunction_cm2 = new double[refLength];
	
	for (int i = 0; i < refLength; i++)
	{
		refX_cm[i] = x_cm[i];
		refAreaFunction_cm2[i] = area_cm2[i];
	}
	hasRefAreaFunction = true;

	return true;
}

// *******************************************************************************
/// Get the reference area function.
// *******************************************************************************
bool OneDimAreaFunction::getRefAreaFunction(int *length, double **x_cm, double **area_cm2)
{
	if (!hasRefAreaFunction){ return false; }

	*length = refLength;
	*x_cm = new double[refLength];
	*area_cm2 = new double[refLength];

	for (int i = 0; i < refLength; i++)
	{
		*x_cm[i] = refX_cm[i];
		*area_cm2[i] = refAreaFunction_cm2[i];
	}

	return true;
}
