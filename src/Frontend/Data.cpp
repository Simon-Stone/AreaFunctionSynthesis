#include <wx/filename.h>
#include <wx/stdpaths.h>
#include <wx/progdlg.h>
#include <iostream>
#include <fstream>

#include "Data.h"

// Define a custom event type to be used for command events.

const wxEventType updateRequestEvent = wxNewEventType();
const int REFRESH_PICTURES = 0;
const int UPDATE_PICTURES = 1;
const int REFRESH_PICTURES_AND_CONTROLS = 2;


// ****************************************************************************
// Static data.
// ****************************************************************************

Data *Data::instance = NULL;


// ****************************************************************************
/// Returns the one instance of this class.
// ****************************************************************************

Data *Data::getInstance()
{
  if (instance == NULL)
  {
    instance = new Data();
  }
  return instance;
}


// ****************************************************************************
/// Init the data. This function must be called once after the first call of
/// getInstance().
/// \param arg0 The first string parameter passed to this program.
// ****************************************************************************

void Data::init(wxString arg0)
{
//  int i;

  // ****************************************************************
  // Determine the program path from arg0. The option
  // wxPATH_GET_SEPARATOR makes sure that the path is always
  // terminated with a "\".
  // ****************************************************************

  wxFileName fileName(arg0);
  programPath = fileName.GetPath(wxPATH_GET_VOLUME | wxPATH_GET_SEPARATOR);
  printf("The program path is %s.\n", (const char*) programPath.mb_str());

  // ****************************************************************
  // Init variables.
  // ****************************************************************

  selectionMark_s[0] = 0.0;
  selectionMark_s[1] = 1.0;
  mark_s = 0.0;

	// ****************************************************************
	// Oscillogram variables
	// ****************************************************************

	for (int i = 0; i < NUM_TRACKS; i++)
	{
		track[i] = new Signal16(TRACK_DURATION_S*SAMPLING_RATE);
	}

	synthesizer = new Synthesizer(track[MAIN_TRACK]);
    synthesizer->modelControlMode = Synthesizer::MANUAL;

	tlModel = new TlModel();
	tlModel->tube.setGlottisArea(0.0);
	tlModel->tube.setVelumOpening(0.0);
	tlModel->tube.setAspirationStrength(Tube::DEFAULT_ASPIRATION_STRENGTH_DB);
	synthesizer->vtAreaFunction->calculateOneDimTubeFunction(&tlModel->tube, synthesizer->vtAreaFunction->hasRefAreaFunction);
	
	tlRef = new TlModel();
	tlRef->tube.setGlottisArea(0.0);
	tlRef->tube.setVelumOpening(0.0);
	tlRef->tube.setAspirationStrength(Tube::DEFAULT_ASPIRATION_STRENGTH_DB);
	synthesizer->vtAreaFunction->calculateOneDimTubeFunction(&tlRef->tube);


	primarySpectrum = new ComplexSignal();	

//	formantOptimizationDialog = new FormantOptimizationDialog(NULL, synthesizer->vtAreaFunction);
//	setTargetSequenceDialog = new SetTargetSequenceDialog(NULL, synthesizer->vtAreaFunction,
//		&targetShape[0], &stationary_s[0], &transition_s[0]);

	analysisMode = REAL_TIME;
}

// ****************************************************************************
/// Returns the mean squared difference between the current and the target
/// formant values. If targetF3 == 0, only the first two formants are considered
/// for the error calculation.
// ****************************************************************************

double Data::getFormantError(double currentF1, double currentF2,
	double currentF3, double targetF1, double targetF2, double targetF3)
{
	const double EPSILON = 1.0;
	double e = 0.0;

	// Consider all three formants.

	if (targetF3 > EPSILON)
	{
		e = (1.0 - currentF1 / targetF1)*(1.0 - currentF1 / targetF1) +
			(1.0 - currentF2 / targetF2)*(1.0 - currentF2 / targetF2) +
			(1.0 - currentF3 / targetF3)*(1.0 - currentF3 / targetF3);

		e /= 3.0;    // because we have 3 formants
	}
	else

		// Consider only two formants.
	{
		e = (1.0 - currentF1 / targetF1)*(1.0 - currentF1 / targetF1) +
			(1.0 - currentF2 / targetF2)*(1.0 - currentF2 / targetF2);

		e /= 2.0;    // because we have 2 formants
	}

	e = sqrt(e);
	e = e*100.0;    // The result is in percent

	return e;
}

// ****************************************************************************
/// Calculates the first three formants of the given vocal tract and the 
/// minimum area of the corresponding area function. If the minimum area is 
/// below a certain threshold, the formant values may not be useful.
/// \param tract The vocal tract with adjusted variables.
/// \param F1 Returns the first formant frequency.
/// \param F2 Returns the second formant frequency.
/// \param F3 Returns the third formant frequency.
/// \param minArea_cm2 Returns the minimum area of the area function.
// ****************************************************************************

bool Data::getVowelFormants(OneDimAreaFunction *vtAf, double &F1_Hz, double &F2_Hz, double &F3_Hz, double &minArea_cm2)
{
	const int MAX_FORMANTS = 3;
	double formantFreq[MAX_FORMANTS];
	double formantBw[MAX_FORMANTS];
	int numFormants;
	bool frictionNoise;
	bool isClosure;
	bool isNasal;
	int i;

	// Default return values.
	F1_Hz = 0.0;
	F2_Hz = 0.0;
	F3_Hz = 0.0;
	minArea_cm2 = 0.0;

	// Calculate and set tube function
	Tube tubeFunction;
	vtAf->calculateOneDimTubeFunction(&tubeFunction);
	// Set the current area function for the transmission line model.   
	tlModel->tube = tubeFunction;
	// The velo-pharyngeal port must be closed.
	tlModel->tube.setGlottisArea(0.0);

	// Find the minimum cross-sectional area.
	minArea_cm2 = 10000.0;    // = extremely high
	for (i = Tube::FIRST_PHARYNX_SECTION; i < Tube::NUM_PHARYNX_MOUTH_SECTIONS; i++)
	{
		if (tlModel->tube.section[i]->area_cm2 < minArea_cm2)
		{
			minArea_cm2 = tubeFunction.section[i]->area_cm2;
		}
	}

	// Get the formant data.
	tlModel->getFormants(formantFreq, formantBw, numFormants, MAX_FORMANTS, frictionNoise, isClosure, isNasal);
	if (numFormants < MAX_FORMANTS)
	{
		return false;
	}

	F1_Hz = formantFreq[0];
	F2_Hz = formantFreq[1];
	F3_Hz = formantFreq[2];

	return true;
}

// ****************************************************************************
/// Optimize the parameters of the given vocal tract area function so that the formants
/// match the given formant target values as well as possible.
// ****************************************************************************

void Data::optimizeFormantsVowel(wxWindow*, OneDimAreaFunction *vtAf,
	double targetF1, double targetF2, double targetF3,
	double maxParamChange_cm, double minAdvisedArea_cm2, bool paramFixed[])
{
	const int MAX_RUNS = 100;

	double F1, F2, F3;
	double minArea_cm2;
	double changeStep[OneDimAreaFunction::NUM_AF_PARAMS];
	int stepsTaken[OneDimAreaFunction::NUM_AF_PARAMS];   // Cummulated steps gone by a parameter
	double currParamValue;
	double bestError;
	double currError;
	double newError;
	double bestParamChange;
	int bestParam;
	int i;
	wxString st;
	wxCommandEvent event(updateRequestEvent);


	// ****************************************************************
	// Make sure that all areas are above the given threshold.
	// ****************************************************************

	getVowelFormants(vtAf, F1, F2, F3, minArea_cm2);
	//if (minArea_cm2 < minAdvisedArea_cm2)
	//{
	//	wxString st = wxString::Format(
	//		"There are cross-sectional areas below the threshold of %2.2f cm^2.\n"
	//		"Press OK to let it be corrected before the optimization.",
	//		minAdvisedArea_cm2);

	//	if (wxMessageBox(st, "Areas to small", wxOK | wxCANCEL) == wxOK)
	//	{
	//		createMinVocalTractArea(updateParent, tract, minAdvisedArea_cm2);
	//	}
	//	else
	//	{
	//		return;
	//	}
	//}

	double initialError = getFormantError(F1, F2, F3, targetF1, targetF2, targetF3);
	wxPrintf("\n=== Before vowel optimization ===\n");
	wxPrintf("F1:%d  F2:%d  F3:%d   F1':%d  F2':%d  F3':%d   error:%2.2f percent\n",
		(int)F1, (int)F2, (int)F3, (int)targetF1, (int)targetF2, (int)targetF3, initialError);

	// ****************************************************************
	// Define the steps by which the individual parameters are adjusted
	// incrementally. A step value of 0.0 means that the parameters is
	// not supposed to be adjusted at all.
	// ****************************************************************

	// How many steps are maximally allowed per parameter in one direction ?
	const double STEP_SIZE_CM = 0.05;   // = 1/2 mm
	int maxSteps = (int)((maxParamChange_cm / STEP_SIZE_CM) + 0.5);
	if (maxSteps < 1)
	{
		maxSteps = 1;
	}
	if (maxSteps > 40)
	{
		maxSteps = 40;    // Corresponds to 2.0 cm
	}
	wxPrintf("The area function control point positions may change by at most %1.1f mm.\n", maxSteps*STEP_SIZE_CM*10.0);

	for (i = 0; i < OneDimAreaFunction::NUM_AF_PARAMS; i++)
	{
		changeStep[i] = 0.0;
		stepsTaken[i] = 0;
	}

	changeStep[OneDimAreaFunction::LLAR_CM] = STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::XP_CM] = STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::XC_CM] = STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::XA_CM] = STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::XIN_CM] = STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::LVT_CM] = STEP_SIZE_CM;

	changeStep[OneDimAreaFunction::ALAR_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::AP_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::AC_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::AA_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::AIN_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	changeStep[OneDimAreaFunction::ALIP_CM2] = STEP_SIZE_CM * STEP_SIZE_CM;
	
	changeStep[OneDimAreaFunction::POWLAR] = 0.1;
	changeStep[OneDimAreaFunction::POWP] = 0.1;
	changeStep[OneDimAreaFunction::POWC] = 0.1;
	changeStep[OneDimAreaFunction::POWA] = 0.1;
	
	// Set the change step to zero for parameters that are supposed to be fixed.
	for (i = 0; i < OneDimAreaFunction::NUM_AF_PARAMS; i++)
	{
		if (paramFixed[i])
		{
			changeStep[i] = 0.0;
		}
	}

	// ****************************************************************
	// Init. the progress dialog.
	// ****************************************************************

	wxProgressDialog progressDialog("Please wait", "The formant optimization is running...",
		MAX_RUNS, NULL, wxPD_CAN_ABORT | wxPD_APP_MODAL | wxPD_AUTO_HIDE);

	// ****************************************************************
	// ****************************************************************

	bool paramChanged = false;
	bool doContinue = false;
	int runCounter = 0;
	OneDimAreaFunction::Param currentParameters[OneDimAreaFunction::NUM_AF_PARAMS];

	do
	{
		paramChanged = false;
		getVowelFormants(vtAf, F1, F2, F3, minArea_cm2);
		currError = getFormantError(F1, F2, F3, targetF1, targetF2, targetF3);

		// **************************************************************
		// Find out the improvement of the error when each parameter is 
		// changed individually by a positive changeStep[i] starting 
		// from the current configuration.
		// **************************************************************

		bestError = currError;
		bestParam = -1;
		bestParamChange = 0.0;

		vtAf->getParameters(currentParameters);

		for (i = 0; i < OneDimAreaFunction::NUM_AF_PARAMS; i++)
		{
			currParamValue = currentParameters[i].x;

			if (changeStep[i] > 0.0)
			{
				// Apply a POSITIVE change to parameter i.
				if (stepsTaken[i] < maxSteps)
				{
					currentParameters[i].x = currParamValue + changeStep[i];
					/* Check if order of control points is still monotonously increasing */
					if (i == OneDimAreaFunction::XP_CM && currentParameters[i].x > currentParameters[OneDimAreaFunction::XC_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XC_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XC_CM].x;
					}
					if (i == OneDimAreaFunction::XC_CM && currentParameters[i].x > currentParameters[OneDimAreaFunction::XA_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XA_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XA_CM].x;
					}
					if (i == OneDimAreaFunction::XA_CM && currentParameters[i].x > currentParameters[OneDimAreaFunction::XIN_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XIN_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XIN_CM].x;
					}
					if (i == OneDimAreaFunction::XIN_CM && currentParameters[i].x > currentParameters[OneDimAreaFunction::LVT_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::LVT_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::LVT_CM].x;
					}

					vtAf->setParameters(currentParameters);
					getVowelFormants(vtAf, F1, F2, F3, minArea_cm2);

					// Check that the minimum area stays above the threshold.
					if (minArea_cm2 >= minAdvisedArea_cm2)
					{
						newError = getFormantError(F1, F2, F3, targetF1, targetF2, targetF3);
						if (newError < bestError)
						{
							bestError = newError;
							bestParam = i;
							bestParamChange = changeStep[i];
						}
					}
				}

				// Apply a NEGATIVE change to parameter i.
				if (stepsTaken[i] > -maxSteps)
				{
					currentParameters[i].x = currParamValue - changeStep[i];
					/* Check if order of control points is still monotonously increasing */
					if (i == OneDimAreaFunction::XP_CM && currentParameters[i].x < currentParameters[OneDimAreaFunction::LLAR_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::LLAR_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::LLAR_CM].x;
					}
					if (i == OneDimAreaFunction::XC_CM && currentParameters[i].x < currentParameters[OneDimAreaFunction::XP_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XP_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XP_CM].x;
					}
					if (i == OneDimAreaFunction::XA_CM && currentParameters[i].x < currentParameters[OneDimAreaFunction::XC_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XC_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XC_CM].x;
					}
					if (i == OneDimAreaFunction::XIN_CM && currentParameters[i].x < currentParameters[OneDimAreaFunction::XA_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XA_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XA_CM].x;
					}
					if (i == OneDimAreaFunction::LVT_CM && currentParameters[i].x < currentParameters[OneDimAreaFunction::XIN_CM].x)
					{
						changeStep[i] = currentParameters[OneDimAreaFunction::XIN_CM].x - currentParameters[i].x;
						currentParameters[i].x = currentParameters[OneDimAreaFunction::XIN_CM].x;
					}

					vtAf->setParameters(currentParameters);
					getVowelFormants(vtAf, F1, F2, F3, minArea_cm2);
					// Check that the minimum area stays above the threshold.
					if (minArea_cm2 >= minAdvisedArea_cm2)
					{
						newError = getFormantError(F1, F2, F3, targetF1, targetF2, targetF3);
						if (newError < bestError)
						{
							bestError = newError;
							bestParam = i;
							bestParamChange = -changeStep[i];
						}
					}
				}

				// Set the parameter value back to its original value.
				currentParameters[i].x = currParamValue;
				vtAf->setParameters(currentParameters);
			}
		}

		// **************************************************************
		// Change the parameter with the best error reduction.
		// **************************************************************
		st = wxString("no change");

		if ((bestParam != -1) && (bestError < currError))
		{
			currentParameters[bestParam].x += bestParamChange;
			vtAf->setParameters(currentParameters);
			if (bestParamChange > 0.0)
			{
				stepsTaken[bestParam]++;
				st = wxString::Format("%s up", currentParameters[bestParam].abbr);
			}
			else
			{
				stepsTaken[bestParam]--;
				st = wxString::Format("%s down", currentParameters[bestParam].abbr);
			}

			paramChanged = true;
		}

		wxPrintf("Run %d: %s. Error=%2.2f\n", runCounter, st, bestError);

		doContinue = progressDialog.Update(runCounter);

		runCounter++;
		refreshAreaFunction();

	} while ((paramChanged) && (doContinue) && (runCounter < MAX_RUNS));

	// Hide the progress dialog.
	progressDialog.Update(MAX_RUNS);

	// ****************************************************************
	// The velo-pharyngeal port must be closed.
	// ****************************************************************
	tlModel->tube.setGlottisArea(0.0);


	// ****************************************************************
	// ****************************************************************

	getVowelFormants(vtAf, F1, F2, F3, minArea_cm2);
	double finalError = getFormantError(F1, F2, F3, targetF1, targetF2, targetF3);
	wxPrintf("=== After optimization ===\n");
	wxPrintf("F1:%d  F2:%d  F3:%d   F1':%d  F2':%d  F3':%d   Error:%2.2f percent\n",
		(int)F1, (int)F2, (int)F3, (int)targetF1, (int)targetF2, (int)targetF3, finalError);
	wxPrintf("The error reduced from %2.2f to %2.2f percent.\n", initialError, finalError);
}


// ****************************************************************************
/// \brief Calculates the average frame within the current selection.
///
/// This function calculates the average frame of the current selection by 
/// calculating the arithmetic mean of the sensor data. The average frame is 
/// not returned but instead saved in \ref Data::averageFrame.
// ****************************************************************************

void Data::calcAverageFrame()
{
	if (isValidSelection() == false)
	{
		printf("Invalid selection for calculating the average frame!\n");
		return;
	}

	int firstFrame = (int)(selectionMark_s[0] * (double)Synthesizer::FRAME_RATE_HZ);
	int lastFrame = (int)(selectionMark_s[1] * (double)Synthesizer::FRAME_RATE_HZ);
	int numFrames = lastFrame - firstFrame + 1;
	
	printf("Averaging over %d frames (frame %d to %d).\n", numFrames, firstFrame, lastFrame);

	Synthesizer::SensorFrame *f;
	synthesizer->clearSensorFrame(averageFrame);
	
	int i, k;
	for (i = 0; i < numFrames; i++)
	{
		f = &synthesizer->frameBuffer[firstFrame + i];

		averageFrame.lipOpening_rel += f->lipOpening_rel;
		averageFrame.lipProtrusion_rel += f->lipProtrusion_rel;
		averageFrame.velicOpening_rel += f->velicOpening_rel;

		for (k = 0; k < Synthesizer::MAX_NUM_OPG_SENSORS; k++)
		{
			averageFrame.distance_cm[k] += f->distance_cm[k];
		}
	}
	averageFrame.lipOpening_rel /= (double )numFrames;
	printf("%f\n", averageFrame.lipOpening_rel);
	averageFrame.lipProtrusion_rel /= (double)numFrames;
	printf("%f\n", averageFrame.lipProtrusion_rel);
	averageFrame.velicOpening_rel /= (double)numFrames;
	printf("%f\n", averageFrame.velicOpening_rel);

	for (k = 0; k < Synthesizer::MAX_NUM_OPG_SENSORS; k++)
	{
		averageFrame.distance_cm[k] /= (double)numFrames;
		printf("%f\n", averageFrame.distance_cm[k]);
	}
	printf("Averaging done.\n");
	analysisMode = STATISTICAL;
}

void Data::refreshAreaFunction()
{
	synthesizer->vtAreaFunction->calculateHiResAreaFunction();
	synthesizer->vtAreaFunction->calculateOneDimTubeFunction(synthesizer->tubeFunction, synthesizer->vtAreaFunction->hasRefAreaFunction);
	synthesizer->vtAreaFunction->calculateOneDimTubeFunction(&tlModel->tube, synthesizer->vtAreaFunction->hasRefAreaFunction);
}

// ****************************************************************************
/// Save the spectrum with a resolution of approx. 1 Hz.
// ****************************************************************************

bool Data::saveSpectrum(const char *fileName, ComplexSignal *spectrum)
{
	if (fileName == NULL)
	{
		return false;
	}

	ofstream os(fileName);

	if (!os)
	{
		wxMessageBox(wxString("Could not open ") + wxString(fileName) + wxString(" for writing."),
			wxString("Error!"));
		return false;
	}

	// ****************************************************************
	// Write the samples of the signals.
	// ****************************************************************

	const double AMPLITUDE_SCALE = 1.0;
	double F0 = (double)SAMPLING_RATE / (double)spectrum->N;
	int stepSize = (int)(1.0 / F0);
	if (stepSize < 1)
	{
		stepSize = 1;
	}
	int maxSample = (int)(10000.0 / F0);

	os << "freq_Hz  magnitude  phase_rad" << endl;

	int i;
	double freq_Hz;

	for (i = 0; i < maxSample; i += stepSize)
	{
		freq_Hz = (double)i * F0;
		os << freq_Hz << " "
			<< spectrum->getMagnitude(i) * AMPLITUDE_SCALE << " "
			<< spectrum->getPhase(i) << " "
			<< endl;
	}

	// ****************************************************************
	// Close the file.
	// ****************************************************************

	os.close();

	return true;
}

// ****************************************************************************
/// Returns the index of the echo frame that is under the mark.
/// The returned frame index is always a valid index.
// ****************************************************************************

int Data::getFrameIndexAtMark()
{
  int frameIndex = (int)(mark_s * (double)Synthesizer::FRAME_RATE_HZ);
  if (frameIndex < 0)
  {
    frameIndex = 0;
  }
  if (frameIndex >= Synthesizer::NUM_BUFFER_FRAMES)
  {
    frameIndex = Synthesizer::NUM_BUFFER_FRAMES - 1;
  }

  return frameIndex;
}


// ****************************************************************************
/// \brief Check if current selection markers are valid.
///
/// The selection is only valid if the first selection marker is placed at a
/// lower time index than the second one and both are not below zero.
/// \return true if selection is valid, false if selection is invalid.
// ****************************************************************************

bool Data::isValidSelection()
{
  if ((selectionMark_s[0] < selectionMark_s[1]) &&
    (selectionMark_s[0] >= 0.0) &&
    (selectionMark_s[1] >= 0.0))
  {
    return true;
  }
  else
  {
    return false;
  }
}


// ****************************************************************************
/// Constructor.
// ****************************************************************************

Data::Data()
{
  // Do nothing. Initialization is done in init().
}

// ****************************************************************************
