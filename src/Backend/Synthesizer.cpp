#include "Synthesizer.h"
#include <math.h>
#include "wx/progdlg.h"
//#include "predictAreaFunctionShape.h"
#include "../Frontend/Data.h"

// ****************************************************************************
/// Constructor.
// ****************************************************************************
Synthesizer::Synthesizer(Signal16 *signalBuffer)
{
	currentEpgLayout = EPG32;
	currentOpgLayout = DIST5;

	audioMode = AUDIO_SYNTHESIS;
	dataMode = RAW_DATA;
	modelControlMode = MANUAL;

	mappingOption = REGRESSION;
	kNearestSamples.clear();
	kNearestModels.clear();
	lastBestSampleIdx = 0;

	frameBuffer = new SensorFrame[NUM_BUFFER_FRAMES];
	numRecordedFrames = 0;
	
	audioData.audioBuffer = signalBuffer;
	initAudio();
	
	vtAreaFunction = new OneDimAreaFunction();
	tubeFunction = new Tube();
	vtAreaFunction->calculateOneDimTubeFunction(tubeFunction);

	glottis = new TriangularGlottis();
	double userF0_Hz = 120.0;
	double userPressure_dPa = 10000.0;
	glottis->controlParam[Glottis::FREQUENCY].x = userF0_Hz;
	glottis->controlParam[Glottis::PRESSURE].x = userPressure_dPa;
	glottisParams = new double[glottis->NUM_CONTROL_PARAMS];
	for (int i = 0; i < glottis->NUM_CONTROL_PARAMS; i++)
	{
		glottisParams[i] = glottis->controlParam[i].x;
	}
	baseF0_Hz = 100;
	addF0Wobble = true;

	tdsModel = new TdsModel();
	tdsModel->setTube(tubeFunction);
	newSignal = new double[AUDIO_FRAME_LENGTH];


	outputPressureFilter.createChebyshev(7000.0 / (double)TDS_SAMPLING_RATE, false, 8);
	
	configureSerialPort();
	reset();
	init();
}

// ****************************************************************************
/// Destructor.
// ****************************************************************************

Synthesizer::~Synthesizer()
{
	Pa_Terminate();
  delete[] frameBuffer;
	delete[] newSignal;
}


void Synthesizer::init()
{
	vtAreaFunction->reset();
	vtAreaFunction->calculateOneDimTubeFunction(tubeFunction);
	initialShapesSet = false;
	synthesizeSignalTds(tubeFunction, glottisParams, 1, newSignal);
}

// ****************************************************************************
/// Initialize audio in recording or synthesis mode.
// ****************************************************************************

void Synthesizer::initAudio()
{
	// Close any open stream, just to be safe.
	Pa_CloseStream(&audioStream);

	if (audioMode == AUDIO_SYNTHESIS)
	{
		audioData.audioReadPos = 0;
		audioData.audioWritePos = 0;


		// Initialize portaudio
		Pa_Initialize();

		PaStreamParameters outputParams = { 0 };

		outputParams.device = Pa_GetDefaultOutputDevice();
		outputParams.channelCount = 1;
		outputParams.sampleFormat = paInt16;
		outputParams.suggestedLatency = Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice())->defaultLowOutputLatency;
		outputParams.hostApiSpecificStreamInfo = NULL;

		Pa_OpenStream(&audioStream, NULL, &outputParams, SAMPLING_RATE,
			paFramesPerBufferUnspecified, paNoFlag, &playAudioLoop, &audioData);
	}
	else //audioMode == AUDIO_RECORDING
	{
		audioData.audioReadPos = 0;
		audioData.audioWritePos = 0;

		// Initialize portaudio
		Pa_Initialize();

		PaStreamParameters inputParams = { 0 };

		inputParams.device = Pa_GetDefaultInputDevice();
		inputParams.channelCount = 1;
		inputParams.sampleFormat = paInt16;
		inputParams.suggestedLatency = Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice())->defaultLowOutputLatency;
		inputParams.hostApiSpecificStreamInfo = NULL;

		Pa_OpenStream(&audioStream, &inputParams, NULL, SAMPLING_RATE,
			paFramesPerBufferUnspecified, paNoFlag, &recordAudio, &audioData);
	}
		

}

// ****************************************************************************
/// Calculate three features from the contact sensors
// sensorData: sensor data frame
//
// ****************************************************************************

void Synthesizer::calculateEpgFeatures(SensorFrame *sensorData, double &centerOfGravity, 
	double &meanLaterality, double &sumActivity)
{
	centerOfGravity = 0.0;
	meanLaterality = 0.0;
	sumActivity = 0.0;
	const int NUM_EPG_ROWS = 7;
	const int NUM_EPG_COLS = 6;
	double epgMatrix[NUM_EPG_ROWS][NUM_EPG_COLS] = { 0 };
	epgMatrix[0][2] = sensorData->contact[8];
	epgMatrix[0][3] = sensorData->contact[7];
	epgMatrix[1][1] = sensorData->contact[10];
	epgMatrix[1][2] = sensorData->contact[9];
	epgMatrix[1][3] = sensorData->contact[6];
	epgMatrix[1][4] = sensorData->contact[5];
	epgMatrix[2][1] = sensorData->contact[12];
	epgMatrix[2][2] = sensorData->contact[11];
	epgMatrix[2][3] = sensorData->contact[4];
	epgMatrix[2][4] = sensorData->contact[3];
	epgMatrix[3][0] = sensorData->contact[15];
	epgMatrix[3][1] = sensorData->contact[14];
	epgMatrix[3][2] = sensorData->contact[13];
	epgMatrix[3][3] = sensorData->contact[2];
	epgMatrix[3][4] = sensorData->contact[1];
	epgMatrix[3][5] = sensorData->contact[0];
	epgMatrix[4][0] = sensorData->contact[31];
	epgMatrix[4][1] = sensorData->contact[30];
	epgMatrix[4][2] = sensorData->contact[29];
	epgMatrix[4][3] = sensorData->contact[18];
	epgMatrix[4][4] = sensorData->contact[17];
	epgMatrix[4][5] = sensorData->contact[16];
	epgMatrix[5][1] = sensorData->contact[28];
	epgMatrix[5][2] = sensorData->contact[27];
	epgMatrix[5][3] = sensorData->contact[20];
	epgMatrix[5][4] = sensorData->contact[19];
	epgMatrix[6][0] = sensorData->contact[26];
	epgMatrix[6][1] = sensorData->contact[25];
	epgMatrix[6][2] = sensorData->contact[24];
	epgMatrix[6][3] = sensorData->contact[23];
	epgMatrix[6][4] = sensorData->contact[22];
	epgMatrix[6][5] = sensorData->contact[21];

	// Calculate sum of all activated contacts and normalize to maximum value
	for (int i = 0; i < CONTACT_NUMBER[currentEpgLayout]; i++)
	{
		sumActivity += sensorData->contact[i];
	}
	sumActivity /= CONTACT_NUMBER[currentEpgLayout];

	// Calculate center of gravity(according to the formula in the Articulate Assistant manual)
	double numerator = 0.0;
	double denominator = 7.0 * sumActivity * CONTACT_NUMBER[currentEpgLayout];

	for (int m = 0; m < NUM_EPG_ROWS; m++)
	{
		double sumCols = 0.0;
		for (int c = 0; c < NUM_EPG_COLS; c++)
		{
			sumCols += epgMatrix[m][c];
		}
		numerator += (m - 0.5)*sumCols;
	}
	if (denominator > 0)
	{
		centerOfGravity = 1 - numerator / denominator;
	}
	else
	{
		centerOfGravity = 0.0;
	}
	// Calculate mean Lateral measure(according to formula in the AAA manual)
	numerator = 0.0;
	denominator = 3 * sumActivity * CONTACT_NUMBER[currentEpgLayout];
	for (int n = 0; n < NUM_EPG_COLS; n++)
	{
		double sumRows = 0.0;
		for (int r = 0; r < NUM_EPG_ROWS; r++)
		{
			sumRows += epgMatrix[r][n];
		}
		numerator = abs(n-3.5) * sumRows;
	}
	if (denominator > 0)
	{
		meanLaterality = 1 - numerator / denominator;
	}
	else
	{
		meanLaterality = 0;
	}
}

// ****************************************************************************
// ****************************************************************************
void Synthesizer::reset()
{
	isRunning = false;
	/* Re-open serial port to flush rx buffer */
	openSerialPort();
	numRecordedFrames = 0;
	frameBuffer = new SensorFrame[NUM_BUFFER_FRAMES];
	audioData.audioBuffer->setZero();
	audioData.audioReadPos = 0;
	audioData.audioWritePos = 0;
	audioData.wraparound = false;
	tdsModel->resetMotion();
	glottis->resetMotion();
	glottis->controlParam[Glottis::PRESSURE].x = 10000.0;
	for (int i = 0; i < glottis->NUM_CONTROL_PARAMS; i++)
	{
		glottisParams[i] = glottis->controlParam[i].x;
	}
  outputPressureFilter.resetBuffers();
}

// ****************************************************************************
/// Find the nearest training sample to the current frame.
// ****************************************************************************
int Synthesizer::findNearestTrainingSample(SensorFrame *sensorData)
{
	//double sum = 0.0;
	//double frontness = 0.0;
	//double laterality = 0.0;
	//double asymmetry = 0.0;
	//calculateEpgFeatures(sensorData, asymmetry, frontness, laterality, sum);

	//int sampleIdx = -1;
	//
	//if (kNearestSamples.size() > 0)
	//{
	//	double bestDistance = 0.0;
	//	double distance = 0.0;
	//	// Go through all training samples
	//	for (unsigned int i = 0; i < kNearestSamples.size(); i++)
	//	{
	//		distance = 0.0;
	//		// Calculate the Manhattan distance between the current sensor data frame and the training sample
	//		distance  = fabs(sensorData->lipOpening_rel - kNearestSamples.at(i).at(0));
	//		distance += fabs(sensorData->lipProtrusion_rel - kNearestSamples.at(i).at(1));
	//		for (int j = 0; j < NUM_TONGUE_PALATE_DISTANCES; j++)
	//		{
	//			distance += fabs(sensorData->distance_cm[j] - kNearestSamples.at(i).at(2+j));
	//		}
	//		distance += fabs(asymmetry - kNearestSamples.at(i).at(2 + NUM_TONGUE_PALATE_DISTANCES));
	//		distance += fabs(frontness - kNearestSamples.at(i).at(2 + NUM_TONGUE_PALATE_DISTANCES+1));
	//		distance += fabs(laterality - kNearestSamples.at(i).at(2 + NUM_TONGUE_PALATE_DISTANCES+2));
	//		distance += fabs(sum - kNearestSamples.at(i).at(2 + NUM_TONGUE_PALATE_DISTANCES+3));

	//		if (i > 0)
	//		{
	//			if (distance < bestDistance)
	//			{
	//				bestDistance = distance;
	//				sampleIdx = i;
	//			}
	//		}
	//		else
	//		{
	//			bestDistance = distance;
	//			sampleIdx = lastBestSampleIdx;
	//		}
	//	}
	//	lastBestSampleIdx = sampleIdx;
	//}
	//else 
	//{ 
	//	printf("K Nearest mapping error: No samples loaded!\n"); 
	//}
	//
	//return sampleIdx;
	return 0;
}

// ****************************************************************************
/// Find the set of area function parameters associated with a specific sample.
// ****************************************************************************
void Synthesizer::getParametersFromTrainingSample(int sampleIndex, OneDimAreaFunction::Param *afParams)
{
	
	if (sampleIndex > -1)
	{
		for (unsigned int i = 0; i < OneDimAreaFunction::NUM_AF_PARAMS; i++)
		{
			afParams[i] = kNearestModels.at(sampleIndex).at(i);
		}
	}
	else
	{
		printf("K Nearest mapping error: No associated model found!\n");
	}
}

// ****************************************************************************
/// Calculate the parameters of the vocal tract model using regression models.
// ****************************************************************************

void Synthesizer::calculateAfParametersFromSensorData(SensorFrame *sensorData, OneDimAreaFunction::Param *afParams)
{
	/* Reduce dimensionality by using EPG factors */
	double centerOfGravity = 0.0;
	double meanLaterality = 0.0;
	double sumActivity = 0.0;
	calculateEpgFeatures(sensorData, centerOfGravity, meanLaterality, sumActivity);

	double predictorData[11] = { 0 };
	/* Repack sensor data */
	predictorData[0] = sensorData->lipOpening_rel;
	predictorData[1] = sensorData->lipProtrusion_rel;
	//predictorData[2] = sensorData->velicOpening_rel;
	for (int i = 0; i < DISTANCE_NUMBER[currentOpgLayout]; i++)
	{
		predictorData[2+i] = sensorData->distance_cm[i];
	}
	predictorData[7] = sumActivity;
	predictorData[8] = centerOfGravity;
	predictorData[9] = meanLaterality;

	/* Call MATLAB generated code to predict current shape here */
	double afParamArray[OneDimAreaFunction::NUM_AF_PARAMS] = { 0.0 };
//	predictAreaFunctionShape(predictorData, afParamArray);

	//afParams[OneDimAreaFunction::LLAR_CM].x = 0.2;
	//afParams[OneDimAreaFunction::ALAR_CM2].x = 1.066730;
	//afParams[OneDimAreaFunction::POWLAR].x = 0.107552;
	//afParams[OneDimAreaFunction::XP_CM].x =2.873600;
	//afParams[OneDimAreaFunction::AP_CM2].x = 6.246228;
	//afParams[OneDimAreaFunction::POWP].x = 1.125617;
	//afParams[OneDimAreaFunction::XC_CM].x = 12.583385;
	//afParams[OneDimAreaFunction::AC_CM2].x = 1.465718;
	//afParams[OneDimAreaFunction::POWC].x =3.078920;
	//afParams[OneDimAreaFunction::XA_CM].x = 16.524902;
	//afParams[OneDimAreaFunction::AA_CM2].x = 1.782538;
	//afParams[OneDimAreaFunction::POWA].x = 1.295639;
	//afParams[OneDimAreaFunction::XIN_CM].x = 17.456739;
	//afParams[OneDimAreaFunction::AIN_CM2].x = 2.585221;
	//afParams[OneDimAreaFunction::LVT_CM].x = 18.359900;
	//afParams[OneDimAreaFunction::ALIP_CM2].x = 1.414908;

	afParams[OneDimAreaFunction::LLAR_CM].x = afParamArray[0];
	afParams[OneDimAreaFunction::ALAR_CM2].x = afParamArray[1];
	afParams[OneDimAreaFunction::POWLAR].x = afParamArray[2];
	afParams[OneDimAreaFunction::XP_CM].x = afParamArray[3];
	afParams[OneDimAreaFunction::AP_CM2].x = afParamArray[4]; 
	afParams[OneDimAreaFunction::POWP].x = afParamArray[5];
	afParams[OneDimAreaFunction::XC_CM].x = afParamArray[6];
	afParams[OneDimAreaFunction::AC_CM2].x = afParamArray[7];
	afParams[OneDimAreaFunction::POWC].x = afParamArray[8];
	afParams[OneDimAreaFunction::XA_CM].x = afParamArray[9];
	afParams[OneDimAreaFunction::AA_CM2].x = afParamArray[10];
	afParams[OneDimAreaFunction::POWA].x = afParamArray[11];
	afParams[OneDimAreaFunction::XIN_CM].x = afParamArray[12];
	afParams[OneDimAreaFunction::AIN_CM2].x = afParamArray[13];
	afParams[OneDimAreaFunction::LVT_CM].x = afParamArray[14];
	afParams[OneDimAreaFunction::ALIP_CM2].x = afParamArray[15];

	/* Check parameter bounds */
	//if (afParams[OneDimAreaFunction::LVT_CM].x > 25.0){ afParams[OneDimAreaFunction::LVT_CM].x = 25.0; }
	//if (afParams[OneDimAreaFunction::XIN_CM].x > afParams[OneDimAreaFunction::LVT_CM].x){ afParams[OneDimAreaFunction::XIN_CM].x = afParams[OneDimAreaFunction::LVT_CM].x; }
	//if (afParams[OneDimAreaFunction::XA_CM].x > afParams[OneDimAreaFunction::XIN_CM].x){ afParams[OneDimAreaFunction::XA_CM].x = afParams[OneDimAreaFunction::XIN_CM].x; }
	//if (afParams[OneDimAreaFunction::XC_CM].x > afParams[OneDimAreaFunction::XA_CM].x){ afParams[OneDimAreaFunction::XC_CM].x = afParams[OneDimAreaFunction::XA_CM].x; }
	//if (afParams[OneDimAreaFunction::XP_CM].x > afParams[OneDimAreaFunction::XC_CM].x){ afParams[OneDimAreaFunction::XP_CM].x = afParams[OneDimAreaFunction::XC_CM].x; }
	//if (afParams[OneDimAreaFunction::LLAR_CM].x > afParams[OneDimAreaFunction::XP_CM].x){ afParams[OneDimAreaFunction::LLAR_CM].x = afParams[OneDimAreaFunction::XP_CM].x; }
	//if (afParams[OneDimAreaFunction::LLAR_CM].x < 0.0){ afParams[OneDimAreaFunction::LLAR_CM].x = 0.0; }
	//if (afParams[OneDimAreaFunction::XP_CM].x < afParams[OneDimAreaFunction::LLAR_CM].x){ afParams[OneDimAreaFunction::XP_CM].x = afParams[OneDimAreaFunction::LLAR_CM].x; }
	//if (afParams[OneDimAreaFunction::XC_CM].x < afParams[OneDimAreaFunction::XP_CM].x){ afParams[OneDimAreaFunction::XC_CM].x = afParams[OneDimAreaFunction::XP_CM].x; }
	//if (afParams[OneDimAreaFunction::XA_CM].x < afParams[OneDimAreaFunction::XC_CM].x){ afParams[OneDimAreaFunction::XA_CM].x = afParams[OneDimAreaFunction::XC_CM].x; }
	//if (afParams[OneDimAreaFunction::XIN_CM].x < afParams[OneDimAreaFunction::XA_CM].x){ afParams[OneDimAreaFunction::XIN_CM].x = afParams[OneDimAreaFunction::XA_CM].x; }
	//if (afParams[OneDimAreaFunction::LVT_CM].x < afParams[OneDimAreaFunction::XIN_CM].x){ afParams[OneDimAreaFunction::LVT_CM].x = afParams[OneDimAreaFunction::XIN_CM].x; }
	//if (afParams[OneDimAreaFunction::ALAR_CM2].x < -2.0){ afParams[OneDimAreaFunction::ALAR_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::AP_CM2].x < -2.0){ afParams[OneDimAreaFunction::AP_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::AC_CM2].x < -2.0){ afParams[OneDimAreaFunction::AC_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::AA_CM2].x < -2.0){ afParams[OneDimAreaFunction::AA_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::AIN_CM2].x < -2.0){ afParams[OneDimAreaFunction::AIN_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::ALIP_CM2].x < -2.0){ afParams[OneDimAreaFunction::ALIP_CM2].x = -2.0; }
	//if (afParams[OneDimAreaFunction::POWLAR].x < 0.1){ afParams[OneDimAreaFunction::POWLAR].x = 0.1; }
	//if (afParams[OneDimAreaFunction::POWP].x < 0.1){ afParams[OneDimAreaFunction::POWP].x = 0.1; }
	//if (afParams[OneDimAreaFunction::POWC].x < 0.1){ afParams[OneDimAreaFunction::POWC].x = 0.1; }
	//if (afParams[OneDimAreaFunction::POWA].x < 0.1){ afParams[OneDimAreaFunction::POWA].x = 0.1; }
	//if (afParams[OneDimAreaFunction::POWLAR].x > 10.0){ afParams[OneDimAreaFunction::POWLAR].x = 10.0; }
	//if (afParams[OneDimAreaFunction::POWP].x > 10.0){ afParams[OneDimAreaFunction::POWP].x = 10.0; }
	//if (afParams[OneDimAreaFunction::POWC].x > 10.0){ afParams[OneDimAreaFunction::POWC].x = 10.0; }
	//if (afParams[OneDimAreaFunction::POWA].x > 10.0){ afParams[OneDimAreaFunction::POWA].x = 10.0; }
}

// ****************************************************************************
// ****************************************************************************

bool Synthesizer::toggleRun()
{
	if (isRunning)
	{
		stopRealTimeAudio();
		isRunning = false;
	}
	else
	{
		startRealTimeAudio();
		isRunning = true;
	}	

	return isRunning;
}

void Synthesizer::clearSensorFrame(SensorFrame & f)
{	
	int i;
	f.frameIndex = 0.0;
	f.lipOpening_rel = 0.0;
	f.lipProtrusion_rel = 0.0;
	f.velicOpening_rel = 0.0;
	for (i = 0; i < MAX_NUM_OPG_SENSORS; i++)
	{
		f.distance_cm[i] = 0.0;
	}
	for (i = 0; i < MAX_NUM_CONTACT_SENSORS; i++)
	{
		f.contact[i] = false;
	}
	f.f0_Hz = 100.0;
	f.voicing = true;
	f.breathiness_rel = 0.0;
	f.pressure_Pa = 800;	 
}

void Synthesizer::startAudioStream()
{
	if (audioMode == AUDIO_SYNTHESIS) 
	{
		wxPrintf("Starting audio playback...\n");
	}
	else // audioMode == AUDIO_RECORDING
	{
		wxPrintf("Starting audio recording...\n");
	}

	Pa_StartStream(audioStream);
}

void Synthesizer::startRealTimeAudio()
{	
	stopAudioStream();
	initAudio();
	reset();
	//wxPrintf("Starting %i Hz data update timer...\n", AUDIO_UPDATE_RATE_HZ);
	//this->Start(1 / AUDIO_UPDATE_RATE_HZ * 1000);
	startAudioStream();
}

void Synthesizer::stopRealTimeAudio()
{
	//wxPrintf("Stopping %i Hz data update timer.\n", AUDIO_UPDATE_RATE_HZ);
	//this->Stop();
	stopAudioStream();
	
}

void Synthesizer::stopAudioStream()
{
	Pa_StopStream(audioStream);
	if (audioMode == AUDIO_SYNTHESIS)
	{
		wxPrintf("Audio playback stopped.\n");
	}
	else // audioMode == AUDIO_RECORDING
	{
		wxPrintf("Audio recording stopped.\n");
	}

}

// ****************************************************************************
/// Synthesize an incremental part of the signal.
/// numNewSamples is the number of new samples in the final audio signal 
/// (22050 Hz).
/// The vector audio receives the new numNewSamples samples generated here,
/// if it is != NULL. The range of samples values is [-1; +1].
// ****************************************************************************

void Synthesizer::synthesizeSignalTds(Tube *newTube, double *newGlottisParams, int numNewSamples, double *newSignal)
{
	
	int i, k;

	int numGlottisParams = glottis->controlParam.size();

	if (initialShapesSet == false)
	{
		prevTube = *newTube;
		for (i = 0; i < numGlottisParams; i++)
		{
			prevGlottisParams[i] = newGlottisParams[i];
		}

		initialShapesSet = true;
		return;
	}

	// ****************************************************************
	// The actual number of samples to be calculated is twice the
	// number of samples in the final signal (44100 Hz vs. 22050 Hz).
	// ****************************************************************

	numNewSamples *= TDS_SR_FACTOR;

	if (numNewSamples < 1)
	{
		numNewSamples = 1;
	}

	double ratio, ratio1;
	double length_cm[Tube::NUM_GLOTTIS_SECTIONS];
	double area_cm2[Tube::NUM_GLOTTIS_SECTIONS];
	bool filtering;
	double pressure_dPa[4];
	double mouthFlow_cm3_s;

	// ****************************************************************
	// Run through all new samples.
	// ****************************************************************

	for (i = 0; i < numNewSamples; i++)
	{
		ratio = (double)i / (double)numNewSamples;
		ratio1 = 1.0 - ratio;

		// ****************************************************************
		// Interpolate the tube.
		// ****************************************************************
		
		tube.interpolate(&prevTube, newTube, ratio);

		// ****************************************************************
		// Interpolate the glottis geometry.
		// ****************************************************************

		for (k = 0; k < numGlottisParams; k++)
		{
			glottis->controlParam[k].x = ratio1*prevGlottisParams[k] + ratio*newGlottisParams[k];
		}
		glottis->calcGeometry();
		glottis->getTubeData(length_cm, area_cm2);

		tube.setGlottisGeometry(length_cm, area_cm2);
		tube.setAspirationStrength(glottis->getAspirationStrength_dB());

		// ****************************************************************
		// Do the acoustic simulation.
		// ****************************************************************

		if (tdsModel->getSampleIndex() == 0)
		{
			filtering = false;
		}
		else
		{
			filtering = true;
		}

		tdsModel->setTube(&tube, filtering);
		tdsModel->setFlowSource(0.0, -1);
		tdsModel->setPressureSource(glottis->controlParam[Glottis::PRESSURE].x, Tube::FIRST_TRACHEA_SECTION);

		// Get the four relevant pressure values for the glottis model:
		// subglottal, lower glottis, upper glottis, supraglottal.

		pressure_dPa[0] = tdsModel->getSectionPressure(Tube::LAST_TRACHEA_SECTION);
		pressure_dPa[1] = tdsModel->getSectionPressure(Tube::LOWER_GLOTTIS_SECTION);
		pressure_dPa[2] = tdsModel->getSectionPressure(Tube::UPPER_GLOTTIS_SECTION);
		pressure_dPa[3] = tdsModel->getSectionPressure(Tube::FIRST_PHARYNX_SECTION);
		
		
		// Increment the time/sample number
		glottis->incTime(1.0 / (double)TDS_SAMPLING_RATE, pressure_dPa);

		mouthFlow_cm3_s = tdsModel->proceedTimeStep();      // The next flow sample


		int pos = tdsModel->getSampleIndex();
		k = pos & TDS_BUFFER_MASK;
		outputFlow[k] = mouthFlow_cm3_s;
		outputPressure[k] = (outputFlow[k] - outputFlow[(k - 1) & TDS_BUFFER_MASK]) / tdsModel->timeStep;
		filteredOutputPressure[k] = outputPressureFilter.getOutputSample(outputPressure[k]);

		if ((pos % TDS_SR_FACTOR) == 0)
		{
			double s = filteredOutputPressure[k] * 0.004;
			if (newSignal != NULL)
			{
				newSignal[i / TDS_SR_FACTOR] = s / SHRT_MAX;
			}
		}

	}

	// ****************************************************************

	prevTube = *newTube;
	for (i = 0; i < numGlottisParams; i++)
	{
		prevGlottisParams[i] = newGlottisParams[i];
	}

}

// ****************************************************************************
/// Synthesize a long or short vowel using the frequency domain synthesis.
// ****************************************************************************

void Synthesizer::synthesizeVowelLf(TlModel *tlModel, LfPulse &lfPulse, int startPos, bool isLongVowel)
{
	reset();
	const int NUM_F0_NODES = 4;
	const int NUM_AMP_NODES = 4;
	const int BUFFER_LENGTH = 2048;
	const int BUFFER_MASK = 2047;

	TimeFunction ampTimeFunction;
	TimeFunction f0TimeFunction;
	double duration_ms;
	Signal singlePulse;
	Signal pulseSignal(BUFFER_LENGTH);
	Signal pressureSignal(BUFFER_LENGTH);
	int i, k;
	int nextPulsePos = 10;   // Get the first pulse shape at sample number 10
	int pulseLength;
	double t_s, t_ms;
	double sum;
	double filteredValue;

	// Memorize the pulse params to restore them at the end of the function
	LfPulse origLfPulse = lfPulse;

	// ****************************************************************
	// Init the time functions for F0 and glottal pulse amplitude.
	// ****************************************************************

	if (isLongVowel)
	{
		duration_ms = 650.0;
		double max = lfPulse.F0;

		TimeFunction::Node f0[NUM_F0_NODES] =
		{
			{ 0.0, 0.9*max },
			{ 300.0, 1.00*max },
			{ 450.0, 0.8*max },
			{ 600.0, 0.7*max }
		};

		TimeFunction::Node amp[NUM_AMP_NODES] =
		{
			{ 0.0, 0.0 },
			{ 40.0, 500.0 },
			{ 400.0, 450.0 },
			{ 600.0, 0.0 }
		};

		ampTimeFunction.setNodes(amp, NUM_AMP_NODES);
		f0TimeFunction.setNodes(f0, NUM_F0_NODES);
	}
	else
	{
		duration_ms = 350.0;
		double max = lfPulse.F0;

		TimeFunction::Node f0[NUM_F0_NODES] =
		{
			{ 0.0, 0.9*max },
			{ 125.0, 1.0*max },
			{ 126.0, 1.0*max },
			{ 300.0, 0.82*max }
		};
		TimeFunction::Node amp[NUM_AMP_NODES] =
		{
			{ 0.0, 0.0 },
			{ 20.0, 500.0 },
			{ 200.0, 450.0 },
			{ 300.0, 0.0 }
		};

		ampTimeFunction.setNodes(amp, NUM_AMP_NODES);
		f0TimeFunction.setNodes(f0, NUM_F0_NODES);
	}

	// The length in samples

	int length = (int)((duration_ms / 1000.0)*(double)SAMPLING_RATE);

	// Init the low-pass filter

	const double NUM_LOWPASS_POLES = 6;
	IirFilter filter;
	filter.createChebyshev((double)SYNTHETIC_SPEECH_BANDWIDTH_HZ / (double)SAMPLING_RATE, false, (int)NUM_LOWPASS_POLES);

	// ****************************************************************
	// Calc. the vocal tract impulse response.
	// ****************************************************************

	const int IMPULSE_RESPONSE_EXPONENT = 9;
	const int IMPULSE_RESPONSE_LENGTH = 1 << IMPULSE_RESPONSE_EXPONENT;
	Signal impulseResponse(IMPULSE_RESPONSE_LENGTH);

	tlModel->getImpulseResponse(&impulseResponse, IMPULSE_RESPONSE_EXPONENT);

	// Reduce amplitude with increasing F0.
	impulseResponse *= 80.0 / lfPulse.F0;

	// ****************************************************************
	// Calc. the speech signal samples.
	// ****************************************************************

	for (i = 0; i < length; i++)
	{
		t_s = (double)i / (double)SAMPLING_RATE;
		t_ms = t_s*1000.0;

		// **************************************************************
		// Is a new glottal pulse starting?
		// **************************************************************

		if (i == nextPulsePos)
		{
			// Get the pulse amplitude and F0.

			lfPulse.AMP = ampTimeFunction.getValue(t_ms);
			lfPulse.F0 = f0TimeFunction.getValue(t_ms);

			// Simulate "flutter".

			lfPulse.F0 += 0.5*(lfPulse.F0 / 100.0)*(sin(2.0*M_PI*12.7*t_s) + sin(2.0*M_PI*7.1*t_s) + sin(2.0*M_PI*4.7*t_s));

			// Get and set the new glottal pulse.

			pulseLength = (int)((double)SAMPLING_RATE / lfPulse.F0);
			lfPulse.getPulse(singlePulse, pulseLength, false);
			for (k = 0; k < pulseLength; k++)
			{
				pulseSignal.x[(i + k) & BUFFER_MASK] = singlePulse.getValue(k);
			}

			nextPulsePos += pulseLength;
		}

		// **************************************************************
		// Do the convolution.
		// **************************************************************

		sum = 0.0;
		for (k = 0; k < IMPULSE_RESPONSE_LENGTH; k++)
		{
			sum += impulseResponse.x[k] * pulseSignal.x[(i - k) & BUFFER_MASK];
		}

		pressureSignal.x[i & BUFFER_MASK] = sum;

		// Set the final sample.

		filteredValue = 2000.0 * filter.getOutputSample(pressureSignal.getValue(i));
		audioData.audioBuffer->setValue(startPos + i, (short)filteredValue);
	}
	audioData.audioWritePos = startPos + i;

	// Restore the pulse params
	lfPulse = origLfPulse;
}

// ****************************************************************************
/// Configures a COM port for communication with EOS.
// ****************************************************************************
bool Synthesizer::configureSerialPort(wxString portName, int baudRate)
{
//	strcpy(comSettings.portName, portName);
//	comSettings.baudRate = baudRate;
	return true;
}

// ****************************************************************************
/// Opens a COM port for communication with EOS hardware with current settings.
// ****************************************************************************

bool Synthesizer::openSerialPort()
{
//	if (serialPort.open(comSettings.portName, comSettings.baudRate))
//	{
//		wxPrintf("Serial port for palate data successfully opened!\n");
//		return true;
//	}
//	else
//	{
//		wxPrintf("Failed to open serial port for palate data!\n");
		return false;
//	}
}

// ****************************************************************************
/// Closes the current COM port.
// ****************************************************************************

bool Synthesizer::closeSerialPort()
{
//	if (serialPort.isOpen())
//	{
//		serialPort.close();
//		wxPrintf("Serial port for palate data successfully closed!\n");
//		return true;
//	}
//	else
//	{
//		wxPrintf("Failed to close serial port for palate data!\n");
		return false;
//	}
}

// ****************************************************************************
// ****************************************************************************

void Synthesizer::sensorDataToTube(SensorFrame *sensorData, Tube *tube)
{
	Data *data = Data::getInstance();
	OneDimAreaFunction::Param afParams[OneDimAreaFunction::NUM_AF_PARAMS];
	if (isUsingSensorData)
	{
		switch (mappingOption)
		{
			int sampleIdx;
		case REGRESSION:
			calculateAfParametersFromSensorData(sensorData, afParams);
			break;
		case K_NEAREST:
			sampleIdx = findNearestTrainingSample(sensorData);
			getParametersFromTrainingSample(sampleIdx, afParams);
			break;
			default:
				break;
		}
		vtAreaFunction->setParameters(afParams);
		vtAreaFunction->calculateOneDimTubeFunction(tube);
		//tube->setVelumOpening(sensorData->velicOpening_rel * 2.0);
	}
	else
	{
		vtAreaFunction->calculateOneDimTubeFunction(tube);
	}				
	data->refreshAreaFunction();
	vtAreaFunction->getParameters(afParams);
	tube->teethPosition_cm = afParams[OneDimAreaFunction::XIN_CM].x;
}


// ****************************************************************************
// ****************************************************************************

void Synthesizer::sensorDataToGlottisParams(SensorFrame *sensorData, double *glottisParams)
{
	for (int i = 0; i < glottis->NUM_CONTROL_PARAMS; i++)
	{
		glottisParams[i] = glottis->controlParam[i].x;
	}

	double f0_Hz;
	//if (isUsingSensorData)
	//{
	//	f0_Hz = sensorData->f0_Hz;
	//	glottisParams[Glottis::PRESSURE] = sensorData->pressure_Pa * 10.0;
	//}
	//else
	//{
		f0_Hz = baseF0_Hz;
		glottisParams[Glottis::PRESSURE] = 8000.0;
	//}
	if (addF0Wobble)
	{
		double t = (double)audioData.audioReadPos / SAMPLING_RATE;
		double microVariationFactor = (1 + 0.25 / 50 * 1 / 100 *
			(sin(2 * M_PI * 12.7 * t)
			+ sin(2 * M_PI * 7.1 * t)
			+ sin(2 * M_PI * 4.7 * t)));
		double macroVariationFactor = (1 + 0.25 / 50 *
			(sin(2 * M_PI * 2.7 * t)
			+ sin(2 * M_PI * 0.8 * t)
			+ sin(2 * M_PI * 1.3 * t)));
		f0_Hz = f0_Hz * microVariationFactor * macroVariationFactor;
	}	
	glottisParams[Glottis::FREQUENCY] = f0_Hz;
}


// ****************************************************************************
/// This function is called to get synthesize a new speech segment.
// ****************************************************************************

void Synthesizer::synthesizeSegment()
{
	if (audioMode == AUDIO_SYNTHESIS)
	{	
		
		if (
			(
			audioData.wraparound == false // buffer wrap around has not occurred
			&& 
			(signed int) audioData.audioWritePos - (signed int) audioData.audioReadPos < 3 * AUDIO_FRAME_LENGTH // The output is not lagging too much
			)
			|| 
			(
			audioData.wraparound == true // buffer wrap around has occurred
			&& 
			(NUM_BUFFER_SAMPLES - (signed int)audioData.audioReadPos + (signed int) audioData.audioWritePos) <  3 * AUDIO_FRAME_LENGTH  // The output is not lagging too much
			)
			)
		{
		fetchSensorData();

		sensorDataToTube(&sensorFrame, tubeFunction);

		sensorDataToGlottisParams(&sensorFrame, glottisParams);

		synthesizeSignalTds(tubeFunction, glottisParams, AUDIO_FRAME_LENGTH, newSignal);

			for (int i = 0; i < AUDIO_FRAME_LENGTH; i++)
			{
				audioData.audioBuffer->setValue(audioData.audioWritePos, (double)newSignal[i] * SHRT_MAX);
				if (newSignal[i] > 1.0)
				{
					audioData.audioBuffer->setValue(audioData.audioWritePos, SHRT_MAX);
				}
				if (newSignal[i] < -1.0)
				{
					audioData.audioBuffer->setValue(audioData.audioWritePos, SHRT_MIN);
				}
				audioData.audioWritePos++;
				/* Check buffer overflow */
				if (audioData.audioWritePos >= NUM_BUFFER_SAMPLES)
				{
					audioData.audioWritePos = 0;
					audioData.wraparound = true;
				}
			}		

		
		}	
		else
		{
			// Do nothing and wait for playback to catch up with synthesis
		}
	}
	else // audioMode == AUDIO_RECORDING
	{
		fetchSensorData();
	}
}

// ****************************************************************************
/// Set and get the audio mode to use.
// ****************************************************************************

int Synthesizer::getAudioMode()
{
	return this->audioMode;
}

int Synthesizer::setAudioMode(AudioMode newAudioMode)
{
	this->audioMode = newAudioMode;
	return this->audioMode;
}

// ****************************************************************************
/// Set the mapping technique to use.
// ****************************************************************************
int Synthesizer::setMappingOption(MappingOption mappingOption)
{
	this->mappingOption = mappingOption;
	
	return this->mappingOption;
}


// ****************************************************************************
// ****************************************************************************

unsigned char bigBuffer[64000];
// ****************************************************************************
/// Get a frame of sensor data from COM port.
// ****************************************************************************
void Synthesizer::fetchSensorData()
{
//	SerialPort *port = &serialPort;
//
//	if (port->isOpen() == false)
//	{
//		//printf("Sensor port is not open.\n");
//		return;
//	}
//
//	const int FRAME_LENGTH = 34;  // Frame length without the 2 sync bytes
//	enum State { WAIT_START_BYTE_0, WAIT_START_BYTE_1, WAIT_DATA_BYTE };
//	static State state = WAIT_START_BYTE_0;
//	const int START_BYTE_0 = 0x5A;
//	const int START_BYTE_1 = 0xA5;
//	unsigned char ch;
//
//	// ****************************************************************
//	// Interpret the received bytes one by one and try to extract a
//	// frame.
//	// ****************************************************************
//
//	int newBytes = port->readData((char*)bigBuffer, port->bytesAvailable());
//
//	for (int k = 0; k < newBytes; k++)
//	{
//		ch = bigBuffer[k];
//
//		switch (state)
//		{
//			// ************************************************************
//
//		case WAIT_START_BYTE_0:
//			if ((unsigned char)ch == START_BYTE_0)
//			{
//				state = WAIT_START_BYTE_1;
//			}
//			break;
//
//			// ************************************************************
//
//		case WAIT_START_BYTE_1:
//			if ((unsigned char)ch == START_BYTE_1)
//			{
//				numRxBytes = 0;
//				state = WAIT_DATA_BYTE;
//			}
//			else
//			if (ch != START_BYTE_1)
//			{
//				state = WAIT_START_BYTE_0;
//			}
//			break;
//
//			// ************************************************************
//
//		case WAIT_DATA_BYTE:
//			if (numRxBytes < MAX_BUFFER_BYTES)
//			{
//				rxBuffer[numRxBytes++] = ch;
//				if (numRxBytes >= FRAME_LENGTH)
//				{
//					// Check the check sum
//					int i;
//					int checkSum = 0;
//
//					for (i = 0; i < FRAME_LENGTH - 1; i++)
//					{
//						checkSum += rxBuffer[i];
//					}
//
//					if ((unsigned char)checkSum == rxBuffer[FRAME_LENGTH - 1])
//					{
//						assignSensorRxData();
//					}
//					else
//					{
//						wxPrintf("Mismatch of checksum for sensor data frame!\n");
//					}
//
//					state = WAIT_START_BYTE_0;
//				}
//			}
//			else
//			{
//				wxPrintf("Buffer overflow!\n");
//				state = WAIT_START_BYTE_0;
//			}
//			break;
//		}
//	}
}

// ****************************************************************************
/// Assign the received sensor data to their respective containers.
// ****************************************************************************
void Synthesizer::assignSensorRxData()
{
	unsigned char *buffer = rxBuffer;

	// Get the frame index
	sensorFrame.frameIndex = ((unsigned int)buffer[1] << 8) + buffer[0];

	/* The data can be pre-processed or raw */ 
	if (dataMode == PROCESSED_DATA) // data is already processed (distances and lip factors)
	{
		
		// Gather all other sensor data and scale to the correct dimensions
		int pos = 2;    // Position of next byte to interpret.
		int value = 0;

		// Lip opening
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.lipOpening_rel = (double)value * 0.001; // Data is transmitted as an integer between 0 and 1000

		// Lip protrusion
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.lipProtrusion_rel = (double)value * 0.001; // Data is transmitted as an integer between 0 and 1000

		// velic opening 
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.velicOpening_rel = (double)value * 0.001; // Data is transmitted as an integer between 0 and 1000

		// Tongue-palate distances
		for (int i = 0; i < MAX_NUM_OPG_SENSORS; i++)
		{
			value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
			pos += 2;
			sensorFrame.distance_cm[i] = (double)value * 0.001; // Floating point data was multiplied by 1000 and transmitted as an integer
		}

		// EPG data is a 64 bit wide bit field (Each unsigned int holds 26 bits of contact data)
		unsigned int values[2];
		unsigned int bitPos = 0;
		for (int i = 0; i < 2; i++)
		{
			values[i] = (((unsigned int)buffer[pos + 3]) << 24)
				+ (((unsigned int)buffer[pos + 2]) << 16)
				+ (((unsigned int)buffer[pos + 1]) << 8)
				+ ((unsigned int)buffer[pos]);

			for (int j = 0; j < MAX_NUM_CONTACT_SENSORS / 2; j++)
			{
				sensorFrame.contact[bitPos] = (unsigned char)((values[i] >> j) & 0x1);
				bitPos++;
			}
			pos += 4;
		}

		// F0 data
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.f0_Hz = (double)value * 0.01; // Floating point data was multiplied by 100 and transmitted as an integer

		// Voicing
		#pragma warning(suppress: 4800) // Performance warning because integer is converted to bool
		sensorFrame.voicing = (bool)buffer[pos];
		pos++;

		// Breathiness
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.breathiness_rel = (double)value * 0.001; // Data is transmitted as an integer between 0 and 1000

		// Pressure
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.pressure_Pa = (double)value * 0.1; // Floating point data was multiplied by ten and transmitted as an integer

	}
	else if (dataMode == RAW_DATA) //data are raw ADC values
	{
		// Gather all other sensor data and scale to the correct dimensions
		int pos = 2;    // Position of next byte to interpret.
		int value = 0;

		// Optical sensor 8
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.lipOpening_rel = (double)value; // Data is transmitted as an integer between 0 and 4095

		// Optical sensor 6
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.lipProtrusion_rel = (double)value; // Data is transmitted as an integer between 0 and 4095

		// Optical sensor 7
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.velicOpening_rel = (double)value; // Data is transmitted as an integer between 0 and 4095

		// Tongue-palate distance sensor values (optical sensors 1 to 5)
		for (int i = 0; i < MAX_NUM_OPG_SENSORS; i++)
		{
			value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
			pos += 2;
			sensorFrame.distance_cm[i] = (double)value; // Data is transmitted as an integer between 0 and 4095
		}

		// EPG data is a 64 bit wide bit field (Each unsigned int holds 26 bits of contact data)
		unsigned int values[2];
		unsigned int bitPos = 0;
		for (int i = 0; i < 2; i++)
		{
			values[i] = (((unsigned int)buffer[pos + 3]) << 24)
				+ (((unsigned int)buffer[pos + 2]) << 16)
				+ (((unsigned int)buffer[pos + 1]) << 8)
				+ ((unsigned int)buffer[pos]);
			
			for (int j = 0; j < MAX_NUM_CONTACT_SENSORS / 2; j++)
			{
				sensorFrame.contact[bitPos] = (unsigned char)((values[i] >> j) & 0x1);
				bitPos++;
			}				
			pos += 4;
		}
		// F0 data
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.f0_Hz = (double)value * 0.1; // Floating point data was multiplied by 10 and transmitted as an integer

		// Voicing
		#pragma warning(suppress: 4800) // Performance warning because integer is converted to bool
		sensorFrame.voicing = (bool)buffer[pos];
		pos++;

		// Breathiness
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.breathiness_rel = (double)value * 0.001; // Data is transmitted as an integer between 0 and 1000

		 // Pressure
		value = (int)((unsigned int)buffer[pos + 1] << 8) + buffer[pos];
		pos += 2;
		sensorFrame.pressure_Pa = (double)value * 0.1; // Floating point data was multiplied by ten and transmitted as an integer
	}

	///* Low-pass filter distance and lip data */
	//const double COEFF_FILTER = 0.6;


	//if (numRecordedFrames > 0)
	//{
	//	for (int i = 0; i < NUM_TONGUE_PALATE_DISTANCES; i++)
	//	{
	//		sensorFrame.distance_cm[i] = COEFF_FILTER * sensorFrame.distance_cm[i] + (1 - COEFF_FILTER) * frameBuffer[numRecordedFrames - 1].distance_cm[i];
	//	}
	//	sensorFrame.lipOpening_rel = COEFF_FILTER * sensorFrame.lipOpening_rel + (1 - COEFF_FILTER) * frameBuffer[numRecordedFrames - 1].lipOpening_rel;
	//	sensorFrame.lipProtrusion_rel = COEFF_FILTER * sensorFrame.lipProtrusion_rel + (1 - COEFF_FILTER) * frameBuffer[numRecordedFrames - 1].lipProtrusion_rel;
	//}
	

	if (numRecordedFrames < NUM_BUFFER_FRAMES)
	{
		frameBuffer[numRecordedFrames] = sensorFrame;
		numRecordedFrames++;
	}
}

// ****************************************************************************
/// Interpolate between two sets of parameters to get the parameters at tX.
// ****************************************************************************
void Synthesizer::interpolateParameters(OneDimAreaFunction::Param *p0,
	OneDimAreaFunction::Param *p1, OneDimAreaFunction::Param *pX,
	double t0, double t1, double tX)
{
	for (int i = 0; i < OneDimAreaFunction::NUM_AF_PARAMS; i++)
	{
		pX[i].x = (p1[i].x - p0[i].x) / 2 * cos((t1 - tX) / (t1 - t0) * M_PI) + (p1[i].x + p0[i].x) / 2;
	}
}

// ****************************************************************************
/// Synthesize the samples corresponding to a sequence of target shapes.
// ****************************************************************************
void Synthesizer::playTargetSequence(OneDimAreaFunction::Shape *targetShape,
	double *stationary_s, double *transition_s)
{
	stopAudioStream();

	this->targetShape = targetShape;
	this->stationary_s = stationary_s;
	this->transition_s = transition_s;
	reset();
	init();	

	double f0_Hz[4] = { 100, 115, 105, 80 };
	

	sensorDataToGlottisParams(&sensorFrame, glottisParams);
	double steadyStateLungPressure = glottisParams[Glottis::PRESSURE];

	OneDimAreaFunction::Param currentParams[OneDimAreaFunction::NUM_AF_PARAMS];
	
	double boundary_s[7] = {
		stationary_s[0],
		stationary_s[0] + transition_s[0],
		stationary_s[0] + transition_s[0] + stationary_s[1],
		stationary_s[0] + transition_s[0] + stationary_s[1] + transition_s[1],
		stationary_s[0] + transition_s[0] + stationary_s[1] + transition_s[1] + stationary_s[2],
		stationary_s[0] + transition_s[0] + stationary_s[1] + transition_s[1] + stationary_s[2] + transition_s[2],
		stationary_s[0] + transition_s[0] + stationary_s[1] + transition_s[1] + stationary_s[2] + transition_s[2] + stationary_s[3],
	};
	double totalTime_s = boundary_s[6];
	int numSamples = SAMPLING_RATE * totalTime_s;
	
	for (int i = 0; i < numSamples; i++)
	{
		/* Fade in lung pressure */
		if (i < (double) 0.1 * (double) SAMPLING_RATE)
		{
			/* First 50 ms of silence to avoid initial click */
			if (i < (double) 0.05 * (double)SAMPLING_RATE)
			{
				glottisParams[Glottis::PRESSURE] = 0.0;
			}
			/* Next 50 ms cosine fade-in */
			else
			{
				glottisParams[Glottis::PRESSURE] = steadyStateLungPressure / 2 * cos((0.1*SAMPLING_RATE - i) / (0.05*SAMPLING_RATE) * M_PI) + steadyStateLungPressure / 2;
			}
			
		}
		/* Create more natural F0 contour */
		if (i < (double) boundary_s[1] * SAMPLING_RATE)
		{
			glottisParams[Glottis::FREQUENCY] = (f0_Hz[0] + f0_Hz[1]) / 2 + (f0_Hz[1] - f0_Hz[0]) / 2 * cos((boundary_s[1] * SAMPLING_RATE-i) / (boundary_s[1] * SAMPLING_RATE) * M_PI);
		}
		else if (i < (double) boundary_s[3] * SAMPLING_RATE)
		{
			glottisParams[Glottis::FREQUENCY] = (f0_Hz[2] + f0_Hz[1]) / 2 + (f0_Hz[2] - f0_Hz[1]) / 2 * cos((boundary_s[3] * SAMPLING_RATE-i) / ((boundary_s[3] - boundary_s[1])*SAMPLING_RATE) * M_PI);
		}
		else
		{
			glottisParams[Glottis::FREQUENCY] = (f0_Hz[3] + f0_Hz[2]) / 2 + (f0_Hz[3] - f0_Hz[2]) / 2 * cos((boundary_s[6] * SAMPLING_RATE-i) / ((boundary_s[6] - boundary_s[3])*SAMPLING_RATE) * M_PI);
		}


		if (i <= boundary_s[0] * SAMPLING_RATE)
		{
			for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
			{
				currentParams[j] = targetShape[0].param[j];
			}
		}
		else if (i <= boundary_s[1] * SAMPLING_RATE)
		{
			interpolateParameters(targetShape[0].param, targetShape[1].param, &currentParams[0],
				boundary_s[0] * SAMPLING_RATE, boundary_s[1] * SAMPLING_RATE, i);
		}
		else if (i <= boundary_s[2] * SAMPLING_RATE)
		{
			for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
			{
				currentParams[j] = targetShape[1].param[j];
			}
		}
		else if (i <= boundary_s[3] * SAMPLING_RATE)
		{
			interpolateParameters(targetShape[1].param, targetShape[2].param, &currentParams[0],
				boundary_s[2] * SAMPLING_RATE, boundary_s[3] * SAMPLING_RATE, i);
		}
		else if (i <= boundary_s[4] * SAMPLING_RATE)
		{
			for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
			{
				currentParams[j] = targetShape[2].param[j];
			}
		}
		else if (i <= boundary_s[5] * SAMPLING_RATE)
		{
			interpolateParameters(targetShape[2].param, targetShape[3].param, &currentParams[0],
				boundary_s[4] * SAMPLING_RATE, boundary_s[5] * SAMPLING_RATE, i);
		}
		else if (i <= boundary_s[6] * SAMPLING_RATE)
		{
			for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
			{
				currentParams[j] = targetShape[3].param[j];
			}
		}

		/* Fade out lung pressure */
		if (i >(double) (totalTime_s - 0.1) * (double) SAMPLING_RATE)
		{
			glottisParams[Glottis::PRESSURE] = -steadyStateLungPressure / 2 * cos((totalTime_s * SAMPLING_RATE - i) / (totalTime_s - 0.1*SAMPLING_RATE) * M_PI) + steadyStateLungPressure / 2;
		}

		
		vtAreaFunction->setParameters(currentParams);
		vtAreaFunction->calculateOneDimTubeFunction(tubeFunction);
		tubeFunction->teethPosition_cm = currentParams[OneDimAreaFunction::XIN_CM].x;
		synthesizeSignalTds(tubeFunction, glottisParams, 1, newSignal);
		audioData.audioBuffer->setValue(audioData.audioWritePos, (double)newSignal[0] * SHRT_MAX);
		audioData.audioWritePos++;
	}
	
	startAudioStream();
}



// ****************************************************************************
/// Portaudio callback function for continuous closed-loop playback.
// ****************************************************************************
int  Synthesizer::playAudioLoop(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo*,
	PaStreamCallbackFlags,
	void *userData)
{
	(void)inputBuffer; /* Prevent unused variable warning. */
	signed short* out = (signed short*) outputBuffer;
	paData *data = (paData*) userData;

	if ((data->wraparound == false) && data->audioReadPos < data->audioWritePos - framesPerBuffer)
	{
		for (unsigned int i = 0; i < framesPerBuffer; i++)
		{
			*out++ = data->audioBuffer->x[data->audioReadPos++];
		}
	}
	else if (data->wraparound == true)
	{
		unsigned int samplesLeft = Synthesizer::NUM_BUFFER_SAMPLES - data->audioReadPos;
		if (samplesLeft > framesPerBuffer)
		{
			for (unsigned int i = 0; i < framesPerBuffer; i++)
			{
				*out++ = data->audioBuffer->x[data->audioReadPos++];
			}
		}
		else
		{
			for (unsigned int i = 0; i < samplesLeft; i++)
			{
				*out++ = data->audioBuffer->x[data->audioReadPos++];
			}
			data->audioReadPos = 0;
			if (data->audioWritePos - data->audioReadPos > framesPerBuffer - samplesLeft)
			{
				for (unsigned int i = 0; i < framesPerBuffer - samplesLeft; i++)
				{
					*out++ = data->audioBuffer->x[data->audioReadPos++];
				}				
			}
			data->wraparound = false;
		}
	}
	else
	{
		for (unsigned int i = 0; i < framesPerBuffer; i++)
		{
			*out++ = 0;
		}
	}


	return 0;
}

// ****************************************************************************
/// Portaudio callback function for playback of selected audio signal only.
// ****************************************************************************
int  Synthesizer::playAudioSegment(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo*,
	PaStreamCallbackFlags,
	void *userData)
{
	(void)inputBuffer; /* Prevent unused variable warning. */
	signed short* out = (signed short*)outputBuffer;
	paData *data = (paData*)userData;
	unsigned int i;
	int finished;

	unsigned int framesLeft = data->audioWritePos - data->audioReadPos;

	if (framesLeft < framesPerBuffer)
    {
		/* final buffer... */
		for (i = 0; i<framesLeft; i++)
		{
			*out++ = data->audioBuffer->x[data->audioReadPos++];
		}
		for (; i<framesPerBuffer; i++)
		{
			*out++ = 0;
		}
		 finished = paComplete;
	}
	else
	{
		for (i = 0; i<framesPerBuffer; i++)
		{
			*out++ = data->audioBuffer->x[data->audioReadPos++];
		}
		finished = paContinue;
	}
	return finished;
}

// ****************************************************************************
/// Play the selected part of the buffer.
// ****************************************************************************

void Synthesizer::playSelection(int startPos, int stopPos)
{
	// Close any open stream, just to be safe.
	Pa_CloseStream(&audioStream);
	audioData.audioReadPos = startPos;
	audioData.audioWritePos = stopPos;
	PaStreamParameters outputParams = { 0 };
	outputParams.device = Pa_GetDefaultOutputDevice();
	outputParams.channelCount = 1;
	outputParams.sampleFormat = paInt16;
	outputParams.suggestedLatency = Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice())->defaultLowOutputLatency;
	outputParams.hostApiSpecificStreamInfo = NULL;
	printf("\n=== Playing selected audio signal... ===\n"); fflush(stdout);
	Pa_OpenStream(&audioStream, NULL, &outputParams, SAMPLING_RATE,
		paFramesPerBufferUnspecified, paNoFlag, &playAudioSegment, &audioData);
	Pa_StartStream(audioStream);
	printf("Waiting for playback to finish.\n"); fflush(stdout);
	while (Pa_IsStreamActive(audioStream) == 1) Pa_Sleep(100);
	printf("Done.\n"); fflush(stdout);
}

// ****************************************************************************
/// Portaudio callback function to record audio buffer.
// ****************************************************************************
int  Synthesizer::recordAudio(const void *inputBuffer, void *outputBuffer,
	unsigned long framesPerBuffer,
	const PaStreamCallbackTimeInfo*,
	PaStreamCallbackFlags,
	void *userData)
{
	(void)outputBuffer; /* Prevent unused variable warning. */
	signed short* in = (signed short*)inputBuffer;
	paData *data = (paData*)userData;

	for (unsigned int i = 0; i < framesPerBuffer; i++)
	{
		data->audioBuffer->x[data->audioWritePos++]  = *in++;
	}

	return 0;
}

// ****************************************************************************

