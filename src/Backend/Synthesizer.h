#ifndef __SYNTHESIZER_H__
#define __SYNTHESIZER_H__

#include "wx/timer.h"
#include "wx/wx.h"
#include "Constants.h"
#include "OneDimAreaFunction.h"
#include "TdsModel.h"
#include "TlModel.h"
#include "LfPulse.h"
#include "Tube.h"
#include "TriangularGlottis.h"
#include "Dsp.h"
#include "IirFilter.h"
//#include "SerialPort.h"
#include <vector>
#include <queue>
#include "wx/msw/winundef.h" // Prevents accidental redefinition by windows files
#include "portaudio.h"


using namespace std;


// ****************************************************************************
/// With this class, the user can synthesize a speech signal by incrementally
/// calculating short pieces of the signal in which the current vocal tract and 
/// glottis settings are interpolated towards new settings.
// ****************************************************************************

class Synthesizer : public wxTimer
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
	static const int MAX_NUM_CONTACT_SENSORS = 64;
	enum EpgLayout
	{
		EPG52,
		EPG32,
		NUM_EPG_LAYOUTS
	};
	inline constexpr static int CONTACT_NUMBER[NUM_EPG_LAYOUTS] = { 52, 32 };
	EpgLayout currentEpgLayout;

	static const int MAX_NUM_OPG_SENSORS = 5;
	enum OpgLayout
	{
		DIST5,
		NUM_OPG_LAYOUTS
	};
	inline constexpr static int DISTANCE_NUMBER[NUM_OPG_LAYOUTS] = { 5 };
	OpgLayout currentOpgLayout;

	static const int MAX_OUTLINE_POINTS = 64;

	static const int FRAME_RATE_HZ = 100;
	static const int NUM_BUFFER_SECONDS = 600;
	static const int NUM_BUFFER_SAMPLES = NUM_BUFFER_SECONDS * SAMPLING_RATE;
	static const int NUM_BUFFER_FRAMES = NUM_BUFFER_SECONDS * FRAME_RATE_HZ;
	static const int AUDIO_UPDATE_RATE_HZ = 20;
	static const int AUDIO_FRAME_LENGTH = SAMPLING_RATE / AUDIO_UPDATE_RATE_HZ;


	static const int NUM_TUBE_SECTIONS = Tube::NUM_PHARYNX_MOUTH_SECTIONS;
	/// Coupling section for the nasal cavity.
	static const int NUM_PHARYNX_SECTIONS = Tube::NUM_PHARYNX_SECTIONS;

  // Palate configuration.
  struct Palate
  {
    string fileName;
    int numOutlinePoints;
    double outlineX_cm[MAX_OUTLINE_POINTS];
    double outlineY_cm[MAX_OUTLINE_POINTS];
    double sensorX_cm[MAX_NUM_OPG_SENSORS];
    double sensorY_cm[MAX_NUM_OPG_SENSORS];
    double sensorAngle_deg[MAX_NUM_OPG_SENSORS];
  };

  // A frame of sensor data.
  struct SensorFrame
  {
    int frameIndex;
    
    double lipOpening_rel;    // The "unit" rel means in the interval [0, 1]
    double lipProtrusion_rel;
    double velicOpening_rel;
    double distance_cm[MAX_NUM_OPG_SENSORS];
    unsigned char contact[MAX_NUM_CONTACT_SENSORS];
    
    // Source parameters.
    double f0_Hz;
    bool voicing;
    double breathiness_rel;
    double pressure_Pa;
  };
  
	// ****************************************************************

	struct TubeSection
	{
		double area;
		double circ;
		double pos;
		double length;
		Tube::Articulator articulator;
		double laterality;
	};

	typedef struct
	{
		Signal16 *audioBuffer;
		unsigned int audioWritePos;
		unsigned int audioReadPos;
		bool wraparound;
	}
	paData;

	enum AudioMode
	{
		AUDIO_RECORDING,
		AUDIO_SYNTHESIS,
		AUDIO_PLAYBACK,
		NUM_AUDIO_MODE_OPTIONS
	};

	enum MappingOption
	{
		REGRESSION,
		K_NEAREST,
		NUM_MAPPING_OPTIONS
	};

	// Variables.	
	SensorFrame *frameBuffer;
	int numRecordedFrames;

  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
	Synthesizer(Signal16 *signalBuffer);
	~Synthesizer();
	void reset();
	void synthesizeSegment();
	int getAudioMode();
	int setAudioMode(AudioMode newAudioMode);
	int setMappingOption(MappingOption mappingOption);
	void startRealTimeAudio();
	void startAudioStream();
	void stopRealTimeAudio();
	void stopAudioStream();
	
	void playSelection(int startPos, int stopPos);
	void playTargetSequence(OneDimAreaFunction::Shape *targetShape,
		double *duration_s, double *transition_s);

	void synthesizeSignalTds(Tube *newTube, double *newGlottisParams, 
		int numNewSamples, double *newSignal);
	void synthesizeVowelLf(TlModel *tlModel, LfPulse &lfPulse, 
		int startPos, bool isLongVowel);
	
	bool toggleRun();
	void clearSensorFrame(SensorFrame &f);


	bool configureSerialPort(wxString portName = "COM1", int baudRate = 115200);
	bool openSerialPort();
	bool closeSerialPort();

	void sensorDataToTube(SensorFrame *sensorData, Tube *tube);
	void sensorDataToGlottisParams(SensorFrame *sensorData, double *glottisParams);


  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
//	SerialPort serialPort;
//	SerialPort::PortSettings comSettings;
	// The analysis mode
	enum DataMode
	{
		// Data frames from hardware contain already processed values (distances and relative openings and so forth)
		PROCESSED_DATA,
		// Data frames from hardware contain raw ADC values
		RAW_DATA,
		NUM_DATA_MODES
	};
	std::string DATA_MODE_NAME[NUM_DATA_MODES] = { "Processed data", "Raw data" };
	DataMode dataMode;
	// Use sensor data to control the area function or allow manual manipulation of the control points
	enum ModelControlMode
	{
		SENSOR_DATA,
		MANUAL,
		NUM_CONTROL_MODES
	};
	ModelControlMode modelControlMode;
	// The currently processed frame
	SensorFrame sensorFrame;
	// Model area function of the vocal tract
	OneDimAreaFunction *vtAreaFunction;
	// Tube function of vocal tract
	Tube *tubeFunction;
	// Time domain simulation model of the vocal tract
	TdsModel *tdsModel;
	// Triangular glottis model
	TriangularGlottis *glottis;
	// Base F0 for synthesis
	double baseF0_Hz;
	// If true, a time-varying wobble will be superimposed on the F0
	bool addF0Wobble;
	// Flag to switch between manual and sensor driven area function manipulation
	bool isUsingSensorData;

	// Audio data object for Portaudio output stream
	paData audioData;
	
	// List of training samples for K Nearest Neighbors classification
	vector<vector<double>> kNearestSamples;
	// List of area function parameters associated with each training sample
	vector<vector<OneDimAreaFunction::Param>> kNearestModels;
	
	bool isRunning;

	// **************************************************************************
	// Private functions.
	// **************************************************************************

private:

	void assignSensorRxData();
	void calculateEpgFeatures(SensorFrame *sensorData, double &frontness, double &backness, double &sum);
	void calculateAfParametersFromSensorData(SensorFrame *sensorData, 
		OneDimAreaFunction::Param *afParams);
	int findNearestTrainingSample(SensorFrame *sensorData);
	void getParametersFromTrainingSample(int sampleIndex, OneDimAreaFunction::Param *afParams);
	void init();	
	void initAudio();
	void interpolateParameters(OneDimAreaFunction::Param *p0,
		OneDimAreaFunction::Param *p1, OneDimAreaFunction::Param *pX,
		double t0, double t1, double tX);
	void fetchSensorData();
	static int  playAudioLoop(const void *inputBuffer, void *outputBuffer,
		unsigned long framesPerBuffer,
		const PaStreamCallbackTimeInfo* timeInfo,
		PaStreamCallbackFlags statusFlags,
		void *userData);
	static int playAudioSegment(const void * inputBuffer, void * outputBuffer, 
		unsigned long framesPerBuffer, 
		const PaStreamCallbackTimeInfo *timeInfo,
		PaStreamCallbackFlags, 
		void * userData);
	static int  recordAudio(const void *inputBuffer, void *outputBuffer,
		unsigned long framesPerBuffer,
		const PaStreamCallbackTimeInfo* timeInfo,
		PaStreamCallbackFlags statusFlags,
		void *userData);

	// **************************************************************************
	// Private data.
	// **************************************************************************

private:
	int playbackMode;
	Tube prevTube;
	Tube tube;
	TubeSection tubeSection[NUM_TUBE_SECTIONS];
	double prevGlottisParams[1024];
	
	static const int TDS_BUFFER_LENGTH = 256;
	static const int TDS_BUFFER_MASK = 255;
	static const int MAX_BUFFER_BYTES = 4096;
	// Important: Buffers must be UNSIGNED char !
	unsigned char rxBuffer[MAX_BUFFER_BYTES];
	int numRxBytes;

	double *glottisParams;
	double *newSignal;
	double outputFlow[TDS_BUFFER_LENGTH];
	double outputPressure[TDS_BUFFER_LENGTH];
	double filteredOutputPressure[TDS_BUFFER_LENGTH];

	OneDimAreaFunction::Shape *targetShape;
	double *stationary_s;
	double *transition_s;

	PaStream *audioStream;
	

	IirFilter outputPressureFilter;
	bool initialShapesSet;

	// Indicates which mapping technique to use from sensor data to vocal tract parameters
	MappingOption mappingOption;
	// Choose between audio recording (no synthesis) for training data recording or synthesis
	AudioMode audioMode;
	// Last found best sample
	int lastBestSampleIdx;


};

#endif

// ****************************************************************************
