#ifndef ___DATA__
#define ___DATA__

#include <wx/wx.h>
#include <string>
#include <vector>
#include "../Backend/Synthesizer.h"

using namespace std;


// ****************************************************************************
// Declare some custom event types.
// The updateRequestEvent is meant to be used by child widgets to request their
// parent to update some region or other child widgets.
// A command event of this type is posted as follows:
//
//   wxCommandEvent event(updateRequestEvent);
//   event.SetInt(REFRESH_PICTURES | REFRESH_PICTURES_AND_CONTROLS);
//   wxPostEvent(receiverWindow, event);
//
// The receiver window must have the following in his event table:
//
//  EVT_COMMAND(wxID_ANY, updateRequestEvent, OnUpdateRequest)
// 
// The function OnUpdateRequest(...) takes one wxCommandEvent parameter.
// ****************************************************************************

extern const wxEventType updateRequestEvent;
extern const int REFRESH_PICTURES;
extern const int UPDATE_PICTURES;
extern const int REFRESH_PICTURES_AND_CONTROLS;


// ****************************************************************************
/// Singleton class containing the data and common methods for the frontend 
/// classes.
// ****************************************************************************

class Data
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  double selectionMark_s[2];
  double mark_s;

  enum TrackType
  {
	  MAIN_TRACK,
	  EGG_TRACK,
	  EXTRA_TRACK,
	  NUM_TRACKS
};

  	// ****************************************************************
	// Oscillogram variables
	// ****************************************************************
	static const int TRACK_DURATION_S = Synthesizer::NUM_BUFFER_SECONDS;                    // Length of the tracks in s
	Signal16 *track[NUM_TRACKS];
	bool showTrack[NUM_TRACKS];
	int selectionMark_pt[2];
	int mark_pt;

	int centerPos_pt;     ///< Sampling index in the middle of the oscillograms
	int oscillogramVisTimeRange_pt;
	double oscillogramAmpZoom;

	// The synthesizer object handles all speech synthesis
	Synthesizer *synthesizer;
	// Transmission line model of the model area function
	TlModel *tlModel;
	// Transmission line model of the reference area function
	TlModel *tlRef;
	
	enum AnalysisMode
	{
		// Data is gathered from sensor device and displayed in real-time
		REAL_TIME,
		// Analysis mode when recorded data is being statistically evaluated (e.g, average frame calculation)
		STATISTICAL,
		NUM_ANALYSIS_MODES
	};
	// The current data analysis mode
	AnalysisMode analysisMode;
	// Holds the averaged frame if calculated
	Synthesizer::SensorFrame averageFrame;

	// Reference tube function
	Tube *refTubeFunction;

	// Target shape parameters for synthesis of manually set sound sequences
	OneDimAreaFunction::Shape targetShape[4];
	double stationary_s[4];
	double transition_s[3];

	ComplexSignal *primarySpectrum;
	ComplexSignal *refSpectrum;
	
//	FormantOptimizationDialog *formantOptimizationDialog;
//	SetTargetSequenceDialog *setTargetSequenceDialog;

  /// Path to the executable file
  wxString programPath;
  
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  static Data *getInstance();
  void init(wxString arg0);

  int getFrameIndexAtMark();
  bool isValidSelection();
  void calcAverageFrame();

  double getFormantError(double currentF1, double currentF2,
	  double currentF3, double targetF1, double targetF2, double targetF3);
  bool getVowelFormants(OneDimAreaFunction *vtAf, double &F1_Hz, double &F2_Hz,
		double &F3_Hz, double &minArea_cm2);
  void optimizeFormantsVowel(wxWindow *updateParent, OneDimAreaFunction *vtAf,
	  double targetF1, double targetF2, double targetF3,
		double maxParamChange_cm, double minAdvisedArea_cm2, bool paramFixed[]);

  void refreshAreaFunction();
  bool saveSpectrum(const char *filename, ComplexSignal *spectrum);

  // **************************************************************************
  // Private data.       
  // **************************************************************************

private:
  static Data *instance;

  // **************************************************************************
  // Private functions.
  // **************************************************************************

private:
  Data();
};

#endif

// ****************************************************************************

