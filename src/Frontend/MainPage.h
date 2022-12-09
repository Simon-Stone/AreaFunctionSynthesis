#ifndef __MAIN_PAGE_H__
#define __MAIN_PAGE_H__

// ****************************************************************************
// ****************************************************************************

#include <vector>

#include <wx/wx.h>
#include <wx/frame.h>
#include <wx/splitter.h>
#include <wx/spinctrl.h>
#include <wx/filename.h>
#include "Data.h"
#include "../Backend/OneDimAreaFunction.h"
#include "AreaFunctionPicture.h"
#include "SynthesisThread.h"


// ****************************************************************************
// ****************************************************************************

class MainPage : public wxPanel
{
public:
	std::vector<OneDimAreaFunction::Shape> afShapeVec;
	wxComboBox *afShapeList;

public:
	MainPage(wxWindow *parent);
	void refreshAfShapeList();
	void updateWidgets();

private:
	Data *data;
	wxSlider *glottalDisplacement;
	wxSlider *baseF0Slider;

	AreaFunctionPicture *areaFunctionPicture;
	wxSpinCtrl *nasalAreaModel_mm2;
	bool isNasalModel;
	SynthesisThread *synthesisThread;

private:
  void initWidgets();

  // Events.
  void OnAddF0Wobble(wxCommandEvent &event);
  void OnAddParam(wxCommandEvent &event);
  void OnBaseF0Change(wxCommandEvent &event);
  void OnGlottalOpeningChange(wxCommandEvent &event);
  void OnNasalArea(wxSpinEvent &event);
  void OnNasality(wxCommandEvent &event);
  void OnShapeSelect(wxCommandEvent &event);
  void OnShowAreaFunc(wxCommandEvent &event);
  void OnShowTubeFunc(wxCommandEvent &event);
  void OnStartStop(wxCommandEvent &event);

  void OnUpdateRequest(wxCommandEvent &event);

  DECLARE_EVENT_TABLE()
};


#endif