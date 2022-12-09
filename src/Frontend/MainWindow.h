#ifndef _MAIN_WINDOW_H_
#define _MAIN_WINDOW_H_

#include <wx/wx.h>
#include <wx/filename.h>
#include <wx/notebook.h>

//#include "Data.h"
#include "MainPage.h"


// ****************************************************************************
// Main window of the application.
// ****************************************************************************

class MainWindow : public wxFrame
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  MainWindow(wxWindow *parent, wxWindowID id = wxID_ANY, const wxString &title = wxEmptyString);
  
  void initWidgets();
  void updateWidgets();
    
  // ****************************************************************************
  // Private data.
  // **************************************************************************
	static const int REF_LENGTH = 129; // Length of area function reference from VTL
private:
  wxMenuBar *menuBar;
  wxNotebook *notebook;
  MainPage *mainPage;

  Data *data;
  wxFileName paramsFileName;

private:
	bool saveAfParamsList(const char* filename);
	bool loadAfParamsList(const char* filename);

  // Window events
	void OnCloseWindow(wxCloseEvent &event);

  // Menu functions
	void OnSaveAfParams(wxCommandEvent &event);
	void OnLoadAfParams(wxCommandEvent &event);

	void OnExit(wxCommandEvent &event);

	void OnAbout(wxCommandEvent &event);

  DECLARE_EVENT_TABLE()
};

#endif
