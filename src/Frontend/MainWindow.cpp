#include "MainWindow.h"

#include <wx/statline.h>
#include <wx/textfile.h>
#include <wx/tokenzr.h>
#include <iostream>
#include <fstream>
#include "Data.h"

using namespace std;

// ****************************************************************************
// IDs.
// ****************************************************************************

// Menu
static const int IDM_LOAD_AUDIO             = 1200;
static const int IDM_SAVE_AUDIO             = 1201;
static const int IDM_EXPORT_SENSOR_DATA		= 1202;
static const int IDM_IMPORT_SENSOR_DATA		= 1203;
static const int IDM_LOAD_REF_AF			= 1204;
static const int IDM_SAVE_AF_COEFF			= 1205;
static const int IDM_LOAD_AF_COEFF			= 1206;

static const int IDM_COM_CONFIG				= 1207;
static const int IDM_PALATE_CONFIG			= 1208;

static const int IDM_SAVE_AF_PARAMS			= 1209;
static const int IDM_LOAD_AF_PARAMS			= 1210;

static const int IDM_LOAD_MAP				= 1211;

// Keys
static const int IDK_CTRL_LEFT              = 1300;
static const int IDK_CTRL_RIGHT             = 1301;

// Others
static const int ID_NOTEBOOK                = 1500;


// ****************************************************************************
// The event table.
// ****************************************************************************

BEGIN_EVENT_TABLE(MainWindow, wxFrame)
	EVT_CLOSE(MainWindow::OnCloseWindow)

	// Menu events
	EVT_MENU(IDM_SAVE_AF_PARAMS, MainWindow::OnSaveAfParams)
	EVT_MENU(IDM_LOAD_AF_PARAMS, MainWindow::OnLoadAfParams)
	EVT_MENU(wxID_EXIT,  MainWindow::OnExit)
  
	EVT_MENU(wxID_ABOUT, MainWindow::OnAbout)

END_EVENT_TABLE()

// ****************************************************************************
// ****************************************************************************

MainWindow::MainWindow(wxWindow *parent, wxWindowID id, const wxString &title)
{
    data = Data::getInstance();

    // ****************************************************************
    // ****************************************************************

    wxFrame::Create(parent, id, title, wxDefaultPosition,
    wxDefaultSize, wxCLOSE_BOX | wxMINIMIZE_BOX |
    wxMAXIMIZE_BOX | wxSYSTEM_MENU | wxCAPTION |
    wxRAISED_BORDER | wxRESIZE_BORDER);

    initWidgets();
    updateWidgets();

    // Make the main window double buffered to avoid any flickering
    // of the child-windows and during resizing.
    this->SetDoubleBuffered(true);
}

// ****************************************************************************
// ****************************************************************************

void MainWindow::initWidgets()
{
  this->SetSize(20, 20, 1024, 830);

  wxMenu *menu = NULL;
  
  menuBar = new wxMenuBar();

  menu = new wxMenu();
  menu->Append(IDM_SAVE_AF_PARAMS, "Save AF parameter list");
  menu->Append(IDM_LOAD_AF_PARAMS, "Load AF parameter list");
  menu->AppendSeparator();
  menu->Append(wxID_EXIT, "Exit");

  menuBar->Append(menu, "File");

  // ****************************************************************
  menu = new wxMenu();
  menu->Append(wxID_ABOUT, "About");
  menuBar->Append(menu, "Help");

  // ****************************************************************
  this->SetMenuBar(menuBar);

  // ****************************************************************
  // Create the page(s).
  // ****************************************************************

  notebook = new wxNotebook(this, ID_NOTEBOOK);

  notebook->SetDoubleBuffered(true);

  mainPage = new MainPage(notebook);
  notebook->AddPage((wxPanel*)mainPage, "Main page", true);
	
  // Load set of vowel parameters from file
  loadAfParamsList(data->programPath + "Default.params");
}

// ****************************************************************************
// ****************************************************************************

void MainWindow::updateWidgets()
{
  mainPage->updateWidgets();
}

// ****************************************************************************
/// \brief Saves the current list of parameters of the area function.
///
/// This function saves the current list of area function parameters to a file.
/// \param filename name of the text file to write to.
/// \return true if file was successfully saved, otherwise false
// ****************************************************************************
bool MainWindow::saveAfParamsList(const char* filename)
{
	// ****************************************************************
	wxTextFile* csvFile = new wxTextFile(filename);

	// Choose separator used in CSV
	wxString separator = "  ";

	if (csvFile->Exists())
	{
		csvFile->Clear();
	}
	else if (!(csvFile->Create()))
	{
		printf("Error in saveCsvFile(): File '%s' already exists.\n", filename);
		return false;
	}	

	// First two lines are general information
	csvFile->AddLine(wxString("Parameters of one dimensional area functions for SecondVoicePc."));
	csvFile->AddLine(wxString("Label Llar_cm  Alar_cm2  powLar  xp_cm  Ap_cm2  powP  xc_cm  Ac_cm2  powC  xa_cm  Aa_cm2  powA  xin_cm  Ain_cm2  Lvt_cm  Alip_cm2"));

	wxString newLine = "";
	for (unsigned int i = 0; i < mainPage->afShapeVec.size(); i++)
	{
		newLine = mainPage->afShapeVec.at(i).name + separator;
		for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
		{
			newLine += wxString::Format(wxT("%f"), mainPage->afShapeVec.at(i).param[j].x);
			if (j < OneDimAreaFunction::NUM_AF_PARAMS - 1)
			{
				newLine += separator;
			}
		}
		// Write parameters to file (one line per set)
		csvFile->AddLine(newLine);
	}
	
	// Make changes to file on disk
	csvFile->Write();
	// Close file
	csvFile->Close();

	wxMessageBox(wxString::Format(wxT("File written to %s"), filename), "Done");

	return true;
}

//
//// ****************************************************************************
///// \brief Loads a list of parameters of the area function.
/////
///// This function loads a list of area function parameters from a file.
///// \param filename name of the text file to read from.
///// \return true if file was successfully loaded, otherwise false
//// ****************************************************************************
bool MainWindow::loadAfParamsList(const char* filename)
{
	wxTextFile *paramFile = new wxTextFile(filename);
	OneDimAreaFunction::Shape afShape;

	if (!paramFile->Open())
	{
		printf("Error in loadAfParamsList(): File '%s' does not exist", filename);
		return false;
	}

	mainPage->afShapeVec.clear();
	int numShapes = paramFile->GetLineCount();
	wxString separator = "  ";
	wxStringTokenizer *tokenizer = new wxStringTokenizer();
	wxString currentLine = "";
	wxString token = "";
	/* First two lines hold no parameters */
	currentLine = paramFile->GetFirstLine();
	currentLine = paramFile->GetNextLine();


	/* Each following line holds one area function shape */
	for (int i = 2; i < numShapes; i++)
	{
		currentLine = paramFile->GetNextLine();
		tokenizer->SetString(currentLine, separator);

		token = tokenizer->GetNextToken();
		afShape.name = token;
		for (int j = 0; j < OneDimAreaFunction::NUM_AF_PARAMS; j++)
		{
			token = tokenizer->GetNextToken();
			token.ToDouble(&afShape.param[j].x);
		}
		mainPage->afShapeVec.push_back(afShape);
	}


	paramFile->Close();

	mainPage->refreshAfShapeList();

	return true;
}

// ****************************************************************************
// ****************************************************************************

void MainWindow::OnCloseWindow(wxCloseEvent&)
{
  if (wxMessageBox("Do you really want to quit?", "Quit", wxYES_NO, this) == wxYES)
  {
    this->Destroy();
    exit(0);
  }
}


// ****************************************************************************
/// \brief Saves the current list of area function parameters as a csv file.
// ****************************************************************************
void MainWindow::OnSaveAfParams(wxCommandEvent&)
{
	wxFileDialog dialog(this, "Save an AF parameter list",
		paramsFileName.GetPath(), paramsFileName.GetFullName(),
		"params files (*.params)|*.params", wxFD_SAVE | wxFD_OVERWRITE_PROMPT);

	if (dialog.ShowModal() == wxID_OK)
	{
		paramsFileName = wxFileName(dialog.GetPath());

		if (!saveAfParamsList(paramsFileName.GetFullPath().c_str()))
		{
			wxMessageBox("Error saving the file.", "Error!");
		}
	}
}

//// ****************************************************************************
///// \brief Loads a list of area function parameters as a csv file.
//// ****************************************************************************
void MainWindow::OnLoadAfParams(wxCommandEvent&)
{
	wxString name = wxFileSelector("Open an AF parameter list", paramsFileName.GetPath(),
		paramsFileName.GetFullName(), "params", "Params-files (*.params)|*.params",
		wxFD_OPEN | wxFD_FILE_MUST_EXIST, this);

	if (name.empty() == false)
	{
		paramsFileName = wxFileName(name);

		if (!loadAfParamsList(paramsFileName.GetFullPath().c_str()))
		{
			wxMessageBox("Error in loading the file.", "Error!");
		}
		data->refreshAreaFunction();
		updateWidgets();
	}
}

// ****************************************************************************
// ****************************************************************************

void MainWindow::OnExit(wxCommandEvent&)
{
	Close(true);
}

// ****************************************************************************
// ****************************************************************************

void MainWindow::OnAbout(wxCommandEvent&)
{
  wxMessageDialog dialog(this, 
    "This is SecondVoicePc, developed by the IAS, TU Dresden, Germany.");
  dialog.ShowModal();
}

// ****************************************************************************
