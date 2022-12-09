#include "MainPage.h"

#include <stdio.h>
#include <math.h>
#include <map>
#include <wx/statline.h>
#include <wx/numdlg.h>
#include <wx/dynarray.h>
#include <thread>

using namespace std;

// ****************************************************************************
// Ids.
// ****************************************************************************


// Buttons
static const int IDB_ADD_PARAM = 1000;
static const int IDB_START_STOP = 1012;

// Checkboxes
static const int IDC_ADD_F0_WOBBLE = 1100;
static const int IDC_NASAL_MODEL = 1103;
static const int IDC_SHOW_AREA_FUNC = 1105;
static const int IDC_SHOW_TUBE_FUNC = 1108;

// Slider and spin controls
static const int IDS_BASE_F0 = 1200;
static const int IDS_GLOTTAL_OPENING = 1201;
static const int IDS_NASAL_AREA_MODEL = 1202;

// Other
static const int IDL_AF_SHAPE_LIST = 1500;


// ****************************************************************************
// The event table.
// ****************************************************************************

BEGIN_EVENT_TABLE(MainPage, wxPanel)
    // Custom event handler for update requests by child widgets.
    // Left side controls
    EVT_BUTTON(IDB_ADD_PARAM, MainPage::OnAddParam)
    EVT_BUTTON(IDB_START_STOP, MainPage::OnStartStop)
    EVT_CHECKBOX(IDC_ADD_F0_WOBBLE, MainPage::OnAddF0Wobble)
    EVT_CHECKBOX(IDC_NASAL_MODEL, MainPage::OnNasality)
    EVT_CHECKBOX(IDC_SHOW_AREA_FUNC, MainPage::OnShowAreaFunc)
    EVT_CHECKBOX(IDC_SHOW_TUBE_FUNC, MainPage::OnShowTubeFunc)
    EVT_COMBOBOX(IDL_AF_SHAPE_LIST, MainPage::OnShapeSelect)
    EVT_SLIDER(IDS_BASE_F0, MainPage::OnBaseF0Change)
    EVT_SLIDER(IDS_GLOTTAL_OPENING, MainPage::OnGlottalOpeningChange)
    EVT_SPINCTRL(IDS_NASAL_AREA_MODEL, MainPage::OnNasalArea)
END_EVENT_TABLE()



// ****************************************************************************
/// Constructor.
// ****************************************************************************


MainPage::MainPage(wxWindow *parent) :
        wxPanel(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxCLIP_CHILDREN) {
    data = Data::getInstance();
    initWidgets();
    updateWidgets();
}


// ****************************************************************************
// ****************************************************************************

void MainPage::updateWidgets() {
    this->Refresh();
}


// ****************************************************************************
/// Init all widgets on this page.
// ****************************************************************************

void MainPage::initWidgets() {
    isNasalModel = false;

    wxSizer *sizer;
    wxStaticText *label;
    wxButton *button;
    wxCheckBox *checkbox;

    // ****************************************************************
    // Left side sizer.
    // ****************************************************************
    wxBoxSizer *leftSizer = new wxBoxSizer(wxVERTICAL);

    leftSizer->AddSpacer(10);

    checkbox = new wxCheckBox(this, IDC_SHOW_AREA_FUNC, "Show area function");
    checkbox->SetValue(1);
    leftSizer->Add(checkbox, 0, wxALL, 2);

    checkbox = new wxCheckBox(this, IDC_SHOW_TUBE_FUNC, "Show tube function");
    leftSizer->Add(checkbox, 0, wxALL, 2);


    // ****************************************************************
    leftSizer->AddSpacer(10);
    label = new wxStaticText(this, wxID_ANY, "Vocal cord displacement [µm]:");
    leftSizer->Add(label, 0, wxALL, 2);
    glottalDisplacement = new wxSlider(this, IDS_GLOTTAL_OPENING, 100, -500, 3000, wxDefaultPosition, wxDefaultSize,
                                       wxSL_MIN_MAX_LABELS | wxSL_VALUE_LABEL);
    leftSizer->Add(glottalDisplacement, 0, wxALL, 2);

    leftSizer->AddSpacer(10);
    label = new wxStaticText(this, wxID_ANY, "Base F0 [Hz]:");
    leftSizer->Add(label, 0, wxALL, 2);
    baseF0Slider = new wxSlider(this, IDS_BASE_F0, 120, 40, 600, wxDefaultPosition, glottalDisplacement->GetSize(),
                                wxSL_MIN_MAX_LABELS | wxSL_VALUE_LABEL);
    leftSizer->Add(baseF0Slider, 0, wxALL, 2);

    leftSizer->AddSpacer(10);
    checkbox = new wxCheckBox(this, IDC_ADD_F0_WOBBLE, "Add F0 wobble");
    checkbox->SetValue(1);
    leftSizer->Add(checkbox, 0, wxALL, 2);

    leftSizer->AddSpacer(10);
    wxFlexGridSizer *nasalSizer = new wxFlexGridSizer(2);
    checkbox = new wxCheckBox(this, IDC_NASAL_MODEL, "Port area in mm²:");
    nasalSizer->Add(checkbox, 1, wxGROW | wxALL, 2);
    nasalAreaModel_mm2 = new wxSpinCtrl(this, IDS_NASAL_AREA_MODEL, "100", wxDefaultPosition, wxSize(50, -1),
                                        wxSP_ARROW_KEYS, 0, 200, 100);
    nasalSizer->Add(nasalAreaModel_mm2, 0, wxALL, 2);
    leftSizer->Add(nasalSizer, 0, wxALIGN_LEFT);

    leftSizer->AddSpacer(10);
    wxStaticLine *horizontalLine = new wxStaticLine(this, wxID_ANY,
                                                    wxDefaultPosition, wxDefaultSize, wxLI_HORIZONTAL);
    leftSizer->Add(horizontalLine, 0, wxGROW);

    leftSizer->AddSpacer(10);

    // ****************************************************************
    // List with pre-defined parameter sets for the 1D area function.
    // ****************************************************************
    wxStaticText *textLabel = new wxStaticText(this, wxID_ANY, wxString("Pre-defined VT shapes"));
    leftSizer->Add(textLabel, 0, wxALL | wxGROW, 5);


    afShapeList = new wxComboBox(this, IDL_AF_SHAPE_LIST, wxEmptyString,
                                 wxDefaultPosition, wxDefaultSize, wxArrayString(), wxCB_SORT);
    leftSizer->Add(afShapeList, 0, wxALL | wxGROW, 5);

    button = new wxButton(this, IDB_ADD_PARAM, "Add current parameters");
    leftSizer->Add(button, 0, wxALL | wxGROW, 5);

    leftSizer->AddStretchSpacer();

    button = new wxButton(this, IDB_START_STOP, "\nStart/Stop\n");
    leftSizer->Add(button, 0, wxALL | wxGROW, 5);


    // ****************************************************************
    // Static line to separate the left and right part.
    // ****************************************************************

    wxStaticLine *verticalLine = new wxStaticLine(this, wxID_ANY, wxDefaultPosition,
                                                  wxDefaultSize, wxLI_VERTICAL);

    // ****************************************************************
    // Panel on the right side
    // ****************************************************************
    wxPanel *rightPanel = new wxPanel(this, wxID_ANY, wxDefaultPosition, wxDefaultSize);

    wxBoxSizer *rightSizer = new wxBoxSizer(wxHORIZONTAL);
    rightPanel->SetSizer(rightSizer);

    areaFunctionPicture = new AreaFunctionPicture(rightPanel);
    areaFunctionPicture->SetMinSize(wxSize(-1, 25));
    rightSizer->Add(areaFunctionPicture, 1, wxGROW | wxALL, 2);

    // ****************************************************************
    // Top level sizer.
    // ****************************************************************

    wxBoxSizer *topLevelSizer = new wxBoxSizer(wxHORIZONTAL);

    topLevelSizer->Add(leftSizer, 0, wxGROW | wxALL, 5);
    topLevelSizer->Add(verticalLine, 0, wxGROW | wxALL, 2);
    topLevelSizer->Add(rightPanel, 1, wxGROW | wxALL, 2);

    this->SetSizerAndFit(topLevelSizer);
}

// ****************************************************************************
/// Refresh the list of currently loaded parameters.
// ****************************************************************************

void MainPage::refreshAfShapeList() {
    afShapeList->Clear();

    for (unsigned int i = 0; i < afShapeVec.size(); i++) {
        afShapeList->AppendString(afShapeVec.at(i).name);
    }
}

// ****************************************************************************
// ****************************************************************************

void MainPage::OnBaseF0Change(wxCommandEvent &) {
    double baseF0_Hz = baseF0Slider->GetValue();
	data->synthesizer->baseF0_Hz = baseF0_Hz;
}

// ****************************************************************************
// ****************************************************************************

void MainPage::OnStartStop(wxCommandEvent &) {
    if (data->synthesizer->toggleRun()) {
        data->analysisMode = Data::REAL_TIME;
        synthesisThread = new SynthesisThread(this, data->synthesizer);
        synthesisThread->Run();
    } else {

    }
}

// ****************************************************************************
// ****************************************************************************
void MainPage::OnGlottalOpeningChange(wxCommandEvent &) {
	double displacement_cm = (double) glottalDisplacement->GetValue() / (double) 10000;
	data->synthesizer->glottis->controlParam[TriangularGlottis::REST_DISP_1].x = displacement_cm;
	data->synthesizer->glottis->controlParam[TriangularGlottis::REST_DISP_2].x = displacement_cm;
}

// ****************************************************************************
/// Toggle the visibility of the continuous function.
// ****************************************************************************
void MainPage::OnShowAreaFunc(wxCommandEvent &event) {
    areaFunctionPicture->showAreaFunction = event.IsChecked();
}

// ****************************************************************************
/// Toggle the visibility of the tube function.
// ****************************************************************************
void MainPage::OnShowTubeFunc(wxCommandEvent &event) {
    areaFunctionPicture->showTubeFunction = event.IsChecked();
    updateWidgets();
}

void MainPage::OnAddF0Wobble(wxCommandEvent &event) {
    data->synthesizer->addF0Wobble = event.IsChecked();
}

// ****************************************************************************
/// Load parameter selected from dropdown list.
// ****************************************************************************
void MainPage::OnShapeSelect(wxCommandEvent &) {
    wxString name = afShapeList->GetValue();
    for (unsigned int i = 0; i < afShapeVec.size(); i++) {
        if (afShapeVec.at(i).name == name) {
            data->synthesizer->vtAreaFunction->setParameters(afShapeVec.at(i).param);
        }
    }
    data->refreshAreaFunction();
    updateWidgets();
}

// ****************************************************************************
/// Add the currently selected parameters of the 1D area function to the list.
// ****************************************************************************
void MainPage::OnAddParam(wxCommandEvent &) {
    OneDimAreaFunction::Shape newShape;

    data->synthesizer->vtAreaFunction->getParameters(&newShape.param[0]);

    // Ask user for a label
    wxTextEntryDialog *nameDlg;
    nameDlg = new wxTextEntryDialog(this, "Enter a name", "Name for current set");
    if (nameDlg->ShowModal() == wxID_OK) {
        newShape.name = nameDlg->GetValue();
        // Add parameter set to the list
        afShapeVec.push_back(newShape);
        refreshAfShapeList();
    }

}

// ****************************************************************************
/// Open or close the nasal port.
// ****************************************************************************
void MainPage::OnNasality(wxCommandEvent &event) {
    if (event.GetId() == IDC_NASAL_MODEL && event.IsChecked()) {
        isNasalModel = true;
        data->tlModel->tube.section[Tube::FIRST_NOSE_SECTION]->area_cm2 = nasalAreaModel_mm2->GetValue() / 100.0;
        data->synthesizer->tubeFunction->setVelumOpening(nasalAreaModel_mm2->GetValue() / 100.0);
    } else if (event.GetId() == IDC_NASAL_MODEL) {
        isNasalModel = false;
        data->tlModel->tube.section[Tube::FIRST_NOSE_SECTION]->area_cm2 = 0.0;
        data->synthesizer->tubeFunction->setVelumOpening(0.0);
    }
}

// ****************************************************************************
/// Change the opening area of the velopharyngeal port.
// ****************************************************************************
void MainPage::OnNasalArea(wxSpinEvent &event) {
    if (event.GetId() == IDS_NASAL_AREA_MODEL && isNasalModel) {
        data->tlModel->tube.section[Tube::FIRST_NOSE_SECTION]->area_cm2 = nasalAreaModel_mm2->GetValue() / 100.0;
        data->synthesizer->tubeFunction->setVelumOpening(nasalAreaModel_mm2->GetValue() / 100.0);
    }
}

// ****************************************************************************
