#include <wx/wx.h>
#include "MainWindow.h"

// ****************************************************************************
// ****************************************************************************

class Application : public wxApp
{
public:
  virtual bool OnInit();
  #if defined(__linux__) || defined(__APPLE__)
  virtual int OnExit();
  #endif
  void createConsole();
};

DECLARE_APP(Application);
IMPLEMENT_APP(Application);

// ****************************************************************************
