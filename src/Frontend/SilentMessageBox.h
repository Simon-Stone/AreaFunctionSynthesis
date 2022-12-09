#ifndef __SILENT_MESSAGE_BOX__
#define __SILENT_MESSAGE_BOX__

#include <wx/wx.h>
#include <wx/dialog.h>

// ****************************************************************************
/// This dialog displays a message box dialog analogue to ::wxMessageBox(...),
/// but WITHOUT playing a sound when it pops up.
// ****************************************************************************

class SilentMessageBox : public wxDialog
{
  // **************************************************************************
  // Public functions.
  // **************************************************************************

public:
  SilentMessageBox(const wxString& message, const wxString& caption = "Message", 
    wxWindow *parent = NULL);
};

#endif

// ****************************************************************************
