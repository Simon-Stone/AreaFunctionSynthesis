#ifndef __SYNTHESIS_THREAD_H__
#define __SYNTHESIS_THREAD_H__

#include <wx/thread.h>
#include <wx/wx.h>
#include "../Backend/Synthesizer.h"

// Define the IDs for events that are sent to the GUI treat.
static const int SYNTHESIS_THREAD_EVENT = 54321;

static wxMutex *synthesizerProtection = new wxMutex();

class SynthesisThread : public wxThread
{
  // **************************************************************************
  // **************************************************************************
public:
	SynthesisThread(wxWindow *window, Synthesizer *synthesizer);
	// Thread execution starts here
	virtual void *Entry();

	// Called when the thread exits (normally or with Delete()), but not when it is killed
	virtual void OnExit();

	// **************************************************************************
	// **************************************************************************

private:
	wxWindow *window;
	Synthesizer *synthesizer;
	bool canceled;
};

// ****************************************************************************

#endif
