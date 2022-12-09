#include <fstream>
#include "SynthesisThread.h"


// ****************************************************************************
// ****************************************************************************

SynthesisThread::SynthesisThread(wxWindow *window, Synthesizer *synthesizer) : wxThread()
{
  this->window = window;
  this->synthesizer = synthesizer;
}


// ****************************************************************************
// ****************************************************************************

void *SynthesisThread::Entry()
{
	while (synthesizer->isRunning)
	{
		synthesizerProtection->Lock();
		synthesizer->synthesizeSegment();
		synthesizerProtection->Unlock();
		// ****************************************************************
		// Send an event to the GUI thread that a new chunk of data is done
		// ****************************************************************
		wxCommandEvent event(wxEVT_COMMAND_MENU_SELECTED, SYNTHESIS_THREAD_EVENT);
		event.SetInt(-1); // that's all
		wxPostEvent(window, event);
		wxMilliSleep(1); // Otherwise the GUI thread won't be fast enough
	}		

	return NULL;
}


// ****************************************************************************
// ****************************************************************************

void SynthesisThread::OnExit()
{

}

// ****************************************************************************
