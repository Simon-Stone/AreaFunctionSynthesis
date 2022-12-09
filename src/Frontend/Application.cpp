#include "Application.h"

#ifdef WIN32
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#endif

#include <iostream>
#include <fstream>
#include <stdio.h>

#include "Data.h"

using namespace std;

// ****************************************************************************
// ****************************************************************************

bool Application::OnInit()
{
  createConsole();
  printf("=== Console output for AreaFunctionSynthesis ===\n\n");

  // Init the data class at the very beginning.  
  Data *data = Data::getInstance();

  if (argc < 1)
  {
    printf("Error: At least the default command line argument is expected.\n");
  }
  data->init(argv[0]);

  // Create and show the main window.

  MainWindow *mainWindow = new MainWindow(NULL, wxID_ANY, "AreaFunctionSynthesis");
  SetTopWindow(mainWindow);
  mainWindow->Show();
  
  return true;
}

// ****************************************************************************
// ****************************************************************************

void Application::createConsole()
{

#ifdef WIN32

    // maximum mumber of lines the output console should have

  static const WORD MAX_CONSOLE_LINES = 500;
  int hConHandle;
  HANDLE lStdHandle;

  CONSOLE_SCREEN_BUFFER_INFO coninfo;

  // ****************************************************************

  // allocate a console for this app
  AllocConsole();

  // set the screen buffer to be big enough to let us scroll text
  GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &coninfo);

  coninfo.dwSize.Y = MAX_CONSOLE_LINES;
  SetConsoleScreenBufferSize(GetStdHandle(STD_OUTPUT_HANDLE), coninfo.dwSize);

  // ****************************************************************
  // redirect unbuffered STDOUT to the console
  // e.g. printf("...", ...);
  // ****************************************************************

  lStdHandle = GetStdHandle(STD_OUTPUT_HANDLE);
  hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);

  freopen("CONOUT$", "w", stdout);
  setvbuf(stdout, NULL, _IONBF, 0);

  // ****************************************************************
  // redirect unbuffered STDIN to the console
  // ****************************************************************

  lStdHandle = GetStdHandle(STD_INPUT_HANDLE);
  hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);

  freopen("CONOUT$", "w", stdin);
  setvbuf(stdin, NULL, _IONBF, 0);

  // ****************************************************************
  // redirect unbuffered STDERR to the console
  // ****************************************************************

  lStdHandle = GetStdHandle(STD_ERROR_HANDLE);
  hConHandle = _open_osfhandle((intptr_t)lStdHandle, _O_TEXT);
  freopen("CONOUT$", "w", stderr);
  setvbuf(stderr, NULL, _IONBF, 0);

  // ****************************************************************
  // make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog
  // point to console as well
  // ****************************************************************

  //ios::sync_with_stdio();

#endif
}

// ****************************************************************************

#if defined(__linux__) || defined(__APPLE__)
int Application::OnExit()
  {
    return fclose(stderr);
  }
#endif
