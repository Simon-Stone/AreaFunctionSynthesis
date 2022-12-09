#include "SerialPort.h"

#include <stdio.h>
#include <string.h>

// Disable warnings, if deprecated functions are used (e.g. printf() instead of s_printf()).
#pragma warning(disable : 4996)

// ****************************************************************************
/// Constructor. Variables are initialized here.
// ****************************************************************************

SerialPort::SerialPort()
{
  strcpy(settings.portName, "COM1");
  settings.baudRate = 115200;
  settings.parity = NO_PARITY;
  settings.stopBits = ONE_STOP_BIT;
  settings.hardwareFlowControl = false;
  settings.timeout_ms = 0;
  settings.timeout_s = 0;

  // Windows specific variables
  handle = INVALID_HANDLE_VALUE;
}

// ****************************************************************************
/// Destructor.
// ****************************************************************************

SerialPort::~SerialPort()
{
  if (isOpen())
  {
    close();
  }
}

// ****************************************************************************
/// Function variant to open the serial port.
// ****************************************************************************

bool SerialPort::open(const char *portName, int baudRate, Parity parity, 
  StopBits stopBits, bool hardwareFlowControl)
{
  PortSettings s;

  strcpy(s.portName, portName);
  s.baudRate = baudRate;
  s.parity = parity;
  s.stopBits = stopBits;
  s.hardwareFlowControl = hardwareFlowControl;
  s.timeout_ms = 0;
  s.timeout_s = 0;

  return open(s);
}

// ****************************************************************************
/// Opens the serial port with the given settings.
// ****************************************************************************

bool SerialPort::open(PortSettings s)
{
  if (isOpen())
  {
    close();
  }

  handle = INVALID_HANDLE_VALUE;
  settings = s;

  // ****************************************************************

  unsigned long confSize = sizeof(COMMCONFIG);
  commConfig.dwSize = confSize;

  handle = CreateFileA(
    s.portName, 
    GENERIC_READ | GENERIC_WRITE, 
    0,
    0,
    OPEN_EXISTING, 
    0, 
    NULL
  );

  if (handle == INVALID_HANDLE_VALUE)
  {
    return false;
  }
  
  // Set the reccommended queue size for buffering

  const DWORD IN_QUEUE_SIZE = 32768;
  const DWORD OUT_QUEUE_SIZE = 32768;
  SetupComm(handle, IN_QUEUE_SIZE, OUT_QUEUE_SIZE);

  // Prepare the settings structures

  GetCommConfig(handle, &commConfig, &confSize);
  GetCommState(handle, &(commConfig.dcb));

  commConfig.dcb.fBinary = TRUE;
  commConfig.dcb.fInX    = FALSE;
  commConfig.dcb.fOutX   = FALSE;
  commConfig.dcb.fAbortOnError = FALSE;
  commConfig.dcb.fNull   = FALSE;

  // ****************************************************************
  // Set the given baud rate.
  // ****************************************************************

  switch (s.baudRate) 
  {
    case 4800: 
      commConfig.dcb.BaudRate = CBR_4800; 
      break;
    case 9600: 
      commConfig.dcb.BaudRate = CBR_9600; 
      break;
    case 19200: 
      commConfig.dcb.BaudRate = CBR_19200; 
      break;
    case 38400: 
      commConfig.dcb.BaudRate = CBR_38400; 
      break;
    case 57600: 
      commConfig.dcb.BaudRate = CBR_57600; 
      break;
    case 115200: 
      commConfig.dcb.BaudRate = CBR_115200; 
      break;
    case 921600: 
      commConfig.dcb.BaudRate = 921600;   // Must use literals for the value here.
      break;
    
    default:
      // Invalid baud rate
      close();
      return false;
      break;
  }

  // ****************************************************************
  // Set the number of data bits to 8.
  // ****************************************************************
  
  commConfig.dcb.ByteSize = 8;

  // ****************************************************************
  // Set the given number of stop bits.
  // ****************************************************************

  switch (s.stopBits) 
  {
    case ONE_STOP_BIT:
      commConfig.dcb.StopBits = ONE_STOP_BIT;
      break;

    case TWO_STOP_BITS:
      commConfig.dcb.StopBits = TWO_STOP_BITS;
      break;
  }

  // ****************************************************************
  // Set the given parity type.
  // ****************************************************************
  
  switch (s.parity) 
  {
    case NO_PARITY:
      commConfig.dcb.Parity = NOPARITY;
      commConfig.dcb.fParity = FALSE;
      break;

    case EVEN_PARITY:
      commConfig.dcb.Parity = EVENPARITY;
      commConfig.dcb.fParity = TRUE;
      break;

    case ODD_PARITY:
      commConfig.dcb.Parity = ODDPARITY;
      commConfig.dcb.fParity = TRUE;
      break;
  }

  // ****************************************************************
  // Disable flow control.
  // ****************************************************************

  commConfig.dcb.fOutxCtsFlow = FALSE;
  if (s.hardwareFlowControl)
  {
    commConfig.dcb.fRtsControl = DTR_CONTROL_ENABLE;
  }
  else
  {
    commConfig.dcb.fRtsControl = RTS_CONTROL_DISABLE;
  }
  commConfig.dcb.fInX = FALSE;
  commConfig.dcb.fOutX = FALSE;

  // ****************************************************************
  // Set the timeouts such that a read function returns immediately
  // with the content of the buffer.
  // ****************************************************************

  COMMTIMEOUTS timeouts = { 0 };

  timeouts.ReadIntervalTimeout = MAXDWORD;
  timeouts.ReadTotalTimeoutConstant = 0;
  timeouts.ReadTotalTimeoutMultiplier = 0;
  timeouts.WriteTotalTimeoutConstant = 50;
  timeouts.WriteTotalTimeoutMultiplier = 10;

  if (!SetCommTimeouts(handle, &timeouts))
  {
    printf("Error in SetCommTimeouts()\n");
  }


  // ****************************************************************
  // Apply the changed configuration.
  // ****************************************************************

  return (bool)(SetCommConfig(handle, &commConfig, sizeof(COMMCONFIG)) != 0); //return (bool)SetCommConfig(handle, &commConfig, sizeof(COMMCONFIG)); // threw warning C4800: 'BOOL': forcing value to bool 'true' or 'false' (performance warning)

//  return true;
}

// ****************************************************************************
/// Closes the serial port.
// ****************************************************************************

void SerialPort::close()
{
  if (handle != INVALID_HANDLE_VALUE)
  {
    CloseHandle(handle);
  }
  handle =INVALID_HANDLE_VALUE;
}

// ****************************************************************************
// ****************************************************************************

bool SerialPort::isOpen()
{
  return (handle != INVALID_HANDLE_VALUE);
}

// ****************************************************************************
/// Returns the number of bytes available for reading in the RX buffer.
// ****************************************************************************

int SerialPort::bytesAvailable()
{
  // Determine the number of bytes in the RX buffer of the device.
  COMSTAT comStat;
  DWORD errorMask = 0;

  // Get the COM port status.
  ClearCommError(handle, &errorMask, &comStat);
  
  // The number of bytes received by the serial provider but not yet 
  // read by a ReadFile operation.
  int numBytes = comStat.cbInQue;

  return numBytes;
}


// ****************************************************************************
/// Tries to read the given number of bytes from the RX buffer. The actual 
/// number of read bytes is returned.
// ****************************************************************************

int SerialPort::readData(char *data, int numBytes)
{
  DWORD bytesRead;
  ReadFile(handle, (void*)data, numBytes, &bytesRead, NULL);

  return (int)bytesRead;
}

// ****************************************************************************
/// Writes data out on the serial port.
/// \return The actual number of bytes written.
// ****************************************************************************

int SerialPort::writeData(const char *data, int numBytes)
{
  if (isOpen() == false)
  {
    return 0;
  }

  DWORD bytesWritten;
    
  if (WriteFile(handle, (void*)data, (DWORD)numBytes, &bytesWritten, NULL)) 
  {
    // Flush the bytes
    FlushFileBuffers(handle);
    return (int)bytesWritten;
  }
  else 
  {
    return 0;
  }
}


// ****************************************************************************
