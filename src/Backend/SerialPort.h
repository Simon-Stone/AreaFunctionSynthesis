#ifndef _SERIAL_PORT_H_
#define _SERIAL_PORT_H_

#include <windows.h>

// ****************************************************************************
/// \brief This class represents the serial COM port used to communicate with the EOS
/// hardware.
// ****************************************************************************

class SerialPort
{
  // **************************************************************************
  // Public data.
  // **************************************************************************

public:
  enum Parity
  {
    NO_PARITY,
    ODD_PARITY,
    EVEN_PARITY
  };

  enum StopBits
  {
    ONE_STOP_BIT,
    TWO_STOP_BITS
  };
	/// \brief This struct represents the settings of the EOS COM port.
  struct PortSettings 
  {
    char portName[256];
    int baudRate;
    Parity parity;
    StopBits stopBits;
    bool hardwareFlowControl;
    unsigned long timeout_s;
    unsigned long timeout_ms;
};

  // **************************************************************************
  // **************************************************************************

public:
  SerialPort();
  ~SerialPort();

  bool open(const char *portName, int baudRate, Parity parity = NO_PARITY, 
    StopBits stopBits = ONE_STOP_BIT, bool hardwareFlowControl = false);
  bool open(PortSettings s);
  void close();
  bool isOpen();

  int bytesAvailable();
  int readData(char *data, int numBytes);
  int writeData(const char *data, int numBytes);

  // **************************************************************************
  // Private data.
  // **************************************************************************

private:
  PortSettings settings;

  // Windows specific variables
  HANDLE handle;
  COMMCONFIG commConfig;
  COMMTIMEOUTS commTimeouts;

};

#endif