// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "IirFilter.h"
#include "Dsp.h"


// ****************************************************************************
/// Constructor. 
// ****************************************************************************


IirFilter::IirFilter()
{
  resetBuffers();
  createUnityFilter();
}

// ****************************************************************************
/// Resets the buffers for the input and output signals.
/// \param initialValue The initial valueto be assumed at the beginning of the
/// signal: set it either to zero or to the first sample value of the input
/// signal.
// ****************************************************************************

void IirFilter::resetBuffers(double initialValue)
{
  int i;
  for (i=0; i < IIR_BUFFER_LENGTH; i++)
  {
    inputBuffer[i] = initialValue;
    outputBuffer[i] = initialValue;
  }

  pos = 0;
}


// ****************************************************************************
/// Calculates the next sample of the output signal.
// ****************************************************************************
 
double IirFilter::getOutputSample(double nextInputSample)
{
  int i;

  inputBuffer[pos & IIR_BUFFER_MASK] = nextInputSample;

  double sum = a[0]*nextInputSample;

  for (i=1; i <= order; i++)
  {
    sum+= a[i]*inputBuffer[(pos-i) & IIR_BUFFER_MASK];
    sum+= b[i]*outputBuffer[(pos-i) & IIR_BUFFER_MASK];
  }

  outputBuffer[pos & IIR_BUFFER_MASK] = sum;
  pos++;

  return sum;
}

// ******************************************************************************
// Gibt den komplexen Wert der Übertragungsfunktion zurück, wobei freqRatio=f/SR
    
ComplexValue IirFilter::getFrequencyResponse(double freqRatio)
{
  int i;
  ComplexValue z = std::exp(ComplexValue(0.0, 2.0*M_PI*freqRatio));
  ComplexValue factor = 1.0;

  ComplexValue numerator = a[0];
  ComplexValue denominator = 1.0;

  for (i=1; i <= order; i++)
  {
    factor/= z;
    numerator+= a[i]*factor;
    denominator-= b[i]*factor;
  }

  return numerator / denominator;
}

// ****************************************************************************
// Gibt den Frequenzgang des Filters zurück. **********************************

void IirFilter::getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength)
{
  int i;

  spectrum->reset(spectrumLength);

  for (i=0; i <= spectrumLength/2; i++)
  {
    spectrum->setValue(i, getFrequencyResponse((double)i / (double)spectrumLength));
  }

  generateNegativeFrequencies(spectrum);
}

// ****************************************************************************
// Gibt den Frequenzgang des Filters zurück. **********************************

void IirFilter::getFrequencyResponse(ComplexSignal *spectrum, int spectrumLength, int SR, double F0)
{
  int i;

  spectrum->reset(spectrumLength);

  for (i=0; i <= spectrumLength/2; i++)
  {
    spectrum->setValue(i, getFrequencyResponse((double)i*F0 / (double)SR));
  }

  generateNegativeFrequencies(spectrum);
}


// ****************************************************************************
// Multipliziert alle a-Koeffizienten mit gain. *******************************

void IirFilter::setGain(double gain)
{
  int i;

  for (i=0; i <= order; i++) { a[i]*= gain; }
}

// ****************************************************************************
// Setzt die Filterkoeffizienten auf die übergebenen Werte ********************
    
void IirFilter::setCoefficients(const double *A, const double *B, const int newOrder)
{
  int i;

  clearCoefficients();
  order = newOrder;
  if (order > MAX_IIR_ORDER) { order = MAX_IIR_ORDER; }

  for (i=0; i <= order; i++)
  {
    a[i] = A[i];
    b[i] = B[i];
  }
}

// ****************************************************************************
// Kombiniert diesen Filter mit dem übergebenen Filter als Reihenschaltung oder
// Parallelschaltung. *********************************************************

bool IirFilter::combineWithFilter(const IirFilter *f, bool cascade)
{
  if ((order > MAX_IIR_ORDER/2) || (f->order > MAX_IIR_ORDER/2)) { return false; }

  double a2[MAX_IIR_ORDER+1];
  double b2[MAX_IIR_ORDER+1];
  double a3[MAX_IIR_ORDER+1];
  double b3[MAX_IIR_ORDER+1];
  int i, j;

  // Umwandlung der Rekursionskoeffizienten in eine Übertragungsfunktion
  for (i=0; i <= MAX_IIR_ORDER; i++)
  {
    a2[i] = f->a[i];
    b2[i] = f->b[i];

    b[i]  = -b[i];
    b2[i] = -b2[i];
  }
  b[0] = 1.0;
  b2[0] = 1.0;

  // Multiplikation der Polynome durch Faltung **********************

  for (i=0; i <= MAX_IIR_ORDER; i++)
  {
    a3[i] = 0.0;
    b3[i] = 0.0;

    for (j=0; j <= MAX_IIR_ORDER/2; j++)
    {
      if ((i-j >= 0) && (i-j <= MAX_IIR_ORDER/2))
      {
        if (cascade)
        {
          a3[i]+= a[j]*a2[i-j];
        }
        else
        {
          a3[i]+= a[j]*b2[i-j] + a2[j]*b[i-j];
        }

        b3[i]+= b[j]*b2[i-j];
      }
    }
  }

  // Dem aktuellen Filter die neuen Rekursionskoeffizienten zuweisen 

  for (i=0; i <= MAX_IIR_ORDER; i++)
  {
    a[i] = a3[i];
    b[i] = b3[i];
  }

  // Übertragungsfunktion in Rekursionskoeffizienten verwandeln
  
  for (i=0; i <= MAX_IIR_ORDER; i++) { b[i] = -b[i]; }
  b[0] = 0.0;

  order+= f->order;

  return true;
}

// ****************************************************************************
// Erzeugt einen Filter mit der konstanten Übertragungsfunktion H(z) = 1 ******
    
void IirFilter::createUnityFilter()
{
  clearCoefficients();
  a[0] = 1.0;
  b[0] = 1.0;
  order = 0;
}

// ****************************************************************************
// ****************************************************************************

void IirFilter::createSinglePoleLowpass(double cutoffFreqRatio)
{
  clearCoefficients();
  double x = exp(-2.0*M_PI*cutoffFreqRatio);

  order = 1;
  a[0] = 1.0 - x;
  b[1] = x;
}

// ****************************************************************************
// ****************************************************************************

void IirFilter::createSinglePoleHighpass(double cutoffFreqRatio)
{
  clearCoefficients();
  double x = exp(-2.0*M_PI*cutoffFreqRatio);

  order = 1;
  a[0] = 0.5*(1.0 + x);
  a[1] = -a[0];
  b[1] = x;
}

// ****************************************************************************
// Ein konjugiert komplexes Polpaar mit der Polgüte Q. Transformation mit der  
// Bilinearen z-Transformation.
// ****************************************************************************

void IirFilter::createSecondOrderLowpass(double freqRatio, double Q)
{
  clearCoefficients();

  order = 2;

  double K = tan(M_PI*freqRatio);
  if (Q == 0.0) { Q = 0.000001; }
  double denominator = K*K + K/Q + 1.0;
  
  b[0] = 1.0;
  b[1] = -2.0*(K*K-1.0) / denominator;
  b[2] = -(K*K - K/Q + 1.0) / denominator;

  a[0] = K*K / denominator;
  a[1] = 2.0*a[0];
  a[2] = a[0];
}

// ****************************************************************************
// Berechnet die Koeffizienten für einen Chebyshev-Filter (Hoch- oder Tiefpass)

void IirFilter::createChebyshev(double cutoffFreqRatio, bool isHighpass, int numPoles)
{
  const double percentRipple = 0.5;    // in %
  double ta[MAX_IIR_ORDER+1];
  double tb[MAX_IIR_ORDER+1];
  double a0, a1, a2, b1, b2;    // Koeff., die für jede 2-Pol-Stufe des Filters berechnet werden
  double t, w, m, d, k;
  double x0, x1, x2, y1, y2;
  double kx, vx;
  double temp;
  double polReal, polImag;
  double es;
  int i;

  // Die Anzahl der Pole muss gerade sein ***************************

  if ((numPoles & 1) != 0)
  {
    numPoles++;
  }
  if (numPoles > MAX_IIR_ORDER) { numPoles = MAX_IIR_ORDER; }
  
  order = numPoles;     // Die Ordnung des Filters in die Klassenglobale Variable übernehmen

  // Koeffizienten initialisieren ***********************************

  for (i=0; i <= MAX_IIR_ORDER; i++)
  {
    a[i] = 0.0;
    b[i] = 0.0;
  }
  a[2] = 1.0;
  b[2] = 1.0;

  // ****************************************************************
  // numPoles/2 mal die Hauptschleife durchlaufen *******************

  int pol;

  for (pol=1; pol <= numPoles/2; pol++)
  {
    
    // Unterprogramm, welches für jede 2-Pol-Stufe des ****
    // Filters ausgeführt wird ****************************
    
    // Polposition auf dem Einheitskreis berechnen
    polReal = -cos(M_PI/(2.0*numPoles) + (M_PI*(pol-1)) / (double)numPoles);
    polImag =  sin(M_PI/(2.0*numPoles) + (M_PI*(pol-1)) / (double)numPoles);

    if (percentRipple != 0.0)
    {
      // von einem Kreis zu einer Ellipse warpen
      temp = 100.0 / (100.0 - percentRipple);
      es = sqrt(temp*temp - 1.0);
      
      vx = (1.0/(double)numPoles)*log((1.0/es) + sqrt((1.0/(es*es)) + 1));
      kx = (1.0/(double)numPoles)*log((1.0/es) + sqrt((1.0/(es*es)) - 1));
      kx = 0.5*(exp(kx) + exp(-kx));

      polReal = polReal*(0.5*(exp(vx) - exp(-vx))) / kx;
      polImag = polImag*(0.5*(exp(vx) + exp(-vx))) / kx;
    }

    // Überführung vom s- in den z-Bereich
    t = 2.0*tan(0.5);
    w = 2.0*M_PI*cutoffFreqRatio;
    m = polReal*polReal + polImag*polImag;
    d = 4.0 - 4.0*polReal*t + m*t*t;
    
    x0 = (t*t) / d;
    x1 = (2.0*t*t) / d;
    x2 = (t*t) / d;
    y1 = (8.0 - 2.0*m*t*t) / d;
    y2 = (-4.0 - 4.0*polReal*t - m*t*t) / d;

    // lowpass -> lowpass oder lowpass -> highpass Transformation
    if (isHighpass) 
      { k = -cos(0.5*w + 0.5) / cos(0.5*w - 0.5); }
    else
      { k =  sin(0.5 - 0.5*w) / sin(0.5 + 0.5*w); }
    
    d = 1.0 + y1*k - y2*k*k;
    a0 = (x0 - x1*k + x2*k*k) / d;
    a1 = (-2.0*x0*k + x1 + x1*k*k - 2.0*x2*k) / d;
    a2 = (x0*k*k - x1*k + x2) / d;
    b1 = (2.0*k + y1 + y1*k*k - 2.0*y2*k) / d;
    b2 = (-(k*k) - y1*k + y2) / d;

    if (isHighpass) { a1 = -a1; b1 = -b1; }

    // Ende der Subroutine ********************************
    // ****************************************************

    // Die berechneten Koeffizienten zur Kaskade addieren
    for (i=0; i <= MAX_IIR_ORDER; i++)
    {
      ta[i] = a[i];
      tb[i] = b[i];
    }

    for (i=2; i <= MAX_IIR_ORDER; i++)
    {
      a[i] = a0*ta[i] + a1*ta[i-1] + a2*ta[i-2];
      b[i] =    tb[i] - b1*tb[i-1] - b2*tb[i-2];
    }
  }

  
  // ******************************************************
  // Die Kombination der Koeffizienten abschliessen *******

  b[2] = 0.0;
  for (i=0; i <= MAX_IIR_ORDER-2; i++)
  {
    a[i] =  a[i+2];
    b[i] = -b[i+2];
  }

  // Die Werte normalisieren ******************************
  double sumA = 0.0;
  double sumB = 0.0;

  for (i=0; i <= MAX_IIR_ORDER-2; i++)
  {
    if (!isHighpass) 
    { 
      sumA+= a[i]; 
      sumB+= b[i]; 
    }
    else
    {
      if ((i & 1) == 0)
      {
        sumA+= a[i];
        sumB+= b[i];
      }
      else
      {
        sumA-= a[i];
        sumB-= b[i];
      }
    }
  }

  double gain = sumA / (1.0 - sumB);
  for (i=0; i <= MAX_IIR_ORDER-2; i++) { a[i]/= gain; }
}


// ****************************************************************************
// Setzt alle Filterkoeffizienten (bis auf b0) auf Null. **********************

void IirFilter::clearCoefficients()
{
  int i;

  order = 0;
  for (i=0; i <= MAX_IIR_ORDER; i++)
  {
    a[i] = 0.0;
    b[i] = 0.0;
  }

  b[0] = 1;
}

// ****************************************************************************
