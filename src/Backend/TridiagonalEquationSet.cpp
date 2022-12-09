// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "TridiagonalEquationSet.h"
#include <string.h>


// ****************************************************************************
// Konstruktor ****************************************************************

TridiagonalEquationSet::TridiagonalEquationSet()
{
  N = 0;
}

// ****************************************************************************
// Destruktor *****************************************************************

TridiagonalEquationSet::~TridiagonalEquationSet()
{
  clearMemory();
}

// ****************************************************************************
// Allokiert den Speicher f�r die Koeff. und den L�sungsvektor ****************

void TridiagonalEquationSet::allocMemory(int numEquations)
{
  if (N == numEquations) { return; }
  clearMemory();                      // alten Speicher freigeben

  N = numEquations;
  a = new double[N];
  b = new double[N];
  c = new double[N];

  d = new double[N];
}

// ****************************************************************************
// Gibt den Speicher f�r die Koeffizienten und den L�sungsvektor frei *********

void TridiagonalEquationSet::clearMemory()
{
  if (N > 0)
  {
    delete[] a;
    delete[] b;
    delete[] c;
    delete[] d;
  }
}

// ****************************************************************************
// Setzt die Koeffizientenwerte und den L�sungsvektor *************************

void TridiagonalEquationSet::setup(const double *leftDiagonal, const double *mainDiagonal, const double *rightDiagonal, 
                                   const double *solutionVector, int numEquations)
{
  allocMemory(numEquations);

  memcpy((void*)a, (void*)leftDiagonal, N*sizeof(double));
  memcpy((void*)b, (void*)mainDiagonal, N*sizeof(double));
  memcpy((void*)c, (void*)rightDiagonal, N*sizeof(double));
  memcpy((void*)d, (void*)solutionVector, N*sizeof(double));
}

// ****************************************************************************
// L�st das Gleichungssystem (Systemmatrix �ndert sich dabei !!) **************

void TridiagonalEquationSet::solve()
{
  if (N < 1) { return; }

  int i;

  for (i=0; i < N-1; i++)
  {
    c[i]/= b[i];
    d[i]/= b[i];
    b[i+1]-= c[i]*a[i+1];
    d[i+1]-= d[i]*a[i+1];
  }
  // in der letzten Zeile nur die Division durchf�hren
  d[N-1]/= b[N-1];

  // Koeffizienten in der rechten Nebendiagonale eleminieren
  for (i=N-2; i >= 0; i--)
  {
    d[i]-= d[i+1]*c[i];
  }
}

// ****************************************************************************
// Gibt den gesuchten L�sungsvektor des LGS zur�ck ****************************

void TridiagonalEquationSet::getResult(double *result)
{
  memcpy((void*)result, (void*)d, N*sizeof(double));
}

// ****************************************************************************
