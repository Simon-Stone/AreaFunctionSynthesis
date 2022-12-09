// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __TRIDIGONAL_EQUATION_SET_H__
#define __TRIDIGONAL_EQUATION_SET_H__

#include <cmath>

// ****************************************************************************
/// This class represents a system of linear equations with a tridiagonal
/// coefficient matrix.
// ****************************************************************************

class TridiagonalEquationSet
{
  public:
    TridiagonalEquationSet();
    ~TridiagonalEquationSet();

    void setup(const double *leftDiagonal, const double *mainDiagonal, const double *rightDiagonal, 
               const double *solutionVector, int numEquations);
    void solve();
    void getResult(double *result);

  private:
    void allocMemory(int numEquations);
    void clearMemory();

    int N;        // Anzahl der Gleichungen

    // Koeffizienten um die Hauptdiagonale
    double *a;
    double *b;
    double *c;

    double *d;    // L�sungsvektor
};


#endif
