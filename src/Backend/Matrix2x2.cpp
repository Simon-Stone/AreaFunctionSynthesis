// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "Matrix2x2.h"

// Initialisierung ************************************************************

Matrix2x2::Matrix2x2(Complex a, Complex b, Complex c, Complex d)
{
  A = a; 
  B = b; 
  C = c; 
  D = d;
}

// ****************************************************************************

Matrix2x2::Matrix2x2(void)
{
  A = Complex(0, 0);
  B = Complex(0, 0);
  C = Complex(0, 0);
  D = Complex(0, 0);
}

// Einheitsmatrix generieren **************************************************

void Matrix2x2::unitMatrix(void)
{
  A = Complex(1, 0);
  B = Complex(0, 0);
  C = Complex(0, 0);
  D = Complex(1, 0);
}

// Invertiert die Matrix (Vorraussetzung ist, dass die Determinante = 1) ******

void Matrix2x2::invert()
{
  Complex t;

  t = A;
  A = D;
  D = t;
  B = -B;
  C = -C;
}

// Addition *******************************************************************

Matrix2x2 &Matrix2x2::operator+=(const Matrix2x2 &x)
{
  A+= x.A;
  B+= x.B;
  C+= x.C;
  D+= x.D;
  return(*this);
}

// ****************************************************************************

Matrix2x2 operator+(const Matrix2x2 x, const Matrix2x2 y)
{
  Matrix2x2 z = x;
  return(z+=y);
}

// Multiplikation *************************************************************

Matrix2x2 &Matrix2x2::operator*=(const Matrix2x2 &x)
{
  Complex a, b, c, d;

  a = A*x.A + B*x.C;
  b = A*x.B + B*x.D;
  c = C*x.A + D*x.C;
  d = C*x.B + D*x.D;

  A = a; B = b; C = c; D = d;
  
  return(*this);
}

// ****************************************************************************

Matrix2x2 operator*(const Matrix2x2 x, const Matrix2x2 y)
{
  Matrix2x2 z = x;
  return(z*=y);
}

// Zuweisung ******************************************************************

void Matrix2x2::operator =(const Matrix2x2 &x)
{
  A = x.A;
  B = x.B;
  C = x.C;
  D = x.D;
}

// ****************************************************************************
