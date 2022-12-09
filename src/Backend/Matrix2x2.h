// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#ifndef __MATRIX2X2_H__
#define __MATRIX2X2_H__

#include <complex>

typedef std::complex<double> Complex;

// *********************************************************************************
// 2 x 2 matrices with complex numbers.
// *********************************************************************************

class Matrix2x2
{
  public:
    Complex A, B, C, D;

    Matrix2x2(Complex a, Complex b, Complex c, Complex d);
    Matrix2x2(void);
    void unitMatrix(void);
    void invert();
    Matrix2x2 &operator*= (const Matrix2x2 &x);
    Matrix2x2 &operator+= (const Matrix2x2 &x);
    void operator=(const Matrix2x2 &x);
};

Matrix2x2 operator*(const Matrix2x2 x, const Matrix2x2 y);
Matrix2x2 operator+(const Matrix2x2 x, const Matrix2x2 y);


#endif
