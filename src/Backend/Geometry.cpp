// ****************************************************************************
// This file is part of VocalTractLab.
// Copyright (C) 2008, Peter Birkholz, Hamburg, Germany
// www.vocaltractlab.de
// author: Peter Birkholz
// ****************************************************************************

#include "Geometry.h"
#include <cmath>

// ----------------------------------------------------------------------------
// 2D-Punkt
// ----------------------------------------------------------------------------

// ****************************************************************************

bool Point2D::isRightFrom(Vector2D& V)
{
  return (V.v.x*(V.P.y - y) - V.v.y*(V.P.x - x)) >= 0;
}

// ****************************************************************************

bool Point2D::isLeftFrom(Vector2D& V)
{
  return (V.v.x*(V.P.y - y) - V.v.y*(V.P.x - x)) <= 0;
}

// ****************************************************************************

Point3D Point2D::toPoint3D()
{
  return Point3D(x, y, 0.0);
}

// ****************************************************************************
// Die Strecke vom Punkt A (Anker) zu diesem Punkt soll mit diesem Punkt an die
// Gerade "angelehnt" werden, die durch den Vektor V beschrieben wird. ********

void Point2D::leanOn(Vector2D V, Point2D A)
{
  Point2D S[2];         // Die 2 möglichen neuen Punkte
  Point2D s[2];         // Die Vektoren von A aus zu diesen neuen Punkten
  Point2D  w = V.P - A;
  Point2D& v = V.v;

  Point2D r(x - A.x, y - A.y);          // Der Vektor der "alten" Linie

  // c ist der Abstand zw. diesem Punkt und A zum quadrat *
  double c = r.x*r.x + r.y*r.y;
  
  double denominator = v.x*v.x + v.y*v.y;
  if (denominator == 0.0) { denominator = 0.0; }

  double p = 2.0*(v.x*w.x + v.y*w.y) / denominator;
  double q = (w.x*w.x + w.y*w.y - c) / denominator;

  double radicant = 0.25*p*p - q;

  // Ist der Wert unter der Wurzel > 0 ? ******************

  if (radicant >= 0.0)
  {
    double root = sqrt(radicant);
    S[0] = V.P + v*(-0.5*p + root);
    S[1] = V.P + v*(-0.5*p - root);

    // Welcher der beiden Ergebnispunkte S[0] und S[1] wird ausgewählt ? ******
    s[0] = S[0] - A;
    s[1] = S[1] - A;

    if (s[0].x*r.x + s[0].y*r.y > s[1].x*r.x + s[1].y*r.y)
    {
      x = S[0].x;
      y = S[0].y;
    }
    else
    {
      x = S[1].x;
      y = S[1].y;
    }
    
  }
  
}

// ****************************************************************************

Point2D Point2D::normalize()
{
  double l = sqrt(x*x + y*y);
  if (l != 0.0)
  {
    x/= l;
    y/= l;
  }
  return Point2D(x, y);
}

// ****************************************************************************
// Dreht den Vektor (x,y) um 90° nach rechts **********************************

Point2D Point2D::turnRight()
{
  double temp = x;
  x = y;
  y = -temp;
  return Point2D(x, y);
}

// Dreht den Vektor (x,y) um 90° nach links ***********************************
// ****************************************************************************

Point2D Point2D::turnLeft()
{
  double temp = x;
  x = -y;
  y = temp;
  return Point2D(x, y);
}

// ****************************************************************************
// Dreht den Vektor um angle rad im math. pos. Sinne. *************************

void Point2D::turn(double angle)
{
  double S = sin(angle);
  double C = cos(angle);
  double nx = x*C - y*S;
  double ny = y*C + x*S;
  x = nx;
  y = ny;
}

// ****************************************************************************
// Gibt den Abstand dieses Punktes vom Vektor V zurück. Bei negativen Abständen
// befindet sich der Punkt links von V, sonst rechts. *************************

double Point2D::getDistanceFrom(Vector2D V)
{
  Point2D w(V.P.x - x, V.P.y - y);
  Point2D &v = V.v;
  double denominator = v.x*v.x + v.y*v.y;

  if (denominator == 0.0) { denominator = 0.0001; }

  return (v.x*w.y - v.y*w.x) / denominator;
}

// ****************************************************************************
// Gibt das Skalarprodukt zweier Punkte zurück. *******************************

double scalarProduct(const Point2D &P, const Point2D &Q)
{ 
  return P.x*Q.x + P.y*Q.y; 
}


// ----------------------------------------------------------------------------
// 3D-Punkt
// ----------------------------------------------------------------------------

// ****************************************************************************

bool Point3D::isRightFrom(Vector2D& V)
{
  return (V.v.x*(V.P.y - y) - V.v.y*(V.P.x - x)) >= 0;
}

// ****************************************************************************

Point2D Point3D::toPoint2D()
{
  return Point2D(x, y);
}

// ****************************************************************************

Point3D Point3D::normalize()
{
  double l = sqrt(x*x + y*y + z*z);
  if (l != 0.0)
  {
    x/= l;
    y/= l;
    z/= l;
  }
  return Point3D(x, y, z);
}

// ****************************************************************************

double scalarProduct(const Point3D &P, const Point3D &Q)
{
  return (P.x*Q.x + P.y*Q.y + P.z*Q.z);

}

// ****************************************************************************

Point3D crossProduct(const Point3D &P, const Point3D &Q)
{
  return Point3D(P.y*Q.z - P.z*Q.y, P.z*Q.x - P.x*Q.z, P.x*Q.y - P.y*Q.x);
}
 
// ----------------------------------------------------------------------------
// 2D-Vektor
// ----------------------------------------------------------------------------

// ****************************************************************************
// Konstruktor ****************************************************************

Vector2D::Vector2D()
{
  P.set(0.0, 0.0);
  v.set(0.0, 0.0);
}

// ****************************************************************************
// Konstruktor ****************************************************************

Vector2D::Vector2D(Point2D _P, Point2D _v)
{
  P = _P;
  v = _v;
}

// ****************************************************************************
// Setzt den Ankerpunkt und den Vektoranteil **********************************

void Vector2D::set(Point2D _P, Point2D _v)
{
  P = _P;
  v = _v;
}

// ****************************************************************************
// Gibt den Schnittpunkt des Vektors mit dem Vektor V zurück ******************

Point2D Vector2D::getIntersection(Vector2D V, double& t)
{
  Point2D& B = V.P;
  Point2D  w = V.v;

  double denominator = v.x*w.y - v.y*w.x;
  if (denominator == 0.0) { denominator = 0.0001; }

  t = (w.x*(P.y - B.y) - w.y*(P.x - B.x)) / denominator;
  return getPoint(t);
}

// ****************************************************************************
// Gibt den Schnittpunkt des Vektors mit der Linie L zurück *******************

Point2D Vector2D::getIntersection(Line2D L, double& t, bool& ok)
{
  Point2D& B = L.P[0];
  Point2D  w = L.P[1] - L.P[0];

  double denominator = v.x*w.y - v.y*w.x;
  if (denominator == 0.0) { denominator = 0.0001; }

  // Liegt der Schnittpunkt auch AUF der Linie ?
  double m = (v.x*(P.y - B.y) - v.y*(P.x - B.x)) / denominator;
  if ((m < -0.01) || (m > 1.01) || (L.P[0] == L.P[1]))
  {
    t = 0.0;
    ok = false;
    return P;
  }

  ok = true;
  t = (w.x*(P.y - B.y) - w.y*(P.x - B.x)) / denominator;
  return getPoint(t);
}

// ****************************************************************************
// Gibt den Schnittpunkt der 3D-Linie mit der Fläche zurück, die durch diesen *
// Vektor und den Vektor (0, 0, -1) aufgespannt wird **************************

Point3D Vector2D::getIntersection(Line3D L, double& t, bool& ok)
{
  Point3D& B = L.P[0];
  Point3D  w = L.P[1] - L.P[0];

  double denominator = v.x*w.y - v.y*w.x;
  if (denominator == 0.0) { denominator = 0.0001; }

  // Liegt der Schnittpunkt auch AUF der Linie ?
  double m = (v.x*(P.y - B.y) - v.y*(P.x - B.x)) / denominator;
  if ((m < -0.01) || (m > 1.01) || ((w.x == 0) && (w.y == 0)))
  {
    t = 0.0;
    ok = false;
    return Point3D(P.x, P.y, 0.0);
  }

  ok = true;
  t = (w.x*(P.y - B.y) - w.y*(P.x - B.x)) / denominator;
  return L.getPoint(m);
}

// ****************************************************************************
// Normalisiert den Vektoranteil **********************************************

void Vector2D::normalize()
{
  double argument = v.x*v.x + v.y*v.y;
  if ((argument != 1.0) && (argument != 0.0))
  {
    argument = sqrt(argument);
    v.x/= argument;
    v.y/= argument;
  }
}

// ****************************************************************************
// Gibt den PUnkt an P+n*v zurück *********************************************

Point2D Vector2D::getPoint(double t)
{
  return P + v*t;
}

// ****************************************************************************
// Gibt die Länge des mit t gestreckten Vektoranteils zurück ******************

double Vector2D::getLength(double t)
{
  return sqrt(v.x*v.x + v.y*v.y)*t;
}

// ****************************************************************************
// Gibt true zurück, wenn der Vektoranteil ungleich dem Nullvektor ist ********

bool Vector2D::isNotNull()
{
  return (v.x != 0.0) || (v.y != 0.0);
}

// ****************************************************************************
// Zuweisung ******************************************************************

void Vector2D::operator=(Vector2D Q)
{
  P = Q.P;
  v = Q.v;  
}


// ----------------------------------------------------------------------------
// 2D-Linie 
// ----------------------------------------------------------------------------

// ****************************************************************************
// Konstruktor ****************************************************************

Line2D::Line2D()
{
  P[0].set(0.0, 0.0);
  P[1].set(0.0, 0.0);
}

// ****************************************************************************
// Konstruktor ****************************************************************

Line2D::Line2D(Point2D P0, Point2D P1)
{
  P[0] = P0;
  P[1] = P1;
}

// ****************************************************************************
// Setzt die Endpunkte ********************************************************

void Line2D::set(Point2D P0, Point2D P1)
{
  P[0] = P0;
  P[1] = P1;
}

// ****************************************************************************
// Gibt den Schnittpunkt dieser Linie mit dem Vektor V zurück *****************

Point2D Line2D::getIntersection(Vector2D V, double& t, bool& ok)
{
  Point2D& A = P[0];
  Point2D  v = P[1] - P[0];
  Point2D& B = V.P;
  Point2D  w = V.v;
  
  double origDen = v.x*w.y - v.y*w.x;
  double denominator = origDen;
  if (denominator == 0.0) { denominator = 0.0001; }

  t = (w.x*(A.y - B.y) - w.y*(A.x - B.x)) / denominator;
  if ((t > -0.01) && (t < 1.01) && (origDen != 0)) { ok = true; } else { ok = false; }
  
  return getPoint(t);
}

// ****************************************************************************
// Gibt den Schnittpunkt dieser Linie mit der Linie L zurück ******************

Point2D Line2D::getIntersection(Line2D L, double& t, bool& ok)
{
  Point2D& A = P[0];
  Point2D  v = P[1] - P[0];
  Point2D& B = L.P[0];
  Point2D  w = L.P[1] - L.P[0];
  
  double origDen = v.x*w.y - v.y*w.x;
  double denominator = origDen;
  if (denominator == 0.0) { denominator = 0.0001; }

  t = (w.x*(A.y - B.y) - w.y*(A.x - B.x)) / denominator;
  double m = (v.x*(A.y - B.y) - v.y*(A.x - B.x)) / denominator;

  if ((t > -0.01) && (t < 1.01) && 
      (m > -0.01) && (m < 1.01) && (origDen != 0.0)) { ok = true; } else { ok = false; }
  
  return getPoint(t);
}

// ****************************************************************************
// Gibt den Punkt an der Stelle t in[0..1] auf der Linie zurück ***************

Point2D Line2D::getPoint(double t)
{
  return P[0] + (P[1]-P[0])*t;
}

// ****************************************************************************
// Gibt die Länge der Linie zurück ********************************************

double Line2D::getLength()
{
  return sqrt((P[1].x-P[0].x)*(P[1].x-P[0].x) + (P[1].y-P[0].y)*(P[1].y-P[0].y));
}

// ****************************************************************************
// Wird Q von dem die Linie umgebenden Rechteck eingeschlossen? ***************

bool Line2D::encloses(Point2D Q)
{
  double top, bottom;
  double left, right;

  if (P[0].x < P[1].x)
  {
    left = P[0].x;
    right = P[1].x;
  }
  else
  {
    left = P[1].x;
    right = P[0].x;
  }

  if (P[0].y < P[1].y)
  {
    bottom = P[0].y;
    top = P[1].y;
  }
  else
  {
    bottom = P[1].y;
    top = P[0].y;
  }

  return ((left-0.01 < Q.x) && (right+0.01 > Q.x) && (bottom-0.01 < Q.y) && (top+0.01 > Q.y));
}


// ----------------------------------------------------------------------------
// 3D-Linie 
// ----------------------------------------------------------------------------

// ****************************************************************************

Line3D::Line3D()
{
  P[0].set(0.0, 0.0, 0.0);
  P[1].set(0.0, 0.0, 0.0);
}

// ****************************************************************************

Line3D::Line3D(Point3D P0, Point3D P1)
{
  P[0] = P0;
  P[1] = P1;
}

// ****************************************************************************

void Line3D::set(Point3D P0, Point3D P1)
{
  P[0] = P0;
  P[1] = P1;
}

// ****************************************************************************
// Gibt den Schnittpunkt der Linie mit der Ebene zurück, die durch den Vektor *
// V in ihrer x-y-Ausrichtung beschrieben wird, und in z-Richtung ragt ********

Point3D Line3D::getIntersection(Vector2D V, double& t, bool& ok)
{
  // Das Ganze auf eine 2D-Schnittpunktsberechnung zurückführen...
  
  Line2D L(Point2D(P[0].x, P[0].y), Point2D(P[1].x, P[1].y));
  Point2D Q = L.getIntersection(V, t, ok);

  return getPoint(t);
}

// ****************************************************************************

Point3D Line3D::getPoint(double t)
{
  return P[0] + (P[1] - P[0])*t;
}

// ****************************************************************************

double Line3D::getLength()
{
  Point3D d = P[1] - P[0];
  return sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}


// ----------------------------------------------------------------------------
// Kreis    
// ----------------------------------------------------------------------------

// Konstruktor ********************************************
// ********************************************************

Circle::Circle()
{
  M.x = 0.0;
  M.y = 0.0;
  r = 1.0;

  arcAngle[0] = 0.0;
  arcAngle[1] = 0.0;
}

// ********************************************************

Circle::Circle(Point2D C, double radius)
{
  M.x = C.x;
  M.y = C.y;
  r = radius;

  arcAngle[0] = 0.0;
  arcAngle[1] = 0.0;
}

// ********************************************************
// Gibt den gültigen Kreisbogenabschnitt durch 2 Winkel an.
// Gültig ist der Abschnitt VON angle0 BIS angle1. ********

void Circle::setValidArc(double angle0, double angle1)
{
  // arcAngle[0] und [1] sollen sich im Bereich 0..2*Pi bewegen
  
  while (angle0 < 0.0) { angle0+= 2.0*M_PI; }
  while (angle0 > 2.0*M_PI) { angle0-= 2.0*M_PI; }
  arcAngle[0] = angle0;

  while (angle1 < 0.0) { angle1+= 2.0*M_PI; }
  while (angle1 > 2.0*M_PI) { angle1-= 2.0*M_PI; }
  arcAngle[1] = angle1;
}

// Gibt die Länge des gültigen Kreisbogens zurück. ********
// ********************************************************

double Circle::getLength()
{
  return getLength(arcAngle[0], arcAngle[1]);
}

// Gibt die Bogenlänge vom Winkel a0 zum Winkel a1 zurück. 
// (mathematisch positiv gesehen) *************************

double Circle::getLength(double a0, double a1)
{
  double l = 0.0;

  while (a0 < 0.0) { a0+= 2.0*M_PI; }
  while (a0 > 2.0*M_PI) { a0-= 2.0*M_PI; }

  while (a1 < 0.0) { a1+= 2.0*M_PI; }
  while (a1 > 2.0*M_PI) { a1-= 2.0*M_PI; }
  
  if (a0 == a1)
    { l = 2.0*M_PI*r; }
  else
  {
    if (a1 > a0)
      { l = r*(a1 - a0); }
    else
      { l = (a1 + (2.0*M_PI - a0))*r; }
  }

  return l;
}

// ********************************************************
// ********************************************************

Point2D Circle::getPoint(double angle)
{
  Point2D P;
  P.x = M.x + r*cos(angle);
  P.y = M.y + r*sin(angle);
  return P;
}

// ****************************************************************************
// Gibt die Normale des Kreises beim Winkel angle zurück **********************

Vector2D Circle::getNormal(double angle)
{
  double C = cos(angle);
  double S = sin(angle);
  Point2D P(M.x + r*C, M.y + r*S);
  Point2D v(C, S);
  return Vector2D(P, v);
}

// ****************************************************************************
// Gibt den Schnittpunkt und den Winkel zurück, an denen **********************
// der Vektor V den gültigen Kreisbogen schneidet. ****************************

Point2D Circle::getIntersection(Vector2D V, double &intersectionAngle, bool &ok)
{
  intersectionAngle = 0.0;
  Point2D intersectionPoint(0.0, 0.0);
  Point2D& P = V.P;
  Point2D& v = V.v;
  
  double denominator = v.x*v.x + v.y*v.y;
  if (denominator == 0.0) { denominator = 0.0001; }

  double p = (2.0*v.x*(P.x - M.x) + 2.0*v.y*(P.y - M.y)) / denominator;
  double q = ((P.x-M.x)*(P.x-M.x) + (P.y-M.y)*(P.y-M.y) - r*r) / denominator;

  double radicant = 0.25*p*p - q;

  // Die Gerade schneidet den Kreis nicht *****************
  
  if (radicant < 0.0)
  {
    ok = false;
  }
  else

  // Der Kreis wird geschnitten ***************************
  {
    ok = true;

    double n[2];
    double root = sqrt(radicant);
    
    n[0] = -0.5*p + root;
    n[1] = -0.5*p - root;
    
    Point2D S[2];       // Die 2 Schnittpunkte

    S[0] = P + v*n[0];
    S[1] = P + v*n[1];

    double angle[2];

    // Die Winkel der beiden Schnittpunkte für den Bereich [0 .. 2*Pi]
    angle[0] = atan2(S[0].y - M.y, S[0].x - M.x);
    if (angle[0] < 0.0) { angle[0]+= 2.0*M_PI; }

    angle[1] = atan2(S[1].y - M.y, S[1].x - M.x);
    if (angle[1] < 0.0) { angle[1]+= 2.0*M_PI; }

    // Liegen die Schnittpunkte im gültigen Bogenabschnitt ?
    bool isInArc[2];

    if (arcAngle[0] == arcAngle[1])
    {
      isInArc[0] = isInArc[1] = true;
    }
    else
    if (arcAngle[1] > arcAngle[0])
    {
      if ((angle[0] >= arcAngle[0]) && (angle[0] <= arcAngle[1])) { isInArc[0] = true; } else { isInArc[0] = false; }
      if ((angle[1] >= arcAngle[0]) && (angle[1] <= arcAngle[1])) { isInArc[1] = true; } else { isInArc[1] = false; }
    }
    else
    {
      if ((angle[0] <= arcAngle[1]) || (angle[0] >= arcAngle[0])) { isInArc[0] = true; } else { isInArc[0] = false; }
      if ((angle[1] <= arcAngle[1]) || (angle[1] >= arcAngle[0])) { isInArc[1] = true; } else { isInArc[1] = false; }
    }

    // Sowohl der 1. als auch der 2. Schnittpunkt liegen im gültigen Kreisabschnitt.
    // Den dem Punkt P an nächsten gelegenen Schnittpunkt auswählen!
    if ((isInArc[0]) && (isInArc[1]))
    {
      if ((S[0].x-P.x)*(S[0].x-P.x) + (S[0].y-P.y)*(S[0].y-P.y) < (S[1].x-P.x)*(S[1].x-P.x) + (S[1].y-P.y)*(S[1].y-P.y))
      {
        intersectionAngle = angle[0];
        intersectionPoint = S[0];
      }
      else
      {
        intersectionAngle = angle[1];
        intersectionPoint = S[1];
      }
    }
    else

    // Nur der 1. Schnittpunkt liegt im gültigen Kreisabschnitt
    if (isInArc[0])
    {
      intersectionAngle = angle[0];
      intersectionPoint = S[0];
    }
    else

    // Nur der 2. Schnittpunkt liegt im gültigen Kreisabschnitt
    if (isInArc[1])
    {
      intersectionAngle = angle[1];
      intersectionPoint = S[1];
    }
    else

    // Kein Schnittpunkt liegt im gültigen Kreisabschnitt
    {
      ok = false;
    }
    
  }

  return intersectionPoint;
}


// ****************************************************************************
// Gibt den Schnittpunkt und den Winkel zurück, an denen **********************
// die Linie V den gültigen Kreisbogen schneidet. *****************************

Point2D Circle::getIntersection(Line2D L, double &intersectionAngle, bool &ok)
{
  intersectionAngle = 0.0;
  Point2D intersectionPoint(0.0, 0.0);
  Point2D& P = L.P[0];
  Point2D  v = L.P[1] - L.P[0];
  
  double denominator = v.x*v.x + v.y*v.y;
  if (denominator == 0.0) { denominator = 0.0001; }

  double p = (2.0*v.x*(P.x - M.x) + 2.0*v.y*(P.y - M.y)) / denominator;
  double q = ((P.x-M.x)*(P.x-M.x) + (P.y-M.y)*(P.y-M.y) - r*r) / denominator;

  double radicant = 0.25*p*p - q;

  // Die Gerade schneidet den Kreis nicht *****************
  
  if (radicant < 0.0)
  {
    ok = false;
  }
  else

  // Der Kreis wird geschnitten ***************************
  {
    ok = true;

    double n[2];
    double root = sqrt(radicant);
    
    n[0] = -0.5*p + root;
    n[1] = -0.5*p - root;
    
    Point2D S[2];       // Die 2 Schnittpunkte
    double angle[2];
    bool isInArc[2];
    bool isInLine[2];
    bool isOk[2];

    S[0] = P + v*n[0];
    S[1] = P + v*n[1];

    // Die Winkel der beiden Schnittpunkte für den Bereich [0 .. 2*Pi]
    angle[0] = atan2(S[0].y - M.y, S[0].x - M.x);
    if (angle[0] < 0.0) { angle[0]+= 2.0*M_PI; }

    angle[1] = atan2(S[1].y - M.y, S[1].x - M.x);
    if (angle[1] < 0.0) { angle[1]+= 2.0*M_PI; }

    // ****************************************************
    // Liegen die Schnittpunkte auf der Linie ? ***********

    if ((n[0] >= -0.01) && (n[0] <= 1.01)) { isInLine[0] = true; } else { isInLine[0] = false; }
    if ((n[1] >= -0.01) && (n[1] <= 1.01)) { isInLine[1] = true; } else { isInLine[1] = false; }

    // ****************************************************
    // Liegen die Schnittpunkte im gültigen Bogenabschnitt ?

    if (arcAngle[0] == arcAngle[1])
    {
      isInArc[0] = isInArc[1] = true;
    }
    else
    if (arcAngle[1] > arcAngle[0])
    {
      if ((angle[0] >= arcAngle[0]) && (angle[0] <= arcAngle[1])) { isInArc[0] = true; } else { isInArc[0] = false; }
      if ((angle[1] >= arcAngle[0]) && (angle[1] <= arcAngle[1])) { isInArc[1] = true; } else { isInArc[1] = false; }
    }
    else
    {
      if ((angle[0] <= arcAngle[1]) || (angle[0] >= arcAngle[0])) { isInArc[0] = true; } else { isInArc[0] = false; }
      if ((angle[1] <= arcAngle[1]) || (angle[1] >= arcAngle[0])) { isInArc[1] = true; } else { isInArc[1] = false; }
    }

    // ****************************************************

    if ((isInLine[0]) && (isInArc[0])) { isOk[0] = true; } else { isOk[0] = false; }
    if ((isInLine[1]) && (isInArc[1])) { isOk[1] = true; } else { isOk[1] = false; }
    
    // ****************************************************
    // Sowohl der 1. als auch der 2. Schnittpunkt liegen im
    // gültigen Kreisabschnitt. Den dem Punkt P an nächsten
    // gelegenen Schnittpunkt auswählen! ******************
    
    if ((isOk[0]) && (isOk[1]))
    {
      if ((S[0].x-P.x)*(S[0].x-P.x) + (S[0].y-P.y)*(S[0].y-P.y) < (S[1].x-P.x)*(S[1].x-P.x) + (S[1].y-P.y)*(S[1].y-P.y))
      {
        intersectionAngle = angle[0];
        intersectionPoint = S[0];
      }
      else
      {
        intersectionAngle = angle[1];
        intersectionPoint = S[1];
      }
    }
    else

    // Nur der 1. Schnittpunkt liegt im gültigen Kreisabschnitt
    if (isOk[0])
    {
      intersectionAngle = angle[0];
      intersectionPoint = S[0];
    }
    else

    // Nur der 2. Schnittpunkt liegt im gültigen Kreisabschnitt
    if (isOk[1])
    {
      intersectionAngle = angle[1];
      intersectionPoint = S[1];
    }
    else

    // Kein Schnittpunkt liegt im gültigen Kreisabschnitt
    {
      ok = false;
    }
    
  }

  return intersectionPoint;
}

// ********************************************************
// ********************************************************

double Circle::getTangentContactAngle(Point2D H, bool uhrzeiger)
{
  double p = H.x - M.x;
  double q = H.y - M.y;
  double radikant = p*p + q*q - r*r;
  double n;

  if (radikant < 0.0) { return 0.0; }
  n = sqrt(radikant);

  if (uhrzeiger) { n = -n; }

  return atan2(p*n+q*r, p*r-q*n);
}

// ********************************************************
// Gibt die Krümmung des Kreises zurück. ******************

double Circle::getBend()
{
  double a = r;
  if (a <= 0.0) { a = 0.0001; }
  return (-1.0/a);
}

    
// ********************************************************
// Gibt true zurück, wenn P innerhalb des Kreises liegt. **
    
bool Circle::isIncluding(Point2D P)
{
  Point2D w = P - M;
  return (w.x*w.x + w.y*w.y <= r*r);
}

// ********************************************************
// Gibt true zurück, wenn L komplett innerhalb des Kreises ist.

bool Circle::isIncluding(Line2D L)
{
  return (isIncluding(L.P[0]) && isIncluding(L.P[1]));
}

// ********************************************************
// Gibt true zurück, wenn es mind. einen Schnittpunkt mit *
// dem GESAMTEN Kreis gibt. *******************************

bool Circle::hasIntersectionWith(Line2D L)
{
  double p, q;
  double root;
  double t[2];
  
  Point2D w = L.P[0] - M;
  Point2D v = L.P[1] - L.P[0];
  v.normalize();
  double length = L.getLength();

  p = 2.0*(w.x*v.x + w.y*v.y);
  q = w.x*w.x + w.y*w.y - r*r;
  root = 0.25*p*p - q;

  if (root <= 0.0) { return false; }

  root = sqrt(root);
  t[0] = -0.5*p + root;
  t[1] = -0.5*p - root;

  if ((t[0] >= 0.0) && (t[0] <= length)) { return true; }
  if ((t[1] >= 0.0) && (t[1] <= length)) { return true; }

  return false;
}

// ********************************************************
// Gibt die gemeinsame Tangente dieses und des Kreises K **
// als Linie zwischen den beiden Berührungspunkten zurück.
// Es wird die Tangente links des Vektors [M -> K.M] zurückgegeben.

Line2D Circle::getCommonLeftTangentWith(Circle K)
{
  Line2D L(Point2D(0.0, 0.0), Point2D(1.0, 1.0));   // Rückgabewert

  double length = (K.M - M).magnitude();
  if (length == 0.0) { length = 0.0001; }
  double cosAngle = (r - K.r) / length;

  if ((cosAngle > -1.0) && (cosAngle < 1.0))
  {
    double angle = acos(cosAngle);    // Gibt immer den pos. Winkel zurück
    Point2D w = K.M - M;
    w.turn(angle);
    w.normalize();

    L.P[0] = M + w*r;
    L.P[1] = K.M + w*K.r;
  }
  
  return L;
}


// ****************************************************************************
// Gibt den Winkel des Punktes auf einem Kreis zur�ck, an dem eine Tangente 
// durch den au�erhalb liegenden Punkt H den Kreis ber�hrt.
// ****************************************************************************

double getCircleTangent(Point2D H, Point2D C, double r, bool uhrzeiger)
{
  double p = H.x - C.x;
  double q = H.y - C.y;
  double radikant = p*p + q*q - r*r;
  double n;

  if (radikant < 0.0) { return 0.0; }
  n = sqrt(radikant);

  if (uhrzeiger) { n = -n; }

  return atan2(p*n+q*r, p*r-q*n);
}

// ----------------------------------------------------------------------------
// Ellipse   
// ----------------------------------------------------------------------------

// ****************************************************************************
// Konstruktor ****************************************************************

Ellipse2D::Ellipse2D()
{
  set(Point2D(0.0, 0.0), 1.0, 1.0);
}
    
// ****************************************************************************
// Konstruktor ****************************************************************

Ellipse2D::Ellipse2D(Point2D _M, double _halfWidth, double _halfHeight)
{
  set(_M, _halfWidth, _halfHeight);
}

// ****************************************************************************
// Setzt die Parameter der Ellipse und berechnet einige zusätzliche Grössen. **

void Ellipse2D::set(Point2D _M, double _halfWidth, double _halfHeight)
{
  M          = _M;
  halfWidth  = _halfWidth;
  halfHeight = _halfHeight;

  long i;
  Point2D propP[NUM_ELLIPSE_PROPS];
  
  // Zu jedem Stützwinkel den entsprechenden Punkt auf der Ellipse berechnen
  
  for (i=0; i < NUM_ELLIPSE_PROPS; i++)
  {
    propAngle[i] = 2.0*M_PI*(double)i / (double)(NUM_ELLIPSE_PROPS-1);
    propP[i] = getPoint(propAngle[i]);
  }

  // Die Entfernung jedes Stützpunktes vom Anfang ermitteln
  propPosition[0] = 0.0;
  
  for (i=1; i < NUM_ELLIPSE_PROPS; i++)
  {
    propPosition[i] = propPosition[i-1] + (propP[i] - propP[i-1]).magnitude(); 
  }

  perimeter = propPosition[NUM_ELLIPSE_PROPS-1];
}

// ****************************************************************************
// Gibt einen Punkt auf der Ellipse zurück ************************************

Point2D Ellipse2D::getPoint(double angle)
{
  return Point2D(cos(angle)*halfWidth, sin(angle)*halfHeight);
}

// ****************************************************************************
// Gibt den Punkt an der Stelle t (in[0..1]) zurück, wenn t ein Maß für die ***
// Bogenlänge vom Punkt beim Winkel 0 aus gesehen ist, und t auf den Umfang ***
// der Ellipse normiert ist. **************************************************

double Ellipse2D::getAngle(double t)
{
  // Rausfinden, in welchem Stützabschnitt sich t befindet...
  double u = t*perimeter;
  if (u < 0.0) { u = 0.0; }
  if (u > perimeter) { u = perimeter; }
  
  long i = 0;

  while ((i < NUM_ELLIPSE_PROPS-2) && (u > propPosition[i+1])) { i++; }

  double a = propPosition[i+1] - propPosition[i];
  if (a == 0.0) { a = 0.0001; }

  return propAngle[i] + ((propAngle[i+1]-propAngle[i])*(u-propPosition[i])) / a;
}

// ****************************************************************************
// Gibt den ungefähren Umfang der Ellipse zurück ******************************

double Ellipse2D::getPerimeter()
{
  return perimeter;
}

// ****************************************************************************
// Returns the angle of a point P on an ellipse, that is given by its midpoint
// C and the half-axes a (horizontal) and b (vertical). The point is the one, 
// for wich HP is a tangent to the ellipse.
// ****************************************************************************

double getEllipseTangent(Point2D H, Point2D C, double a, double b, bool uhrzeiger)
{
  Point2D R = C - H;

  double radikant = R.x*R.x*b*b - b*b*a*a + R.y*R.y*a*a;
  if (radikant < 0.0) { radikant = 0.0; }
  double root = sqrt(radikant);
  
  double deno = R.x*R.x*b*b + R.y*R.y*a*a;
  const double EPSILON = 0.000001;
  if (fabs(deno) < EPSILON) { deno = EPSILON; }

  double alpha0 = atan2(
    -b*(R.y*a*a + R.x*root)/deno,
    -a*(R.x*b*b - R.y*root)/deno);

  double alpha1 = atan2(
    -b*(R.y*a*a - R.x*root)/deno,
    -a*(R.x*b*b + R.y*root)/deno);

  // Vector from H to the point on the ellipse at angle alpha0
  double COS = cos(alpha0);
  double SIN = sin(alpha0);
  Point2D v1 = C + Point2D(a*COS, b*SIN) - H;
  // Tangent to ellipse at angle alpha0
  Point2D v2(-a*SIN, b*COS);    

  double test = v1.x*v2.x + v1.y*v2.y;

  if (((test >= 0.0) && (uhrzeiger)) || ((test < 0.0) && (uhrzeiger == false)))
  {
    return alpha1;
  }
  else
  {
    return alpha0;
  }
}

// ****************************************************************************
