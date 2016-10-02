#include <iostream>
#include "MathVector.h"
#include "MathMatrix.h"

// test logistic regression
// use gradient & newton

// main logistic regression function
void GetLogitCoeff( Vector & betaHat, Matrix & X, Vector & Y );

// sigmoid function
double sigmoid( double x);

// convert Y to +1 and -1. Record the original value.
void convertBinaryTrait(Vector & Y, int & Yplus1, int & Yminus1);

// norm between 2 vectors
double NormVector( Matrix & v1, Matrix & v2);

// determinant of a square matrix
double GetDet( int n, Matrix & m);

// check if Y contains only 2 values
bool isBinaryY( Vector & Y);

