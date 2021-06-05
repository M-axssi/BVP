#ifndef UTILITY_H
#define UTILITY_H
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <vector>
#include <lapacke.h>
#include <cblas.h>

using Real=double;

void GetCoordinate(int i,int &x,int &y,int n);

Real GetValue(int x,int y,Real *F,int N);

Real GetMax(int n,Real *Y);

void GetIndex(int &i,int x,int y,int n);

int GetIndex(int x,int y,int n);

Real Quadratic(Real pos,Real y1,Real y2,Real y3);

Real * Eye(int n);

Real * KroneckerProduct(Real * A,Real *B,int n);

#else
//Do nothing!
#endif
