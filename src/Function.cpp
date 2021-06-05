#include "Function.h"

Real Dim1BoundaryF::operator()(Real x)
{
  if (x==0) return u0;
  if (x==1) return u1;
  return 0;
}

Real Dim1BoundaryF::operator()(Real x,Real y)
{
  return 0;
}
