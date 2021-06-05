#ifndef FUNCTION_H
#define FUNCTION_H
#include "Utility.h"

class Function
{
 public:
  Function(){};
  virtual Real operator()(Real x) =0;
  virtual Real operator()(Real x,Real y)=0;
};

class Dim1BoundaryF:public Function
{
public:
 Dim1BoundaryF():Function(){u0=u1=0;};
 Dim1BoundaryF(Real u0,Real u1):Function(),u0(u0),u1(u1){};
  virtual Real operator()(Real x);
  virtual Real operator()(Real x,Real y);
private:
  Real u0,u1;
};

#else
//Do nothing!
#endif
