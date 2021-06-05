#ifndef BVP_H
#define BVP_H
#include "Utility.h"
#include "Function.h"

class BVP{
 public:
  //  BVP(){a=0;b=1;u0=un=0;}
 BVP(int _Dim,Function &_B,Function &_F):Dim(_Dim),B(_B),F(_F){};
  Real * GenerateF(int n);
  Real * GenerateF1(int n);
  Real * GenerateF2(int n);
  void reset();
  int GetDim(){return Dim;};
  int GetVectorSize(int n);
  Real GetH(int n){return 1.0/n;};
  Real GetBoundaryValue(Real x,Real y);
  ~BVP()
    {
      reset();
    }
 private:
  std::vector<Real *> Blocks;
  Function &F;
  Function &B;
  int Dim=1;
};

#endif
