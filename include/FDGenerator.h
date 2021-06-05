#ifndef FDGENERATOR_H
#define FDGENERATOR_H
#include "Utility.h"
#include "BVP.h"

class FDGenerator{
 public:
  FDGenerator(){};
 FDGenerator(int Dim):Dim(Dim){}
  Real * GenerateA(int n,Real scale);
  Real * GenerateA1(int n,Real scale);
  Real * GenerateA2(int n,Real scale);
  void reset();
  int GetMatrixSize(int n);
  ~FDGenerator(){
    reset();
  };
 private:
  int Dim;
  std::vector<Real *> Blocks;
};


#else
//Do nothing!
#endif
