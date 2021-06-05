#ifndef MULTIGRID_H
#define MULTIGRID_H
#include "Utility.h"
#include "FDGenerator.h"
#include "BVP.h"

class MultiGrid{
 public:
  MultiGrid(){};
 MultiGrid(bool RO,bool IO,bool CY,int MI,int _v1,int _v2,Real _epison):IsFullWeighting(RO),
    IsLinear(IO),IsVcycle(CY),MaxInterations(MI),v1(_v1),v2(_v2),epison(_epison){};

  void ProlongationV(Real *,Real *,int ,int,int);
  void ProlongationVDim1(Real *,Real *,int ,int);
  void ProlongationVDim2(Real *,Real *,int ,int);
  void RestrictionV(Real * ,Real * ,int,int,int);
  void RestrictionVDim1(Real *,Real *,int);
  void RestrictionVDim2(Real *,Real *,int ,int );
  void Solve(FDGenerator & G,BVP & B,Real *U,int n,Real * RU,std::ofstream & doc);
  void Vcycle_Solve(FDGenerator &G, BVP & B,Real * U,int n,Real * RU,std::ofstream& doc);
  void Vcycle(Real *A,Real *F,Real *U,int n,int layer,FDGenerator & G,BVP &B,int N);
  void WeightedJacobi(Real *A,Real *F,Real *U,int n,int count);
  void FullMultigrid_Solve(FDGenerator &G, BVP & B,Real * U,int n,Real * RU,std::ofstream &doc);
  Real GetResidule(Real * A,Real *F,Real *U,int n,Real * Y);
  Real GetError(Real *U,Real *RU,int n);
  void FullMultigrid(Real *A,Real *F,Real *U,int n,int layer,FDGenerator & G,BVP & B,int N);
 private:
  bool IsFullWeighting=1;
  bool IsLinear=1;
  bool IsVcycle=1;
  int MaxInterations=10;
  int v1=1;
  int v2=1;
  Real epison=1e-8;
};

#else
//Do nothing!
#endif
