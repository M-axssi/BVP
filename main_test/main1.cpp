#include <iostream>
#include "MultiGrid.h"
#include "Function.h"

class LapalaceF:public Function
{
public:
  LapalaceF():Function(){};
  virtual Real operator()(Real x);
  virtual Real operator()(Real x,Real y);
};

Real LapalaceF::operator()(Real x)
{
  return std::exp(std::sin(x))*(std::sin(x)-std::cos(x)*std::cos(x));
}

Real LapalaceF::operator()(Real x,Real y)
{
  return 0;
}

int main()
{
  std::ifstream fin("Dim1_Input.txt");
  Real u0,u1;
  u0=u1=0;
  int n;
  bool IsHomogeneous;
  bool IsFullWeighting;
  bool IsLinear;
  bool IsVcycle;
  bool HaveInitialGuess;
  int MaxInterations;
  Real epison;
  fin>>n>>IsHomogeneous;
  if (!IsHomogeneous) {u0=std::exp(std::sin(0)); u1=std::exp(std::sin(1));}
  Dim1BoundaryF Boundary(u0,u1);
  LapalaceF F;
  BVP B(1,Boundary,F);
  FDGenerator G(1);
  Real U[n-1]={0};
  Real RU[n-1];
  fin>>IsFullWeighting>>IsLinear>>IsVcycle>>MaxInterations>>epison>>HaveInitialGuess;
  if (HaveInitialGuess)
    {
      std::ifstream finit("Dim1_InitialGuess.txt");
      for (int i=0;i<n-1;++i)
	{
	  Real x;
	  finit>>x;
	  U[i]=x;
	}
    }
  Real h=1.0/n;
  Real t1,t2;
  t1=std::exp(std::sin(0)); t2=std::exp(std::sin(1));
  for (int i=0;i<n-1;++i)
    if (!IsHomogeneous)
      RU[i]=std::exp(std::sin((i+1)*h));
    else RU[i]=std::exp(std::sin((i+1)*h))+(t1-t2)*(i+1)*h-t1;

  std::ofstream doc("AnswerDim1.out");
  doc<<"n = "<<n<<std::endl;
  doc<<"Boundary conditions : ";
  if (IsHomogeneous) doc<<"homogeneous, ";else doc <<"nonhomogeneous, ";
  doc<<"Restriction operators : ";
  if (IsFullWeighting) doc<<"full-weighting, ";else doc<<"injection, ";
  doc <<"Interpolation operators : ";
  if (IsLinear) doc<<"linear, ";else doc<<"quadratic,";
  doc<<"Cycles : ";
  if (IsVcycle) doc<<"V-cycle"<<std::endl;else doc<<"FMG"<<std::endl;
  doc<<"Stopping criteria : "<<"epison = "<<epison<<" max-iterations = "<<MaxInterations<<std::endl;
  if (HaveInitialGuess) doc<<"Have initial guess"<<std::endl;else doc<<"Have no initial guess"<<std::endl;
  
  MultiGrid  MSolver(IsFullWeighting,IsLinear,IsVcycle,MaxInterations,2,1,epison);
  MSolver.Solve(G,B,U,n,RU,doc);
  
  
  std::ofstream fout("Dim1_Output.m");
  fout<<std::fixed;
  fout.precision(10);
  fout<<"X=["<<0<<",";
  for (int i=0;i<n-1;++i)
    fout<<(i+1)*h<<",";
  fout<<1<<"]"<<std::endl;
  fout<<"U=["<<u0<<",";
  for (int i=0;i<n-1;++i)
    {
      fout<<U[i]<<",";
    }
  fout<<u1<<"]"<<std::endl;
  fout<<"plot(X,U)";
}
