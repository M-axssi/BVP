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
  std::ofstream doc("../result/Dim1NoHomFinalRE.txt");
  bool HaveInitialGuess=false;
  bool IsHomogeneous=false;
  int MaxInterations=100;
  Real epison=2.2e-16;
  
  Real u0,u1;
 u0=std::exp(std::sin(0));
 u1=std::exp(std::sin(1));
  int _n[4]={128,256,512,1024};
  bool _IL[2]={0,1};
  bool _IV[2]={0,1};
  bool _IF[2]={0,1};
  for (int i1=0;i1<=3;++i1)
    for (int i2=0;i2<=1;++i2)
      for (int i3=0;i3<=1;++i3)
	for (int i4=0;i4<=1;++i4)
	{
  
	  int n=_n[i1];
	  bool IsFullWeighting=_IF[i2];
	  bool IsLinear=_IL[i3];
	  bool IsVcycle=_IV[i4];

	  Dim1BoundaryF Boundary(u0,u1);
	  LapalaceF F;
	  BVP B(1,Boundary,F);
	  FDGenerator G(1);
	  Real U[n-1]={0};
	  Real RU[n-1];

	  Real h=1.0/n;
	  
	  Real t1,t2;
	  t1=std::exp(std::sin(0)); t2=std::exp(std::sin(1));
	  for (int i=0;i<n-1;++i)
	    if (!IsHomogeneous)
	      RU[i]=std::exp(std::sin((i+1)*h));
	    else RU[i]=std::exp(std::sin((i+1)*h))+(t1-t2)*(i+1)*h-t1;

  
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
	  doc<<std::endl<<std::endl;
	}  
 
}
