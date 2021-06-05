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
  return 0; 
}

Real LapalaceF::operator()(Real x,Real y)
{
  Real sx=std::sin(x);
  Real sy=std::sin(y);
  Real exy=std::exp(sx*sy);
  return (exy*(2*sx*sy-sx*sx-sy*sy+2*sx*sx*sy*sy));
}

class Dim2BoundaryF:public Function
{
public:
  Dim2BoundaryF():Function(){};
  Dim2BoundaryF(bool _IH):Function(),IsHomogeneous(_IH){};
  virtual Real operator()(Real x);
  virtual Real operator()(Real x,Real y);
private:
  bool IsHomogeneous =0;
};

Real Dim2BoundaryF::operator()(Real x)
{
  return 0;
}

Real Dim2BoundaryF::operator()(Real x,Real y)
{
  if (IsHomogeneous) return 0;
  return(std::exp(std::sin(x)*std::sin(y)));
}

int main()
{
  std::ofstream doc("../result/Dim2NoHomResult.txt");
  bool HaveInitialGuess=false;
  bool IsHomogeneous=false;
  int MaxInterations=100;
  Real epison=1e-8;

  int _n[4]={16,32,48,64};
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
	    
	    Dim2BoundaryF Boundary(IsHomogeneous);
	    LapalaceF F;
	    BVP B(2,Boundary,F);
	    FDGenerator G(2);
	    Real U[(n-1)*(n-1)]={0};
	    
	    Real h=1.0/n;
	    Real RU[(n-1)*(n-1)];
	    for (int i=0;i<(n-1)*(n-1);++i)
	      {
		int x,y;
		GetCoordinate(i,x,y,n-1);
		RU[i]=std::exp(std::sin((x+1)*h)*std::sin((y+1)*h));
	      }
  
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
