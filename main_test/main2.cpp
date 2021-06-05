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
  std::ifstream fin("Dim2_Input.txt");
  int n;
  bool IsHomogeneous;
  bool IsFullWeighting;
  bool IsLinear;
  bool IsVcycle;
  bool HaveInitialGuess;
  int MaxInterations;
  Real epison;
  fin>>n>>IsHomogeneous;
  Dim2BoundaryF Boundary(IsHomogeneous);
  LapalaceF F;
  BVP B(2,Boundary,F);
  FDGenerator G(2);
  Real U[(n-1)*(n-1)]={0};
  fin>>IsFullWeighting>>IsLinear>>IsVcycle>>MaxInterations>>epison>>HaveInitialGuess;
  if (HaveInitialGuess)
    {
      std::ifstream finit("Dim2_InitialGuess.txt");
      for (int i=0;i<(n-1)*(n-1);++i)
	{
	  Real x;
	  finit>>x;
	  U[i]=x;
	}
    }
  Real h=1.0/n;
  Real RU[(n-1)*(n-1)];
  for (int i=0;i<(n-1)*(n-1);++i)
    {
      int x,y;
      GetCoordinate(i,x,y,n-1);
      RU[i]=std::exp(std::sin((x+1)*h)*std::sin((y+1)*h));
    }
  
  std::ofstream doc("AnswerDim2.out");
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
  
  
  std::ofstream fout("Dim2_Output.m");
  fout<<std::fixed;
  fout.precision(10);
  fout<<"[X,Y]=meshgrid(0:"<<h<<":1);"<<std::endl;
  fout<<"U=[";
  for (int i=0;i<=n;++i)
    {
      for (int j=0;j<=n;++j)
	{
	  Real x=i*h;
	  Real y=j*h;
	  Real V;
	  if (i==0 || j==0 || i==n || j==n) V=Boundary(x,y);else V=U[GetIndex(i-1,j-1,n-1)];
	  fout<<V<<" ";
	}
      fout<<";";
    }
  fout<<"]"<<std::endl;
  fout<<"mesh(X,Y,U)";
}
