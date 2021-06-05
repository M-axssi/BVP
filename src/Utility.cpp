#include "Utility.h"

void GetCoordinate(int i,int &x,int &y,int n)
{
  x=i%n;
  y=i/n;
}

void GetIndex(int &i,int x,int y,int n)
{
  i=y*n+x;
}

Real GetMax(int n,Real *Y)
{
  Real ma=-1;
  for (int i=0;i<n;++i)
    if (std::fabs(Y[i]))
      ma=fabs(Y[i]);
  return ma;
}

int GetIndex(int x,int y,int n)
{
  return y*n+x;
}

Real GetValue(int x,int y,Real *F,int N)
{
  if (x>=N  || y>=N) return 0;
  if (x<0 || y<0) return 0;else return F[GetIndex(x,y,N)];
}

Real Quadratic(Real pos,Real y1,Real y2,Real y3)
{
  return (y1-2*y2+y3)/2*pos*pos+(-3.0*y1/2+2*y2-y3/2)*pos+y1;
}

Real * Eye(int n)
{
  Real * E=new Real [n*n]();
  for (int i=0;i<n;++i)
    {
      int index;
      GetIndex(index,i,i,n);
      E[index]=1; 
    }
  return E;
}

Real * KroneckerProduct(Real * A,Real *B,int n)
{
  Real *temp=new Real[n*n*n*n];
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j)
      {
	int index;
	GetIndex(index,i,j,n);
	Real C=A[index];
	Real x,y;
	x=i*n;
	y=j*n;
	for (int k=0;k<n;++k)
	  for (int l=0;l<n;++l)
	    {
	      int index1,index2;
	      GetIndex(index1,x+k,y+l,n*n);
	      GetIndex(index2,k,l,n);
	      temp[index1]=C*B[index2];
	    }
      }
  return temp;
}
