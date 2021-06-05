#include "BVP.h"

Real *BVP::GenerateF(int n)
{
  if (Dim==1)return GenerateF1(n);else return GenerateF2(n);
}

Real * BVP::GenerateF1(int n)
{
  Real h=1.0/n;
  --n;
  Real * temp=new Real [n];
  for (int i=0;i<=n-1;++i)
    {
      temp[i]=F((i+1)*h);
    }
  temp[0]+=B(0)/h/h;
  temp[n-1]+=B(1)/h/h;
  Blocks.push_back(temp);
  return temp;
}

Real BVP::GetBoundaryValue(Real x,Real y)
{
  if (x==0 || x==1 || y==0 || y==1)
    {
      return B(x,y);
    }else return 0;
}

Real * BVP::GenerateF2(int n)
{
  Real h=1.0/n;
  --n;
  Real * temp=new Real [n*n];
  for (int i=0;i<n;++i)
    for (int j=0;j<n;++j)
    {
      int index;
      GetIndex(index,i,j,n);
      Real x=(i+1)*h; Real y=(j+1)*h;
      temp[index]=F(x,y);
      Real ofx[4]={h,-h,0,0};
      Real ofy[4]={0,0,h,-h};
      for (int k=0;k<=3;++k)
	{
	  temp[index]+=GetBoundaryValue(x+ofx[k],y+ofy[k])/h/h;
	}
    }
  Blocks.push_back(temp);
  return temp;
}

int BVP::GetVectorSize(int n)
{
  if (Dim==1)return n-1;else return (n-1)*(n-1);
}

void BVP::reset()
{
  for (auto x:Blocks)
    {
      delete []x;
    }
  Blocks.clear();
}
