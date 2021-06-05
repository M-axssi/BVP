#include "FDGenerator.h"

Real* FDGenerator::GenerateA(int n,Real scale)
{
  if (Dim==1)
    {
      return GenerateA1(n,scale);
    }else
    {
      return GenerateA2(n,scale);
    }
}

Real* FDGenerator::GenerateA1(int n,Real scale)
{
  --n;
  Real * A=new Real [n*n]();
  Blocks.push_back(A);
  for (int i=0;i<=n-1;++i)
    {
      int index;
      GetIndex(index,i,i,n);
      A[index]=2;
    }
  for (int i=0;i<=n-2;++i)
    {
      int index1,index2;
      GetIndex(index1,i+1,i,n);
      GetIndex(index2,i,i+1,n);
      A[index1]=A[index2]=-1;
    }
  cblas_dscal(n*n,scale,A,1);
  return A;
}

Real* FDGenerator::GenerateA2(int n,Real scale)
{
  Real *tempA=GenerateA1(n,scale);
  -- n;
  Real *E=Eye(n);
  Real *A1;
  Real *A2;
  A1=KroneckerProduct(tempA,E,n);
  A2=KroneckerProduct(E,tempA,n);
  cblas_daxpy(n*n*n*n,1,A2,1,A1,1);
  Blocks.push_back(A1);
  delete [] A2;
  delete [] E;
  return A1;
}

int FDGenerator::GetMatrixSize(int n)
{
  if (Dim==1)return n-1;else return (n-1)*(n-1);
}

void FDGenerator::reset()
{
  for (auto x:Blocks)
    {
      delete [] x;
    }
  Blocks.clear();
}
