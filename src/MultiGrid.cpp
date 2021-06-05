#include "MultiGrid.h"

void MultiGrid::Solve(FDGenerator & G,BVP & B,Real *U,int n,Real * RU,std::ofstream & doc)
{
  if (IsVcycle) Vcycle_Solve(G,B,U,n,RU,doc);
  else FullMultigrid_Solve(G,B,U,n,RU,doc);
  G.reset();
  B.reset();
}

Real MultiGrid::GetResidule(Real * A,Real *F,Real *U,int n,Real * Y)
{
  cblas_dgemv(CblasColMajor, CblasNoTrans, n, n, 1, A, n, U,1, 0, Y, 1);
  cblas_daxpy(n,-1,F,1,Y,1);
  return GetMax(n,Y);
}

Real MultiGrid::GetError(Real *U,Real *RU,int n)
{
  Real Y[n];
  cblas_dcopy(n,U,1,Y,1);
  cblas_daxpy(n,-1,RU,1,Y,1);
  return GetMax(n,Y);
}

void MultiGrid::WeightedJacobi(Real *A,Real *F,Real *U,int n,int count)
{
  const Real w=2.0/3;
  Real *invD=new Real[n*n]();
  Real *LAU=new Real[n*n]();
  Real *T=new Real[n*n];
  Real c[n];
  for (int i=0;i<=n*n-1;++i)
    {
      int x,y;
      GetCoordinate(i,x,y,n);
      if (x==y)invD[i]=1/A[i];
      else LAU[i]=A[i];
    }
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,n,n,n,-1,invD,n,LAU,n,0,T,n);
  cblas_dgemv (CblasColMajor,CblasNoTrans,n,n,1,invD,n,F,1,0,c,1);
  for (int i=1;i<=count;++i)
    {
      Real tempu[n];
      cblas_dgemv (CblasColMajor,CblasNoTrans,n,n,1,T,n,U,1,0,tempu,1);
      cblas_daxpy(n,1,c,1,tempu,1);
      cblas_dscal(n,w,tempu,1);
      cblas_daxpy(n,1-w,U,1,tempu,1);
      cblas_dcopy(n,tempu,1,U,1);
    }
  delete invD;
  delete LAU;
  delete T;
}

void MultiGrid::RestrictionV(Real *F2,Real *F,int Dim,int N1,int N2)
{
  if (Dim==1) RestrictionVDim1(F2,F,N1);else
    RestrictionVDim2(F2,F,N1,N2);
}

void MultiGrid::RestrictionVDim2(Real *F2,Real *F,int N1,int N2)
{
  for (int i=0;i<N1;++i)
    for (int j=0;j<N1;++j)
      {
	int index;
	GetIndex(index,i,j,N1);
	int x=2*i+1;
	int y=2*j+1;
	int Cindex;
	GetIndex(Cindex,x,y,N2);
	if (!IsFullWeighting)
	  {
	    F2[index]=F[Cindex];
	    continue;
	  }else
	  F2[index]=F[Cindex]/4;
        int lox[4]={1,-1,0,0};
	int loy[4]={0,0,1,-1};
	for (int i=0;i<4;++i)
	  {
	    int tIndex;
	    GetIndex(tIndex,x+lox[i],y+loy[i],N2);
	    F2[index]+=F[tIndex]/8;
	  }
	int dox[4]={1,1,-1,-1};
	int doy[4]={1,-1,1,-1};
	for (int i=0;i<4;++i)
	  {
	    int tIndex;
	    GetIndex(tIndex,x+dox[i],y+doy[i],N2);
	    F2[index]+=F[tIndex]/16;
	  }
      }
}

void MultiGrid::RestrictionVDim1(Real *F2,Real *F,int N)
{
  for (int i=0;i<N;++i)
    {
      if (IsFullWeighting)
	{
	  F2[i]=1.0/4*(2*F[2*i+1]+F[2*i]+F[2*i+2]);
	}else
	{
	  F2[i]=F[i*2+1];
	}
    }
}

void MultiGrid::ProlongationV(Real *F,Real *F2,int Dim,int N1,int N2)
{
  if (Dim==1) ProlongationVDim1(F,F2,N1,N2);else
    ProlongationVDim2(F,F2,N1,N2);
}

void MultiGrid::ProlongationVDim2(Real *F,Real *F2,int N1,int N2)
{
  for (int i=0;i<N1;++i)
    for (int j=0;j<N1;++j)
      {
	int index;
	GetIndex(index,i,j,N1);
	if (i%2==0 && j%2==0)
	  {
	    Real f0,f1,f2,f3;
	    f0=GetValue(i/2,j/2,F2,N2);
	    f1=GetValue(i/2-1,j/2,F2,N2);
	    f2=GetValue(i/2,j/2-1,F2,N2);
	    f3=GetValue(i/2-1,j/2-1,F2,N2);
	    F[index]=(f0+f1+f2+f3)/4.0;
	  }else
	  if (i%2==0)
	    {
	      Real f0,f1;
	      f0=GetValue(i/2-1,(j-1)/2,F2,N2);
	      f1=GetValue(i/2,(j-1)/2,F2,N2);
	      F[index]=(f0+f1)/2.0;
	    }else
	    if (j%2==0)
	    {
	      Real f0,f1;
	      f0=GetValue((i-1)/2,j/2-1,F2,N2);
	      f1=GetValue((i-1)/2,j/2,F2,N2);
	      F[index]=(f0+f1)/2.0;
	    }else
	      {
		F[index]=GetValue((i-1)/2,(j-1)/2,F2,N2);
	      }
      }
}

void MultiGrid::ProlongationVDim1(Real * F,Real *F2,int N,int N2)
{
  if (IsLinear)
    for (int i=0;i<N;++i)
      {
	if (i==0) F[0]=(F2[0])/2;
	else if (i==N-1) F[i]=(F2[i/2-1])/2;
	else 
	  {
	    if (i%2==1) F[i]=F2[(i-1)/2];
	    else F[i]=(F2[(i-1)/2]+F2[i/2])/2;
	  }
      }else
    {
      for (int i=0;i<N;++i)
	{
	  if (i==0) F[i]=Quadratic(0.5,0,F2[0],F2[1]);
	  else if (i==N-1) F[i]=Quadratic(1.5,F2[N2-2],F2[N2-1],0);
	  else if (i==N-3) F[i]=Quadratic(0.5,F2[N2-2],F2[N2-1],0);
	  else if (i%2==1) F[i]=F2[(i-1)/2];
	  else F[i]=Quadratic(0.5,F2[(i-1)/2],F2[i/2],F2[i/2+1]);
	}
    }
}

void MultiGrid::FullMultigrid(Real *A,Real *F,Real *U,int n,int layer,FDGenerator & G,BVP &B,int N)
{
  if (layer==4 || n<=2)
    {
      Vcycle(A,F,U,n,layer+1,G,B,N);
    }else
    {
      int Dim=B.GetDim();
      int N2=B.GetVectorSize(n/2);
      Real h=B.GetH(n/2);
      Real F2h[N2];
      Real V2h[N2];
      Real *A2h=G.GenerateA(n/2,1.0/h/h);
      RestrictionV(F2h,F,Dim,n/2-1,n-1);
      RestrictionV(V2h,U,Dim,n/2-1,n-1);
      FullMultigrid(A2h,F2h,V2h,n/2,layer+1,G,B,N2);
      ProlongationV(U,V2h,Dim,n-1,n/2-1);
      Vcycle(A,F,U,n,layer,G,B,N);
    }
}

void MultiGrid::Vcycle(Real *A,Real *F,Real *U,int n,int layer,FDGenerator & G,BVP& B,int N)
{
  if (layer==5 || n<=2)
    {
      int ipiv[N];
      Real A_[N*N];
      cblas_dcopy(N*N,A,1,A_,1);
      cblas_dcopy(N,F,1,U,1);
      LAPACKE_dgesv(LAPACK_COL_MAJOR, N, 1,A_, N, ipiv, U, N);
    }else
    {
      int Dim=B.GetDim();
      WeightedJacobi(A,F,U,N,v1);
      int N2=B.GetVectorSize(n/2);
      Real h=B.GetH(n/2);
      Real tempF[N];
      Real F2h[N2];
      Real V2h[N2]={0};
      Real *A2h=G.GenerateA(n/2,1.0/h/h); 
      Real error=GetResidule(A,F,U,N,tempF);
      RestrictionV(F2h,tempF,Dim,n/2-1,n-1);      
      Vcycle(A2h,F2h,V2h,n/2,layer+1,G,B,N2);
      Real tempV[N];
      ProlongationV(tempV,V2h,Dim,n-1,n/2-1);
      cblas_daxpy(N,-1,tempV,1,U,1);
      WeightedJacobi(A,F,U,N,v2);
    }
}

void MultiGrid::Vcycle_Solve(FDGenerator &G, BVP & B,Real * U,int n,Real * RU,std::ofstream &doc)
{
  Real h=B.GetH(n);
  int N=G.GetMatrixSize(n);
  Real RUnorm=GetMax(N,RU);
  Real * A=G.GenerateA(n,1.0/h/h);
  Real * F=B.GenerateF(n);
  int count=1;
  Real Y[N];
  Real presidule=GetResidule(A,F,U,N,Y);
  Real perror=GetError(U,RU,N);
  doc<<"V-cycle : "<<0<<", Residule = "<<presidule<<", Error = "<<perror<<std::endl;
  while (count<=MaxInterations && (perror/RUnorm)>epison)
    {
      Vcycle(A,F,U,n,1,G,B,N);
      Real residule=GetResidule(A,F,U,N,Y);
      Real error=GetError(U,RU,N);
      doc<<"V-cycle : "<<count<<", Residule = "<<residule<<" Ratio = "<<residule/presidule<<" , Error = "<<error<<" Ratio = "<<error/perror<<std::endl;
      presidule=residule; perror=error;
      ++count;
      G.reset();
      B.reset();
      A=G.GenerateA(n,1.0/h/h);
      F=B.GenerateF(n);
    }
}

void MultiGrid::FullMultigrid_Solve(FDGenerator &G, BVP & B,Real * U,int n,Real * RU,std::ofstream &doc)
{
  Real h=B.GetH(n);
  int N=G.GetMatrixSize(n);
  Real RUnorm=GetMax(N,RU);
  Real * A=G.GenerateA(n,1.0/h/h);
  Real * F=B.GenerateF(n);
  int count=1;
  Real Y[N];
  Real presidule=GetResidule(A,F,U,N,Y);
  Real perror=GetError(U,RU,N);
  doc<<"FMG : "<<0<<", Residule = "<<presidule<<", Error = "<<perror<<std::endl;
  Real pdivn=perror/RUnorm;
  while (count<=MaxInterations && pdivn>epison)
    {
      Real tempU[N]={0};
      FullMultigrid(A,Y,tempU,n,1,G,B,N);
      cblas_daxpy(N,-1,tempU,1,U,1);
      Real residule=GetResidule(A,F,U,N,Y);
      Real error=GetError(U,RU,N);
      doc<<"FMG : "<<count<<", Residule = "<<residule<<" Ratio = "<<residule/presidule<<" , Error = "<<error<<" Ratio = "<<error/perror<<std::endl;
      presidule=residule; perror=error;
      ++count;
      G.reset();
      B.reset();
      A=G.GenerateA(n,1.0/h/h);
      F=B.GenerateF(n);
    }
}
