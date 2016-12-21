#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
using namespace std;
double *Archivos(int NStts);

class Files
{
private:
  int NDate;
  double Dtime;
  int NStat;
  double **M;
  double tp;
  int C;
public:
  void Inicie(int NDatos,double Deltat,int NStts);
  void CarVal(double t,double A, int j);
  void CrearMatriz(void);
  void BorrarMatriz(void);
  void Calculo(int i);
  void Avance(double t, double A);
  double ValC(void);
  double Tiempo(void);
  void Mostrar(int j);
};

void Files::CrearMatriz(void)
{
  int i;
  M=new double *[3];
  for(i=0;i<3;i++)
    {
      M[i]=new double [2];
    }
}

void Files::BorrarMatriz(void)
{
  int i;
  
  for(i=0;i<3;i++)
    {
     delete [] M[i];
    }
  delete [] M;
}

double Files::ValC(void)
{
  return C;
}

void Files::Mostrar(int j)
{
  cout<<M[j][0]<<" "<<M[j][1]<<endl;
}

double Files::Tiempo(void)
{
  return tp;
}

void Files::Avance(double t, double A)
{
  M[0][0]=M[1][0];
  M[1][0]=M[2][0];
  M[2][0]=t;

  M[0][1]=M[1][1];
  M[1][1]=M[2][1];
  M[2][1]=A;
}

void Files::Inicie(int NDatos,double Deltat, int NStts)
{
  NDate=NDatos;
  Dtime=Deltat;
  NStat=NStts;
  C=0;
}

void Files::CarVal(double t,double A, int j)
{
  M[j][0]=t;
  M[j][1]=A;
}

void Files::Calculo(int i)
{
  if(C==0)
    {
     
      if((M[1][1]>M[0][1])&&(M[1][1]>M[2][1]))
	{
	  C=1;
	  tp=M[1][0];
	}
      else
	{
	  if((M[1][1]<M[0][1])&&(M[1][1]<M[2][1]))
	    {
	      C=1;
	      tp=M[1][0];
	    }
	  else
	    {
	      
	    }
	}
    }
}

double *Archivos(int NStts)
{
  int i,j,k;
  int NDatos;
  double Deltat,t,A;
  double VC=0;
  double *T;
  T=new double [NStts];
  ifstream InfoStation,Arch;
  Files File[NStts];
  Arch.open("RStations.dat");
  InfoStation.open("InfoStation.dat");
  InfoStation>>NDatos;
  InfoStation>>Deltat;
  
  for(i=0;i<NStts;i++)
    {
      File[i].Inicie(NDatos,Deltat,NStts);
      File[i].CrearMatriz();
    }
  //carga primeras lineas con k para las estaciones
  for(k=0;k<=1;k++)
    {
      for(j=0;j<=2;j++)
	{
	  Arch>>t;
	  for(i=0;i<NStts;i++)
	    {
	      Arch>>A;
	      File[i].CarVal(t,A,j);
	    }
	}
    }
  
  for(j=1;j<=NDatos-3;j++)
    {
      for(i=0;i<NStts;i++)
	{
	  File[i].Calculo(i);
	}
      
      Arch>>t;
      for(i=0;i<NStts;i++)
	{	  
	  Arch>>A;
	  File[i].Avance(t,A);
	}
      
      for(i=0;i<NStts;i++)
	{
	  VC=VC+File[i].ValC();
	}
      if(VC==NStts)
	{
	  for(i=0;i<NStts;i++)
	    {
	      //cout<<File[i].Tiempo()<<endl;
	      T[i]=File[i].Tiempo();
	    }
	  j=NDatos;//para salir del ciclo
	}
      else
	{
	  VC=0;
	}
    }
  
  for(i=0;i<NStts;i++)
    {
      File[i].BorrarMatriz();
    }

  InfoStation.close();
  Arch.close();
  /*
  for(i=0;i<NStts;i++)
    {
      cout<<T[i]<<endl;
    }
  */
  return T;
  delete []T;
}
