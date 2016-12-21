// simulación de la propagación de una onda
//ecuación de onda viscoelástica completa para un cambio de medio con un poro
#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
using namespace std;

class Matriz
{
private:
  double ***E,***UFm1,***UF0,***UF1,***V;
  int Lx;
  double Dx;
  int Ly; 
  double Dy;
  int Lz;
  int Dimx,Dimy,Dimz;
  double Dz;
  int ax; 
  int bx; 
  int ay; 
  int by; 
  int az; 
  int bz; 
  double Dt;
  int mux;
  int muy;
  int  muz;
  double W;
  double T;
  double V1;
  double S0;
  double Amplitud;
 public:
  Matriz(double,double,double,double,int,int,int,double);//constructor
  ~Matriz(void);//destructor
  void CargueParam(double Sourr[], double SS0);
  void Modelo(void);  
  void Velocidad(double VecVel[], int NCubos,int NDiv[]);
  void Inicie(double t);
  void Evolucion(double t);
  double RegistroEstaciones(double t, int Si, int Sj,int Sk);
  void AnimacionInicie(void);
  void Animacion(double t);
  void Cambio(void);  
  //void Imprima(void);
  
};

Matriz::Matriz(double Dimmmx,double Dimmmy, double Dimmmz, double Dttt, int LLx, int LLy, int LLz,double T)
{
  int i,j;
  
  Dimx=Dimmmx;
  Dimy=Dimmmy;
  Dimz=Dimmmz;
  Dt=Dttt;
  Lx=LLx;
  Ly=LLy;
  Lz=LLz;
  W=2*M_PI/T;
  
  UFm1=new double **[Dimx+1]; 
  UF0=new double **[Dimx+1];
  UF1=new double **[Dimx+1];
  V=new double **[Dimx+1];
  E=new double **[Dimx+1];
  for(i=0;i<=Dimx;i++)
    {
      UFm1[i]=new double *[Dimy+1];      
      UF0[i]=new double *[Dimy+1];
      UF1[i]=new double *[Dimy+1];
      V[i]=new double *[Dimy+1];
      E[i]=new double *[Dimy+1];
      for(j=0;j<=Dimy;j++)
	{
	  UFm1[i][j]=new double [Dimz+1];      
	  UF0[i][j]=new double [Dimz+1];
	  UF1[i][j]=new double [Dimz+1];
	  V[i][j]=new double [Dimz+1];
	  E[i][j]=new double [Dimz+1];
	}
    }
}

Matriz::~Matriz(void)
{
  int i,j;
  for(i=0;i<=Dimx;i++)
    {
      for(j=0;j<=Dimy;j++)
	{
	  delete [] UFm1[i][j];
	  delete [] UF0[i][j];
	  delete [] UF1[i][j];
	  delete [] V[i][j];
	  delete [] E[i][j];
	}
      delete [] UFm1[i];
      delete [] UF0[i];
      delete [] UF1[i];
      delete [] V[i];
      delete [] E[i];  
    }
      delete [] UFm1;
      delete [] UF0;
      delete [] UF1;
      delete [] V;
      delete [] E;  
}
void Matriz::CargueParam(double Sourr[],double SS0)
{
  Dx=1.0*Lx/Dimx; //distancia espacial en x. entre punto y punto
  Dy=1.0*Ly/Dimy; //distancia espacial en y. entre punto y punto
  Dz=1.0*Lz/Dimz; //distancia espacial en x. entre punto y punto
  mux=Sourr[0];//posicion de la fuente en km
  muy=Sourr[1];//posicion de la fuente en km
  muz=Sourr[2];//posicion de la fuente en km
  S0=SS0;//2*M_PI/(T*V1/Dx);//es mas o menos K
  Amplitud= S0; //limites en la amplitud para la animacion
  ax=0; //x inicial para la animacion
  bx=Lx; //limites en la amplitud para la animacion
  ay=0; //x inicial para la animacion
  by=Ly; //limites en la amplitud para la animacion
  az=0; //x inicial para la animacion
  bz=Lz; //limites en la amplitud para la animacion
}
void Matriz::Modelo(void)
{
  int i,j,k;
  for(i=0;i<=Dimx;i++)
    {
      for(j=0;j<=Dimy;j++)
	{
	  for(k=0;k<=Dimz;k++)
	    {
	      if(k>=0&&k<=Dimz/4)
		{
		  E[i][j][k]=1;//medio1
		}
	      else
		{
		  if(k>Dimz/4 && k<=Dimz/2)
		    {
		      E[i][j][k]=2;//medio2
		    }
		  else
		    {
		      if(k>Dimz/2 && k<=3*Dimz/4)
			{
			  E[i][j][k]=3;//medio3
			}
		      else
			{
			  E[i][j][k]=4;//medio4
			}
		    }
		}
	    }
	}
    }
}

void Matriz::Velocidad(double VecVel[],int NCubos, int NDiv[])
{ 
  int i,j,k,p=0;
  int ii,jj,kk;
  double Delx,Dely,Delz;
  double ValVel;
  
  Delx=Dimx/NDiv[0];
  Dely=Dimy/NDiv[1];
  Delz=Dimz/NDiv[2];

  //cout<<Delx<<" "<<Dely<<" "<<Delz<<endl;
  
  for(j=0;j<NDiv[1];j++)      
    {
      for(i=0;i<NDiv[0];i++)      
	{
	  for(k=0;k<NDiv[2];k++)	  
	    {
	      ValVel=VecVel[p];
	      for(jj=(j*Dely);jj<(j+1)*Dely;jj++)
		{
		  for(ii=(i*Delx);ii<(i+1)*Delx;ii++)  
		    {
		      for(kk=(k*Delz);kk<(k+1)*Delz;kk++)
			{
			  V[ii][jj][kk]=ValVel;
			} 
		    } 
		}
	      p++;
	    }
	}
    }
  
  for(j=0;j<=Dimy;j++)
    {
      for(k=0;k<=Dimz;k++)
	{
	  V[Dimx][j][k]=V[Dimx-1][j][k];
	}
    }
  for(i=0;i<=Dimx;i++)
    {
      for(k=0;k<=Dimz;k++)
	{
	  V[i][Dimy][k]=V[i][Dimy-1][k];
	}
    }
  for(i=0;i<=Dimx;i++)
    {
      for(j=0;j<=Dimy;j++)
	{
	  V[i][j][Dimz]=V[i][j][Dimz-1];
	}
    }
  /*
  for(j=12;j<=16;j++)
    {
      for(k=0;k<=Dimz;k++)
	{
	  for(i=0;i<=Dimx;i++)
	    {
	      cout<<V[i][j][k]<<" ";
	    }
	  cout<<endl;
	}
      cout<<endl<<endl;
    }
  */
  
}

void Matriz::Inicie(double t)
{ 
  int i,j,k;
  double lam;
  double r=0;
  
  /******************************Cálculo de UFm1********************************/
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)    
	{ 
	  for(k=1;k<=Dimz-1;k++)    
	    {
	      if((i==mux)&&(j==muy)&&(k==muz))
		{
		  UFm1[i][j][k]=S0*sin(W*t);
		}
	      else
		{
		  UFm1[i][j][k]=0;
		}
	    }
	} 
    }
  
  /************************Fronteras para UFm1********/
  // en direccion j
  for(i=1;i<=Dimx-1;i++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  UFm1[i][0][k]=UFm1[i][2][k];
	  UFm1[i][Dimy][k]=UFm1[i][Dimy-2][k];
	}
    }
  //en direccion i
  for(j=1;j<=Dimy-1;j++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  UFm1[0][j][k]=UFm1[2][j][k];
	  UFm1[Dimx][j][k]=UFm1[Dimx-2][j][k];
	}
    }
  //en direccion k
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)    
	{
	  UFm1[i][j][0]=UFm1[i][j][2];
	  UFm1[i][j][Dimz]=UFm1[i][j][Dimz-2];
	}
    }
  
  //bordes 
  for(j=1;j<=Dimy-1;j++)
    {
      UFm1[0][j][0]=0;
      UFm1[Dimx][j][0]=0;
      UFm1[0][j][Dimz]=0;
      UFm1[Dimx][j][Dimz]=0;
    }
  for(i=1;i<=Dimx-1;i++)
    {
      UFm1[i][0][0]=0;
      UFm1[i][Dimy][0]=0;
      UFm1[i][0][Dimz]=0;
      UFm1[i][Dimy][Dimz]=0;
    }
  for(k=1;k<=Dimz-1;k++)
    {
      UFm1[0][0][k]=0;
      UFm1[0][Dimy][k]=0;
      UFm1[Dimx][Dimy][k]=0;
      UFm1[Dimx][0][k]=0;
    }

  //esquinas superior z=0
  UFm1[0][0][0]=0;//
  UFm1[0][Dimy][0]=0;//
  UFm1[Dimx][Dimy][0]=0;//
  UFm1[Dimx][0][0]=0;//
  //esquinas inferiores z=Dimz
  UFm1[0][0][Dimz]=0;//
  UFm1[0][Dimy][Dimz]=0;//
  UFm1[Dimx][Dimy][Dimz]=0;//
  UFm1[Dimx][0][Dimz]=0;//
  

  /******************************Cálculo de UF0 a partir de UFm1****************/
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)
	{
	  for(k=1;k<=Dimz-1;k++)    
	    {
	      if((i==mux)&&(j==muy)&&(k==muz))
		{
		  UFm1[i][j][k]=S0*sin(W*(t+Dt));
		}
	      else
		{
	      lam=Dt*V[i][j][k];
	      UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
		}
	    }
	}
    }
  //tapas laterales en j
  j=0;
  for(i=1;i<=Dimx-1;i++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j+2][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
	}
    }
  
  j=Dimy;
  for(i=1;i<=Dimx-1;i++)
    { 
      for(k=1;k<=Dimz-1;k++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j-2][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
	}
    }
//tapas laterales en i
  i=0;
  for(j=1;j<=Dimy-1;j++)
    { 
      for(k=1;k<=Dimz-1;k++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i+2][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
	}
    }
  i=Dimx;
  for(j=1;j<=Dimy-1;j++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i-2][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));

	}
    }
  //tapas laterales en k
  k=0;
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k+2]))));
	}
    }
  k=Dimz;
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)
	{
	  lam=Dt*V[i][j][k];
	  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k-2]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
	}
    }
  
  ////////////Esquinas arriba
  i=0;j=0,k=0;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i+2][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j+2][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k+2]))));

  i=0;j=Dimy;k=0;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i+2][j][k]))+((1/(Dy*Dy))*(UFm1[i][j-2][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k+2]))));
  
  i=Dimx;j=Dimy;k=0;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i-2][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j-2][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k+2]))));
  
  i=Dimx;j=0;k=0;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i-2][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j+2][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k+2]))));
  ////////////Esquinas abajo
  i=0;j=0;k=Dimz;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i+2][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j+2][k]))+((1/(Dz*Dz))*(UFm1[i][j][k-2]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
  
  i=0;j=Dimy;k=Dimz;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i+2][j][k]))+((1/(Dy*Dy))*(UFm1[i][j-2][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k-2]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
  
  i=Dimx;j=Dimy;k=Dimz;  //esquina 
  lam=Dt*V[i][j][k];
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i-2][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j-2][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k-2]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
  
  i=Dimx;j=0;k=Dimz;  //esquina
  lam=Dt*V[i][j][k]; 
  UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i-2][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j+2][k]))+((1/(Dz*Dz))*(UFm1[i][j][k-2]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));

////bordes
  for(j=1;j<=Dimy-1;j++)
    {
      UF0[0][j][0]=UFm1[0][j][0]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[1][j][0]-(2*UFm1[0][j][0])+UFm1[2][j][0]))+((1/(Dy*Dy))*(UFm1[0][j+1][0]-(2*UFm1[0][j][0])+UFm1[0][j-1][0]))+((1/(Dz*Dz))*(UFm1[0][j][1]-(2*UFm1[0][j][0])+UFm1[0][j][2]))));
      
      UF0[Dimx][j][0]=UFm1[Dimx][j][0]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[Dimx-2][j][0]-(2*UFm1[Dimx][j][0])+UFm1[Dimx-1][j][0]))+((1/(Dy*Dy))*(UFm1[Dimx][j+1][0]-(2*UFm1[Dimx][j][0])+UFm1[Dimx][j-1][0]))+((1/(Dz*Dz))*(UFm1[Dimx][j][1]-(2*UFm1[Dimx][j][0])+UFm1[Dimx][j][2]))));
					
     UF0[0][j][Dimz]=UFm1[0][j][Dimz]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[1][j][Dimz]-(2*UFm1[0][j][Dimz])+UFm1[2][j][Dimz]))+((1/(Dy*Dy))*(UFm1[0][j+1][Dimz]-(2*UFm1[0][j][Dimz])+UFm1[0][j-1][Dimz]))+((1/(Dz*Dz))*(UFm1[0][j][Dimz-2]-(2*UFm1[0][j][Dimz])+UFm1[0][j][Dimz-1]))));
     
     UF0[Dimx][j][Dimz]=UFm1[Dimx][j][Dimz]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[Dimx-2][j][Dimz]-(2*UFm1[Dimx][j][Dimz])+UFm1[Dimx-1][j][Dimz]))+((1/(Dy*Dy))*(UFm1[Dimx][j+1][Dimz]-(2*UFm1[Dimx][j][Dimz])+UFm1[Dimx][j-1][Dimz]))+((1/(Dz*Dz))*(UFm1[Dimx][j][Dimz-2]-(2*UFm1[Dimx][j][Dimz])+UFm1[Dimx][j][Dimz-1]))));
    }

  for(i=1;i<=Dimx-1;i++)
    {
      UF0[i][0][0]=UFm1[i][0][0]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][0][0]-(2*UFm1[i][0][0])+UFm1[i-1][0][0]))+((1/(Dy*Dy))*(UFm1[i][1][0]-(2*UFm1[i][0][0])+UFm1[i][2][0]))+((1/(Dz*Dz))*(UFm1[i][0][1]-(2*UFm1[i][0][0])+UFm1[i][0][2]))));
     
      UF0[i][Dimy][0]=UFm1[i][Dimy][0]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][Dimy][0]-(2*UFm1[i][Dimy][0])+UFm1[i-1][Dimy][0]))+((1/(Dy*Dy))*(UFm1[i][Dimy-2][0]-(2*UFm1[i][Dimy][0])+UFm1[i][Dimy-1][0]))+((1/(Dz*Dz))*(UFm1[i][Dimy][1]-(2*UFm1[i][Dimy][0])+UFm1[i][Dimy][2]))));
      
      UF0[i][0][Dimz]=UFm1[i][0][Dimz]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][0][Dimz]-(2*UFm1[i][0][Dimz])+UFm1[i-1][0][Dimz]))+((1/(Dy*Dy))*(UFm1[i][1][Dimz]-(2*UFm1[i][0][Dimz])+UFm1[i][2][Dimz]))+((1/(Dz*Dz))*(UFm1[i][0][Dimz-2]-(2*UFm1[i][0][Dimz])+UFm1[i][0][Dimz-1]))));
 
      UF0[i][Dimy][Dimz]=UFm1[i][Dimy][Dimz]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][Dimy][Dimz]-(2*UFm1[i][Dimy][Dimz])+UFm1[i-1][Dimy][Dimz]))+((1/(Dy*Dy))*(UFm1[i][Dimy-2][Dimz]-(2*UFm1[i][Dimy][Dimz])+UFm1[i][Dimy-1][Dimz]))+((1/(Dz*Dz))*(UFm1[i][Dimy][Dimz-2]-(2*UFm1[i][Dimy][Dimz])+UFm1[i][Dimy][Dimz-1]))));
    }
    
  for(k=1;k<=Dimz-1;k++)
    {
      UF0[0][0][k]=UFm1[0][0][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[1][0][k]-(2*UFm1[0][0][k])+UFm1[2][0][k]))+((1/(Dy*Dy))*(UFm1[0][1][k]-(2*UFm1[0][0][k])+UFm1[0][2][k]))+((1/(Dz*Dz))*(UFm1[0][0][k+1]-(2*UFm1[0][0][k])+UFm1[0][0][k-1]))));

      UF0[0][Dimy][k]=UFm1[0][Dimy][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[1][Dimy][k]-(2*UFm1[0][Dimy][k])+UFm1[2][Dimy][k]))+((1/(Dy*Dy))*(UFm1[0][Dimy-2][k]-(2*UFm1[0][Dimy][k])+UFm1[0][Dimy-1][k]))+((1/(Dz*Dz))*(UFm1[0][Dimy][k+1]-(2*UFm1[0][Dimy][k])+UFm1[0][Dimy][k-1]))));

      UF0[Dimx][Dimy][k]=UFm1[Dimx][Dimy][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[Dimx-2][Dimy][k]-(2*UFm1[Dimx][Dimy][k])+UFm1[Dimx-1][Dimy][k]))+((1/(Dy*Dy))*(UFm1[Dimx][Dimy-2][k]-(2*UFm1[Dimx][Dimy][k])+UFm1[Dimx][Dimy-1][k]))+((1/(Dz*Dz))*(UFm1[Dimx][Dimy][k+1]-(2*UFm1[Dimx][Dimy][k])+UFm1[Dimx][Dimy][k-1]))));

      UF0[Dimx][0][k]=UFm1[Dimx][0][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[Dimx-2][0][k]-(2*UFm1[Dimx][0][k])+UFm1[Dimx-1][0][k]))+((1/(Dy*Dy))*(UFm1[Dimx][1][k]-(2*UFm1[Dimx][0][k])+UFm1[Dimx][2][k]))+((1/(Dz*Dz))*(UFm1[Dimx][0][k+1]-(2*UFm1[Dimx][0][k])+UFm1[Dimx][0][k-1]))));
    }
    //UF0[i][j][k]=UFm1[i][j][k]+((lam*lam/2.0)*(((1/(Dx*Dx))*(UFm1[i+1][j][k]-(2*UFm1[i][j][k])+UFm1[i-1][j][k]))+((1/(Dy*Dy))*(UFm1[i][j+1][k]-(2*UFm1[i][j][k])+UFm1[i][j-1][k]))+((1/(Dz*Dz))*(UFm1[i][j][k+1]-(2*UFm1[i][j][k])+UFm1[i][j][k-1]))));
}
  
void Matriz::Evolucion(double t)
{ 
  int i,j,k;
  double parte1,parte2,parte3;
  
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)    
	{
	  for(k=1;k<=Dimz-1;k++)
	    {
	      if((i==mux)&&(j==muy)&&(k==muz))
		{
		  UFm1[i][j][k]=S0*sin(W*t);
		}
	      else
		{
		  parte1=(1.0/(Dx*Dx))*(UF0[i+1][j][k]-(2*UF0[i][j][k])+UF0[i-1][j][k]);
		  parte2=(1.0/(Dy*Dy))*(UF0[i][j+1][k]-(2*UF0[i][j][k])+UF0[i][j-1][k]);
		  parte3=(1.0/(Dz*Dz))*(UF0[i][j][k+1]-(2*UF0[i][j][k])+UF0[i][j][k-1]);
		  UF1[i][j][k]=((V[i][j][k]*V[i][j][k]*Dt*Dt)*(parte1 + parte2+ parte3))+(2.0*UF0[i][j][k])-UFm1[i][j][k];
		}
	    }
	}
    }
  
  /************************Fronteras para UF1********/
  // en direccion j
  for(i=1;i<=Dimx-1;i++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  UF1[i][0][k]=UF1[i][2][k];
	  UF1[i][Dimy][k]=UF1[i][Dimy-2][k];
	}
    }
  //en direccion i
  for(j=1;j<=Dimy-1;j++)
    {
      for(k=1;k<=Dimz-1;k++)
	{
	  UF1[0][j][k]=UF1[2][j][k];
	  UF1[Dimx][j][k]=UF1[Dimx-2][j][k];
	}
    }
  //en direccion k
  for(i=1;i<=Dimx-1;i++)
    {
      for(j=1;j<=Dimy-1;j++)    
	{
	  UF1[i][j][0]=UF1[i][j][2];
	  UF1[i][j][Dimz]=UF1[i][j][Dimz-2];
	}
    }
  
  //bordes 
  for(j=1;j<=Dimy-1;j++)
    {
      UF1[0][j][0]=UF1[2][j][2];
      UF1[Dimx][j][0]=UF1[Dimx-2][j][2];
      UF1[0][j][Dimz]=UF1[2][j][Dimz-2];
      UF1[Dimx][j][Dimz]=UF1[Dimx-2][j][Dimz-2];
    }
  for(i=1;i<=Dimx-1;i++)
    {
      UF1[i][0][0]=UF1[i][2][2];
      UF1[i][Dimy][0]=UF1[i][Dimy-2][2];
      UF1[i][0][Dimz]=UF1[i][2][Dimz-2];
      UF1[i][Dimy][Dimz]=UF1[i][Dimy-2][Dimz-2];
    }
  for(k=1;k<=Dimz-1;k++)
    {
      UF1[0][0][k]=UF1[2][2][k];
      UF1[0][Dimy][k]=UF1[2][Dimy-2][k];
      UF1[Dimx][Dimy][k]=UF1[Dimx-2][Dimy-2][k];
      UF1[Dimx][0][k]=UF1[Dimx-2][2][k];
    }

  //esquinas superior z=0
  UF1[0][0][0]=UF1[2][2][2];//
  UF1[0][Dimy][0]=UF1[2][Dimy-2][2];//
  UF1[Dimx][Dimy][0]=UF1[Dimx-2][Dimy-2][2];//
  UF1[Dimx][0][0]=UF1[Dimx-2][2][2];//
  //esquinas inferiores z=Dimz
  UF1[0][0][Dimz]=UF1[2][2][Dimz-2];//
  UF1[0][Dimy][Dimz]=UF1[2][Dimy-2][Dimz-2];//
  UF1[Dimx][Dimy][Dimz]=UF1[Dimx-2][Dimy-2][Dimz-2];//
  UF1[Dimx][0][Dimz]=UF1[Dimx-2][2][Dimz-2];//

}
 
void Matriz::AnimacionInicie(void)
{
  cout<<"set xrange["<<ax<<":"<<bx<<"]"<<endl
      <<"set yrange["<<bz<<":"<<az<<"]"<<endl
      <<"set zrange["<<-1*Amplitud<<":"<<Amplitud<<"]"<<endl
      <<"set xlabel 'Distancia x'"<<endl
      <<"set ylabel 'Profundidad'"<<endl
      <<"set zlabel 'amplitud'"<<endl
      <<"set size square "<<endl
      <<"set ytics auto"<<endl
    // <<"set ytics ('bottom' 100, '' 1, 'top' 0)"<<endl
    //  <<"set origin 0.0, 0.0"<<endl
    //  <<"set pm3d "<<endl
      <<"set pm3d map"<<endl
      <<"set palette defined(-0.05'dark-red',-0.025 'red',0 'black', 0.025 'yellow',0.05 'dark-orange')"<<endl;
  // <<"set palette rgbformulae 22,13,-31"<<endl;
  //<<"set palette gray"<<endl;
}
 
void Matriz::Animacion(double t)
{
  int i=0,j=0,k;
  
  cout<<"set title 'tiempo="<<t<<" '"<<endl;
  cout <<"set cbrange["<<-1*Amplitud<<":"<<Amplitud<<"]"<<endl; 
  cout<<"splot '-' w pm3d"<<endl;
  
  for(k=0;k<=Dimz;k++)
    {
      j=Dimy/2;
      //for(j=0;j<=Dimy;j++)
      //{
	  for(i=0;i<=Dimx;i++)
	    cout<<i*Dx<<"    "<<k*Dz<<"    "<<UF0[i][j][k]<<endl;
	  cout<<endl; 
	  //}
    }
  cout<<"e"<<endl;  
}

void Matriz::Cambio(void)
{
  int i,j,k;
  for(i=0;i<=Dimx;i++)
  {
    for(j=0;j<=Dimy;j++)
      {
	for(k=0;k<=Dimz;k++)
	  {
	    UFm1[i][j][k]=UF0[i][j][k];
	    UF0[i][j][k]=UF1[i][j][k];
	  }
      }
  }
}

double Matriz::RegistroEstaciones(double t, int Si, int Sj,int Sk)
{

  double val;
  val=UFm1[Si][Sj][Sk];
	
  return val;
}

     
