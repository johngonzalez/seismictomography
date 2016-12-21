//Clase para solucionar Ecuaciones Lineales por factorizacion LU

#include <iostream>
#include <cmath>
#include <stdio.h>

using namespace std;

class Sistema{
  
  int N;
  double * a;
  double * anew ;  
  
  double * b;
  double * bnew;
  double * x;
public:
  void Inicie (double*,double*,int);
  void LU();  
  double *Muestre();
  void Termine(void);
};

void Sistema::Inicie(double * a0,double * b0, int N0) {
  N=N0;
  anew=new double [N];
  b=   new double [N];
  
  bnew=new double [N];
  x=   new double [N];
  a=   new double [N];
  
  anew=a0;
  b=b0;
}

void Sistema::LU(){
  
  int i,j,k,p,s,w;
  double suma,l[N*N],u[N*N];
  j=0; s=1; p=1;
  for(i=0;i<N;i++){
    if(s<0){i--;}
    for(j=0;j<N;j++){
      if(i>j){
	suma=0;
	for(k=0;k<j;k++)
	  suma+=l[i+k*N]*u[k+j*N];
	l[i+j*N]=(anew[i+j*N]-suma)/u[j+j*N];
      }
      else if(i==j)
	l[i+j*N]=1;
      else
	l[i+j*N]=0;
    }
    j=0;
    s=1;
    while(j<N && s>0){
      if(i>j)
	{u[i+j*N]=0; j++;} 
      else{
	suma=0;
	for(k=0;k<i;k++)
	  suma+=l[i+k*N]*u[k+j*N];
	u[i+j*N]=anew[i+j*N]-suma;
	if(i==j && u[i+j*N]==0){
	  for(w=0;w<N;w++)
	    {anew[i+p+w*N]=a[i+w*N]; anew[i+w*N]=a[i+p+w*N];} //intercambie filas
	  for(w=0;w<N;w++)
	    {a[i+p+w*N]=anew[i+p+w*N]; a[i+w*N]=anew[i+w*N];}
	  
	  bnew[i+p]=b[i];   bnew[i]=b[i+p];
	  b[i+p]=bnew[i+p]; b[i]=bnew[i];
	  
	  p++;
	  s=-1;
	}
	else{j++; p=1;}
      }
    }
  }
  
  double y[N];
  
  for(i=0;i<N;i++){
    suma=0;
    for(j=0;j<i;j++)
      suma+=l[i+j*N]*y[j];
    y[i]=b[i]-suma;
  }
  
  for(i=N-1;i>-1;i--){
    suma=0;
    for(j=N-1;j>i;j--){
      suma+=u[i+j*N]*x[j];
    }
    x[i]=(y[i]-suma)/u[i+i*N];
  }

  //Imprima Resultados
  //for(i=0;i<Nnew ;i++)
  //cout<<x[i]<<endl;
  
}

double *Sistema::Muestre(){
  int i,j;
  
  //for(i=0;i<N;i++){
  //for(j=0;j<N;j++){
    // cout<<anew[i+j*N]<<"   ";
    //}
    //cout<<endl;
  //}
  
  return x;
  //for(i=0;i<N;i++)
  //cout<<x[i]<<endl;
  

}
void Sistema::Termine(void){

  delete [] anew;
  delete [] a;
  delete [] b;
  delete [] bnew;
  delete [] x;
}
