#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <sstream>

//#include "SL.h"

using namespace std;

double * CambieNumD(double *M, int n,int m){

  double * CM;
  int i,k;
  
  CM=new double [m*n];
  
  for(i=0;i<m*n;i++){
    k=int(i/n);
    CM[i]=M[m*(i-n*k)+k];
  }
  return CM;
  delete [] CM;

}
int * CambieNumI(int *M, int n,int m){

  int * CM;
  int i,k;
  
  CM=new int [m*n];
  
  for(i=0;i<m*n;i++){
    k=int(i/n);
    CM[i]=M[m*(i-n*k)+k];
  }
  return CM;
  delete [] CM;

}
string *CambieNumS(string *M, int n,int m){

  string * CM;
  int i,k;
  
  CM=new string [m*n];
  
  for(i=0;i<m*n;i++){
    k=int(i/n);
    CM[i]=M[m*(i-n*k)+k];
  }
  return CM;
  delete [] CM;

}
double *ConviertaDouble(string * S,int n){
  double *Doub;
  stringstream ins[n];
  Doub=new double [n];
  for(int i=0;i<n;i++){ins[i].str(S[i]); ins[i]>>Doub[i];}
  return Doub;
  delete [] Doub;
}
int *ConviertaInt(string * S,int n){
  int *Int;
  stringstream ins[n];
  Int=new int [n];
  for(int i=0;i<n;i++){ins[i].str(S[i]); ins[i]>>Int[i];}
  return Int;
  delete [] Int;
}

void Muestre(double *M,int n,int m){
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++)
      cout<<M[j*n+i]<<"   ";
    cout<<endl;
  }
  cout<<endl;
}
void MuestreString(string *M,int n,int m){
  
  for(int i=0;i<n;i++){
    for(int j=0;j<m;j++)
      cout<<M[j*n+i]<<"   ";
    cout<<endl;
  }
  cout<<endl;
}

string *AdquieraDatos(const char *Nom,int n,int m,int LineasInicio){
  string *M,palabra;
  int i,j;
  ifstream FILE(Nom);
  char cadena[256];
  
  for(i=0;i<LineasInicio;i++)
    FILE.getline(cadena,256);
  
  M=new string [m*n+1];
  i=0;
  while (!FILE.eof()){
    FILE>>palabra;
    if(palabra=="SISMO") FILE.getline(cadena,256);
    else{
      M[0+i*m]=palabra;
      for(j=1;j<m;j++) FILE>>M[j+i*m];
      i++;
    }
  }
  FILE.close();
  
  M=CambieNumS(M,n,m);
  return M;
  delete [] M;
}
string *AdquieraDatos2(const char *Nom,int n,int m,int LineasInicio){
  string *M;
  int i,j;
  ifstream FILE(Nom);
  char cadena[256];
  
  for(i=0;i<LineasInicio;i++)
    FILE.getline(cadena,256);
  
  M=new string [m*n];
  i=0;
  
  for(i=0;i<n;i++)
    for(j=0;j<m;j++) 
      FILE>>M[j+i*m];
  
  FILE.close();
  
  M=CambieNumS(M,n,m);
  return M;
  delete [] M;
}

double * Transpuesta(double *M,int n,int m){
  //TM[m][n]
  double * TM;
  int i,k;
  
  TM=new double [m*n];
  
  for(i=0;i<m*n;i++){
    k=int(i/m);
    TM[i]=M[n*(i-m*k)+k];
  }
  return TM;
  delete [] TM;
}

double * Multiplique(double *M1,int n1,int m1,double *M2,int n2,int m2){

  //cout<<n1<<" "<<m1<<" "<<n2<<" "<<m2<<endl;
  //cout<<m2<<endl;
  
  if(m1!=n2)
    cout<<"ERROR LAS MATRICES TIENEN MAL LA DIMENSION"<<endl;
  else
    {
      int i,j,k;
      double suma;
      double * M1M2;
      M1M2=new double [n1*m2];
      // cout<<endl;
      //Muestre(M2,n2,m2);
      
      
      for(i=0;i<n1;i++){
	for(j=0;j<m2;j++){
	  suma=0;
	  for(k=0;k<m1;k++)
	    suma=suma+M1[i+k*n1]*M2[k+j*n2];
	  M1M2[i+j*n1]=suma;
	}
      }
      return M1M2;
	delete [] M1M2;
	
    }
}
double * Identidad(int n){
  
  double * I;
  I=new double [n*n];
  
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      if(i==j)
	I[i+j*n]= 1;
      else
	I[i+j*n]= 0;
  return I;
  delete [] I;
}
double * MulEscalar(double a, double * M,int n, int m){
  
  double * aM;
  aM=new double [n*m];
  
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      aM[i+j*n]= a*M[i+j*n];
  return aM;
  delete [] aM;
}
double * Sume(double *M1,double *M2,int n, int m){
  
  double * M1masM2;
  M1masM2=new double [n*m];
  
  for(int i=0;i<n;i++)
    for(int j=0;j<n;j++)
      M1masM2[i+j*n]= M1[i+j*n]+M2[i+j*n];
  return M1masM2;
  delete [] M1masM2;
}

double **CambieFormato(double *M,int n,int m){
  int i,j;
  double **a;
  
  a=   new double *[n];
  for(j=0;j<m;j++)
    a[j]=new double [n];  

  for(i=0;i<n;i++)
    for(j=0;j<m;j++)
      a[i][j]=M[i+j*n];
  return a;
  
  for(i=0;i<n;i++) 
    delete [] a[j];  
  delete [] a;
}
double norma2(double *V,int n){
  double suma=0;
  for(int i=0;i<n;i++)
    suma+=V[i]*V[i];
  return suma;
}

double norma(double *V,int n){
  double suma=0;
  for(int i=0;i<n;i++)
    suma+=V[i]*V[i];
  return sqrt(suma);
}

void VariablesControl(const char *ArchPar,double &e,double &lim, double * &L){
  int Num;
  ifstream FILE(ArchPar);
  FILE>>e; FILE>>lim;
  for(int i=0;i<3*2;i++)
    FILE>>L[i];
  L=CambieNumD(L,3,2);
  FILE.close();
}

double *Posicion(const char *Nom,int n,string Estacion){
  int numero,s,i,j;
  double *Pos;
  Pos=new double[3];
  string *dato,Est;
  stringstream ins[3];
  dato=AdquieraDatos(Nom,n,4,1);
  
  s=1; i=0;
  while(s==1&&i<n){
    Est=dato[i];
    if(Est==Estacion){s=0;numero=i;}
    i++;
  }
 
  for(j=0;j<3;j++){
    ins[j].str(dato[numero+n*(j+1)]);
    ins[j]>>Pos[j];
  }
  return Pos;
  delete [] Pos;
}
int Numero(const char *Arch,int NumL){
  int Num;
  ifstream FILE(Arch);
  for(int i=0;i<NumL;i++)
    FILE>>Num;
  return Num;
  FILE.close();
}
int *NumeroBloques(const char *Arch){
  int *Num;
  Num=new int[3];
  ifstream FILE(Arch);
  FILE>>Num[0]>>Num[1]>>Num[2];
  return Num;
  FILE.close();
}



















int *NumeroEstacionesSismoFase(const char *ArchPic,int NumSis,string fase){
  
  int *i,k;
  ifstream FILE(ArchPic); char cadena[256];
  string palabra1,palabra2;
  
  i=new int [NumSis];
  k=-1;
  while (!FILE.eof()){
    FILE>>palabra1>>palabra2;
    if(palabra1=="SISMO") {k++; if(k!=NumSis-1) i[k]=0; else if (fase=="P") i[k]=0; else i[k]=-1;}
    else
      if(palabra2==fase)  
	i[k]++;
    FILE.getline(cadena,256);
  }
  return i;
  delete [] i;
}


string *NombreEstacionesSismoFase(const char *ArchPic,int NumPicT,string fase){

  int k;
  ifstream FILE(ArchPic); char cadena[256];
  string palabra1,palabra2;
  string *Est;
  
  Est=new string [NumPicT+1];
  k=0;
  while (!FILE.eof()){
    FILE>>palabra1>>palabra2;
    if(palabra1!="SISMO"&& palabra2==fase) {Est[k]=palabra1; k++;} 
    FILE.getline(cadena,256);
  }
  return Est; 
}

int *NumeroFecha(const char *ArchPic,int NumSis){

  int k;
  ifstream FILE(ArchPic); char cadena[256];
  string palabra;
  int *Fecha;
  Fecha =new int [3*NumSis];
  
  k=0;
  while (!FILE.eof()){
    FILE>>palabra;
    if(palabra=="SISMO") {FILE>>Fecha[k+0*NumSis]>>Fecha[k+1*NumSis]>>Fecha[k+2*NumSis]; k++;} 
    FILE.getline(cadena,256);
  }
  return Fecha;
}

int num(int ano,int mes){
  int n;
  if(mes==0)  n=0 ;
  if(mes==1)  n=31;
  if(mes==2)  if(ano%4==0) n=29; else n=28;
  if(mes==3)  n=31;
  if(mes==4)  n=30;
  if(mes==5)  n=31;
  if(mes==6)  n=30;
  if(mes==7)  n=31;
  if(mes==8)  n=31;
  if(mes==9)  n=30;
  if(mes==10) n=31;
  if(mes==11) n=30;
  if(mes==12) n=31;
  return n;
}
int Sdias(int ano,int mes){
  int suma;
  suma=0;
  for(int i=0; i<mes; i++)
    suma+=num(ano,i);
  return suma;
}
void Observacion(const char *ArchPic,int NumSis,int NumPicTP,int NumPicTSP,int NumEstCP, int NumEstCSP,int *NumPicP,int *NumPicSP,string *NomEstCP,string *NomEstCSP,string *NomEstP,string *NomEstSP,double *&TObs,string * &Fase){
  
  string *a,*aP,*aSP,*TObsS,*TObsSP,*TObsSSP,*FaseP,*FaseSP;
  int NumFil,ip,KP,KSP,k,i,j,s,jP,jSP,NumFilP,NumFilSP,NumPicT,NumEstC;
  
  NumPicT  =NumPicTP+NumPicTSP;
  NumEstC  =NumEstCP+NumEstCSP;
  NumFilP  =NumEstCP*NumSis;
  NumFilSP =NumEstCSP*NumSis;
  NumFil   =NumFilP+NumFilSP;

  TObsS    =new string [NumFil];
  TObsSP   =new string [NumFilP];
  TObsSSP  =new string [NumFilSP];
  TObs     =new double [NumFil];
    
  a        =new string [NumPicT*3];
  aP       =new string [NumPicTP*3];
  aSP      =new string [NumPicTSP*3];
  FaseP    =new string [NumFilP];
  FaseSP   =new string [NumFilSP];
  //Fase     =new string [NumFil];
 
  a=AdquieraDatos(ArchPic,NumPicT,3,0);
  
  jP=0; jSP=0;
  for(i=0;i<NumPicT;i++){
    if(a[i+1*NumPicT]=="P") { for(j=0;j<3;j++) aP[jP+j*NumPicTP]   =a[i+j*NumPicT]; jP++;}
    if(a[i+1*NumPicT]=="SP"){ for(j=0;j<3;j++) aSP[jSP+j*NumPicTSP]=a[i+j*NumPicT]; jSP++;}
  }
    
  KP=0; KSP=0;
  for(k=0;k<NumSis;k++){
    
    for(i=0;i<NumEstCP;i++){ 
      ip=0; s=1;
      while (s==1){
	if(NomEstCP[i]==NomEstP[ip+KP])    {TObsSP[i+NumEstCP*k]=aP[ip+KP+2*NumPicTP];          FaseP[i+NumEstCP*k]=aP[ip+KP+1*NumPicTP]; s=0;    }
	else	                           ip++;
	if(ip==NumPicP[k] && s==1)         {TObsSP[i+NumEstCP*k]="0";                           FaseP[i+NumEstCP*k]="N"; s=0;                     }
      }
    }
    KP+=NumPicP[k];
    for(i=0;i<NumEstCSP;i++){
      ip=0; s=1;
      while (s==1){
	if(NomEstCSP[i]==NomEstSP[ip+KSP]) {TObsSSP[i+NumEstCSP*k]=aSP[ip+KSP+2*NumPicTSP];     FaseSP[i+NumEstCSP*k]=aSP[ip+KSP+1*NumPicTSP]; s=0;}
	else	                           ip++;
	if(ip==NumPicSP[k] && s==1)        {TObsSSP[i+NumEstCSP*k]="0";                         FaseSP[i+NumEstCSP*k]="N"; s=0;                    }
      }
    }
    KSP+=NumPicSP[k];
    
    for(i=0;i<NumEstCP+NumEstCSP;i++)
      if(i<NumEstCP)                       {TObsS[i+NumEstC*k]=TObsSP[i+NumEstCP*k];            Fase[i+NumEstC*k]=FaseP[i+NumEstCP*k];             }
      else                                 {TObsS[i+NumEstC*k]=TObsSSP[i-NumEstCP+NumEstCSP*k]; Fase[i+NumEstC*k]=FaseSP[i-NumEstCP+NumEstCSP*k];  }
  }
  
  TObs=ConviertaDouble(TObsS,NumFil);
  
  delete [] TObsS;
  delete [] TObsSP;
  delete [] TObsSSP;
  
    
  delete [] a;
  delete [] aP;
  delete [] aSP;
  //delete [] Fase;
  delete [] FaseP;
  delete [] FaseSP;
}

double *PosicionEstaciones(const char *ArchEst,string *NomEstC,int NumEst,int NumEstC){
  //int NumPic,NumSis;
  double *Pos1,*Pos2;
  //string *a;
  //ifstream FILE(ArchPic);
 
  Pos1=new double[3];
  Pos2=new double[NumEstC*3];
  //a=new string[NumEstC*3];
  
  // a=AdquieraDatos(ArchPic,NumPicT,3,0);
  
  for(int i=0;i<NumEstC;i++){
    Pos1=Posicion(ArchEst,NumEst,NomEstC[i]);
    for(int j=0;j<3;j++) Pos2[i+NumEstC*j]=Pos1[j];
  }
  
  return Pos2;
  delete [] Pos1;
  delete [] Pos2;
  //delete [] a;
}
double *PosicionEstacionesTodas(const char *ArchEst){
  int NumEst;
  double *Pos;
  string *a,*b;
  ifstream FILE(ArchEst);
 
  FILE>>NumEst;
  
  Pos=new double[NumEst*3];
  a=new string[NumEst*4];
  b=new string[NumEst*3];
  
  a=AdquieraDatos(ArchEst,NumEst,4,1);
  
  for(int i=0;i<NumEst;i++){
    for(int j=0;j<3;j++) 
      b[i+NumEst*j]=a[i+NumEst*(j+1)];
  }
  Pos=ConviertaDouble(b,3*NumEst);
  
  return Pos;
  delete [] Pos;
  delete [] a;
  delete [] b;
}
string *NombreEstacionesTodas(const char *ArchEst){
  int NumEst;
  string *NomEst;
  string *a;
  ifstream FILE(ArchEst);
 
  FILE>>NumEst;
  
  NomEst=new string[NumEst];
  a=new string[NumEst*4];
  
  
  a=AdquieraDatos(ArchEst,NumEst,4,1);
  
  for(int i=0;i<NumEst;i++)
    NomEst[i]=a[i];
  
  return NomEst;
  delete [] NomEst;
  delete [] a;
}
string * NombreEstaciones(const char *ArchPic){
  int NumPic,NumSis;
  string *a,*NomEst;
  ifstream FILE(ArchPic);
  FILE>>NumSis>>NumPic;
  a     =new string[NumPic*3];
  NomEst=new string [NumPic];
  
  a=AdquieraDatos(ArchPic,NumPic,3,1);
  for(int i=0;i<NumPic;i++) NomEst[i]=a[i];
  
  return NomEst;  
  delete [] NomEst;
  delete [] a;
}
double *ModeloInicial(const char *Arch,int NumSis,int NumBlo){
  double *Mod;
  string *a,*ModS;
  ifstream FILE(Arch);
  ModS=new string [4*NumSis+NumBlo];
  Mod =new double [4*NumSis+NumBlo];
  a   =new string [4*NumSis+NumBlo];
  
  a=AdquieraDatos(Arch,4*NumSis+NumBlo,1,1);
  Mod=ConviertaDouble(a,4*NumSis+NumBlo);
  
  return Mod;
  delete [] Mod;
  delete [] ModS;
  delete [] a;
  FILE.close();
}
double *ModeloDentroLimites(int NumSis,double *Mod,double *L){
  
  for(int k=0;k<NumSis;k++)
    for(int j=0;j<3;j++){
      if(Mod[4*k+j]<L[j])   Mod[4*k+j]=L[j];
      if(Mod[4*k+j]>L[j+3]) Mod[j]=L[j+3];
    }
  return Mod;
  
}
void NumeroSismosPicadas(const char *ArchPic,int & NumSis, int & NumPicT){
  ifstream FILE(ArchPic);
  string palabra;
  char cadena[256];
  NumSis=0; NumPicT=-1;
  while (!FILE.eof()){
    FILE>>palabra;
    if(palabra=="SISMO")  NumSis++;
    else                  NumPicT++;
    FILE.getline(cadena,256);
  }
}

void NumeroPicadasFase(const char *ArchPic,int NumPicT,int & NumPicTP, int & NumPicTSP){
  ifstream FILE(ArchPic);
  string *a;
  char cadena[256];
  a=new string [NumPicT*3];
  a=AdquieraDatos(ArchPic,NumPicT,3,0);
  
  NumPicTP=0; NumPicTSP=0;
  for(int i=0; i<NumPicT;i++){
    if(a[i+1*NumPicT]=="P")  NumPicTP++;
    if(a[i+1*NumPicT]=="SP") NumPicTSP++;
  }
}
int NumeroEstacionesConsideradasFase(const char *ArchPic,int NumPicT,string fase){
  string *a,*b;
  int cont,i,s,j,NumEstCF;
  
  ifstream FILE(ArchPic);
  a= new string [NumPicT];
  b= new string [NumPicT*3];
  b=AdquieraDatos(ArchPic,NumPicT,3,0);
  
  cont=0; NumEstCF=0;
  
  for(j=0;j<NumPicT;j++){
    if(b[j+1*NumPicT]==fase){
      if(cont==0) {a[0]=b[j];NumEstCF=1;}
      else{
	s=1;i=0;
	while(s==1 && i<NumEstCF){
	  if(b[j]!=a[i])    {s=1;}
	  else              {s=0;}
	  if(i==NumEstCF-1 && s==1) {a[NumEstCF]=b[j]; s=0; NumEstCF++;}
	  i++;
	}}
      cont++;
    }
  }
  
  return NumEstCF;
  delete [] a;
  delete [] b;
}
int NumeroEstacionesConsideradasC(const char *ArchPic,int NumPicT){
  string *a,*b;
  int cont,i,s,j,NumEstCF;
  
  ifstream FILE(ArchPic);
  a= new string [NumPicT];
  b= new string [NumPicT*3];
  b=AdquieraDatos(ArchPic,NumPicT,3,0);
  
  cont=0; NumEstCF=0;
  
  for(j=0;j<NumPicT;j++){
    //if(b[j+1*NumPicT]==fase){
    if(cont==0) {a[0]=b[j];NumEstCF=1;}
    else{
      s=1;i=0;
      while(s==1 && i<NumEstCF){
	if(b[j]!=a[i])    {s=1;}
	else              {s=0;}
	if(i==NumEstCF-1 && s==1) {a[NumEstCF]=b[j]; s=0; NumEstCF++;}
	i++;
      }}
    cont++;
  }
  
  return NumEstCF;
  delete [] a;
  delete [] b;
}

string *EstacionesConsideradasFase(const char *ArchPic,int N,string fase){

  ifstream FILE(ArchPic);
  string palabra1,palabra2,*a;
  char cadena[256];
  int cont,i,j,s;
  
  a= new string [N];
  
  cont=0;
  while (!FILE.eof()){
    FILE>>palabra1>>palabra2;
    
    if(palabra1!="SISMO" && palabra2==fase){     
      if(cont==0) {a[0]=palabra1;j=1;}
      else{
	s=1;i=0;
	while(s==1 && i<j){
	  if(palabra1!=a[i]) {s=1;}
	  else               {s=0;}
	  if(i==j-1 && s==1) {a[j]=palabra1; s=0; j++;}
	  i++;
	}
      }
      cont++;
    }
    
    FILE.getline(cadena,256);
  }
  
  return a;
  delete [] a; 

}
string *EstacionesConsideradasC(const char *ArchPic,int N){
  
  ifstream FILE(ArchPic);
  string palabra1,palabra2,*a;
  char cadena[256];
  int cont,i,j,s;
  
  a= new string [N];
  
  cont=0;
  while (!FILE.eof()){
    FILE>>palabra1>>palabra2;
    
    if(palabra1!="SISMO"){     
      if(cont==0) {a[0]=palabra1;j=1;}
      else{
	s=1;i=0;
	while(s==1 && i<j){
	  if(palabra1!=a[i]) {s=1;}
	  else               {s=0;}
	  if(i==j-1 && s==1) {a[j]=palabra1; s=0; j++;}
	  i++;
	}
      }
      cont++;
    }
    
    FILE.getline(cadena,256);
  }
  
  return a;
  delete [] a; 
}
string *EstacionesConsideradas(int NumEstCP, int NumEstCSP, string *NomEstCP, string *NomEstCSP){
  string *NomEstC;
  NomEstC= new string [NumEstCP+NumEstCSP];

  for(int i=0;i<NumEstCP+NumEstCSP;i++)
    if(i<NumEstCP)
      NomEstC[i]=NomEstCP [i];
    else
      NomEstC[i]=NomEstCSP[i-NumEstCP];
  
  return NomEstC;
  delete [] NomEstC;
  
}
