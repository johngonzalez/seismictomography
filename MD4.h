#include <iostream>
#include <cmath>
#include <stdio.h>
#include "Estructura3.h"
#include "Max.h"


double * CalculeTF(int NEST,int *NumBlo,int NumSis,double *pos,double *Mod,string *FaseP,double *Lim,double &TC){
  double *t;
  int i;
  TC++;
  //cout<<TC<<endl;
  t=  new double [NEST];
  for(i=0;i<NEST;i++){
    if(FaseP[i]=="P")
      t[i]=(1/Mod[4])*sqrt((Mod[0]-pos[0*NEST+i])*(Mod[0]-pos[0*NEST+i])+(Mod[1]-pos[1*NEST+i])*(Mod[1]-pos[1*NEST+i])+(Mod[2]-pos[2*NEST+i])*(Mod[2]-pos[2*NEST+i]))+Mod[3];
    if(FaseP[i]=="SP")
      t[i]=(0.7/Mod[4])*sqrt((Mod[0]-pos[0*NEST+i])*(Mod[0]-pos[0*NEST+i])+(Mod[1]-pos[1*NEST+i])*(Mod[1]-pos[1*NEST+i])+(Mod[2]-pos[2*NEST+i])*(Mod[2]-pos[2*NEST+i]));
    if(FaseP[i]=="N")
      t[i]=0;
  }
  return t;
  delete [] t;
}

double *CalculeTA(int NEST,int *NumBlo,int NumSis,double *pos,double *Mod,string *Fase,double *Lim,double & TC){
  int NStation=0;
  int i,j,ip;
  for(i=0;i<NEST;i++) if(Fase[i]!="N") NStation++;
  
  int Dimmx, Dimmy,Dimmz,NDiv[3],NumBloT;
  int Lx,Ly,Lz;
  
  double Station[NStation][3],*VecVel,MinVel;
  double Source[3];
  double Dtt=0.08;
  double TLim,TLim1,TLim2,Per;
  double S0=0.04;
  double *T,*TTeo;
  
  T     =new double [NStation];
  TTeo  =new double [NEST];
  //NDiv  =new int    [3];
  
  NDiv[0]= NumBlo[0];
  NDiv[1]= NumBlo[1];
  NDiv[2]= NumBlo[2];
  
  NumBloT=NumBlo[0]*NumBlo[1]*NumBlo[2];
  VecVel=new double [NumBloT];
  Per=6;
  //cout<<NStation<<endl;
  //limite de tiempo en segundos  
  Dimmx=int(Lim[3]-Lim[0]);//Numero de divisiones
  Dimmy=int(Lim[4]-Lim[1]);//Numero de divisiones
  Dimmz=int(Lim[5]-Lim[2]);//Numero de divisiones
  
  Lx=int(Lim[3]-Lim[0]);//en Km
  Ly=int(Lim[4]-Lim[1]);//en Km
  Lz=int(Lim[5]-Lim[2]);//en Km
  
  for(i=0;i<NumBloT;i++){
    VecVel[i]=Mod[4+i];
    if     (i==0)             MinVel=VecVel[i];
    else if(MinVel>VecVel[i]) MinVel=VecVel[i];  
  }
  
  ip=0;
  for(i=0;i<NEST;i++){
    if(Fase[i]!="N"){
      if(ip==1) TLim2=TLim1;
      TLim1=(1.0/MinVel)*sqrt((Mod[0]-pos[0*NEST+i])*(Mod[0]-pos[0*NEST+i])+(Mod[1]-pos[1*NEST+i])*(Mod[1]-pos[1*NEST+i])+(Mod[2]-pos[2*NEST+i])*(Mod[2]-pos[2*NEST+i]));
      if( i!=0 && TLim1>TLim2) TLim2=TLim1;
      ip++;
    }
  }
  //TC+=TLim2*5.015;
  TC++;
  //cout<<TLim2*5.015<<"  "<<TC/60<<endl;
  
  cout<<TC<<"  "<<TLim2<<endl;
  TLim=TLim2+3;
  
  
  for(j=0;j<3;j++)
    Source[j]=Mod[j]-Lim[j];    
  
  ip=0;
  for(i=0;i<NEST;i++)
    if(Fase[i]!="N"){
      for(j=0;j<3;j++){Station[ip][j]=pos[i+j*NEST]-Lim[j];}
      ip++;
    }
  //Process(NStation,Station,Source,Dimmx,Dimmy,Dimmz,Dtt,TLim,vel,Lx,Ly,Lz,S0,VecVel,NCubos,NDiv,Per);
  Process(NStation,Station,Source,Dimmx,Dimmy,Dimmz,Dtt,TLim,Lx,Ly,Lz,S0,VecVel,NDiv[0]*NDiv[1]*NDiv[2],NDiv,Per);
  //Process(NStation,Station,Source,Dimmx,Dimmy,Dimmz,Dtt,TLim,vel,Lx,Ly,Lz,S0);
  T=Archivos(NStation);//calcula los tiempos
   
  ip=0;
  for(i=0;i<NEST;i++){
    if(Fase[i]=="P")  {TTeo[i]=T[ip]+Mod[3];      ip++;}
    if(Fase[i]=="SP") {TTeo[i]=(sqrt(3)-1)*T[ip]; ip++;}
    if(Fase[i]=="N")  {TTeo[i]=0;                      }
  }
  
  return TTeo;
  delete []T;
  delete []TTeo;
}

double * CalculeTT(int NumEstC,int *NumBlo,int NumSis,double *Pos,double *Mod,string *Fase, double *L,double & TC) {
  
  int i,j,k,NumBloT;
  double *ModP,*tteoP,*tteo;
  string *FaseP;
  NumBloT=NumBlo[0]*NumBlo[1]*NumBlo[2];
  ModP   =new double [4+NumBloT];
  tteoP  =new double [NumEstC];
  tteo   =new double [NumEstC*NumSis];
  FaseP  =new string [NumEstC];
  
  for(k=0;k<NumSis;k++) {
    for(i=0;i<NumEstC;i++) FaseP[i]=Fase[i+NumEstC*k];
    for(j=0;j<4+NumBloT;j++)
      if (j<4) ModP[j]=Mod[4*k+j];
      else     ModP[j]=Mod[4*(NumSis-1)+j];
      tteoP=CalculeTA(NumEstC,NumBlo,NumSis,Pos,ModP,FaseP,L,TC);
      for(i=0;i<NumEstC;i++) {tteo[NumEstC*k+i]=tteoP[i];}
  }
  
    return tteo;
    delete [] ModP;
    delete [] tteoP;
    delete [] tteo;
    delete [] FaseP;
    
}

double * CalculeG2(int NEST,int *NumBlo,int NumSis,double *pos,double *Mod,string *Fase,double *Lim, double & TC){
  
  int i,j,j1,k,l,y,ki,kj,e,NumFil,NumCol,NumBloT;
  double *a,*c,*g,*Mod1,*Mod2,*DMod,DV,*tteo1,*tteo2,V1,V2;
  string *FaseP;
  NumBloT=NumBlo[0]*NumBlo[1]*NumBlo[2];
  NumFil =NEST*NumSis;
  NumCol =4*NumSis;
  
  a    =new double [4*NEST*NumSis];
  //c    =new double [NEST*NumSis*NumBloT];
  g    =new double [NumFil*NumCol];
  Mod1 =new double [4+NumBloT];
  Mod2 =new double [4+NumBloT];
  DMod =new double [3];
  tteo1=new double [NEST];
  tteo2=new double [NEST];
  FaseP=new string [NEST];
  
  for(j=0;j<3;j++) DMod[j]=4;
  
  
  for (k=0;k<NumSis;k++){  
    for(j=0;j<4;j++){
      
      for(j1=0;j1<4+NumBloT;j1++){
	if(j1<3)
	  if(j1==j)     {Mod1[j1]=Mod[4*k+j1]-DMod[j1]/2;          Mod2[j1]=Mod[4*k+j1]+DMod[j1]/2;         }
	  else          {Mod1[j1]=Mod[4*k+j1];                     Mod2[j1]=Mod[4*k+j1];                    }
	else if (j1==3) {Mod1[j1]=Mod[4*k+j1];                     Mod2[j1]=Mod[4*k+j1];                    }	
	else            {Mod1[j1]=Mod[4*(NumSis-1)+j1];            Mod2[j1]=Mod[4*(NumSis-1)+j1];           }	   
      }
      
      for(i=0;i<NEST;i++) FaseP[i]=Fase[i+NEST*k];
      
      if(j!=3){
	tteo1=CalculeTA(NEST,NumBlo,NumSis,pos,Mod1,FaseP,Lim,TC);
	tteo2=CalculeTA(NEST,NumBlo,NumSis,pos,Mod2,FaseP,Lim,TC);
      }
      
      for(i=0;i<NEST;i++){
	if(j<3 )           a[(4*k+j)*NEST+i]      =(tteo2[i]-tteo1[i])/DMod[j];
	if(j==3){
	  if(Fase[i]=="P") a[(4*k+j)*NEST+i]      =1;
	  else             a[(4*k+j)*NEST+i]      =0;
	}
      }
    }
  }
  
  
  for(y=0;y<NumFil*NumCol;y++){
    i=y%NumFil;
    j=int(y/NumFil);
    ki=int(i/NEST);
    kj=int(j/4);
    e=y%NEST;
    if  (ki!=kj) g[y]=0;
    else         g[y]=a[int((y-e-ki*NEST)/NumSis)+e];
  }
  
  
  
  return g;
  delete [] a;
  delete [] c;
  delete [] g;
  delete [] Mod1;
  delete [] Mod2;  
  delete [] DMod;
  delete [] tteo1;
  delete [] tteo2;
  
}
  
