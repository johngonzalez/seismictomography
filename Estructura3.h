#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <sstream>
#include "Clase3.h"

using namespace std;

void Process(int NSt,double Stts[][3],double Sour[],int Dimmx, int Dimmy, int Dimmz,double Dtt, double TLim, int Lx, int Ly, int Lz, double S0, double VecVel[],int NCubos,int NDiv[],double T);

void Process(int NSt,double Stts[][3],double Sour[],int Dimmx, int Dimmy, int Dimmz,double Dtt, double TLim, int Lx, int Ly, int Lz, double S0, double VecVel[],int NCubos, int NDiv[],double T)
{ 
  
  double t=0, Amp=0;
  int p=0;
  int i;
  
  ofstream RegStation;
  ofstream InfoStat;   
  
  Matriz Mtz(Dimmx,Dimmy,Dimmz,Dtt,Lx,Ly,Lz,T);  
  RegStation.open("RStations.dat");
  InfoStat.open("InfoStation.dat");
  Mtz.CargueParam(Sour,S0);
  Mtz.Modelo();
  Mtz.Velocidad(VecVel,NCubos,NDiv);
  Mtz.Inicie(t);
  
  //Mtz.AnimacionInicie();
  
  for(t=0;t<=TLim;t=t+Dtt)
    {/*
      if(p%10==0)
	{
	  Mtz.Animacion(t);
	}
     */
      RegStation<<t;
      for(i=0;i<NSt;i++)
	{
	  Amp=Mtz.RegistroEstaciones(t,Stts[i][0],Stts[i][1],Stts[i][2]);
	  RegStation<<" "<<Amp;
	}
      RegStation<<endl;
      Mtz.Evolucion(t);
      Mtz.Cambio();
      p++;
    }
  InfoStat<<p<<" "<<Dtt<<endl;
  
  RegStation.close();
  InfoStat.close();
}
