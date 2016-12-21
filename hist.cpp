#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "INV3.h"

using namespace std;

void Histograma(const char *ArchIn,const char *ArchOut,double x1,double x2,int bin){
  int i,*Sum,Nbin;
  double *D,Dbin;
  string *DS;
  ifstream FILEIn (ArchIn);
  ofstream FILEOut(ArchOut);
  
  DS     =new string[53];
  D      =new double[53];
  Sum    =new int [bin];

  DS=AdquieraDatos2(ArchIn,53,1,0);
  D=ConviertaDouble(DS,53);
  
  Dbin=(x2-x1)/bin;
  for(i=0;i<bin;i++) Sum[i]=0;
  for(i=0;i<53;i++){
    Nbin=int(D[i]/Dbin);
    Sum[Nbin]++;
  }
  for(i=0;i<bin;i++) FILEOut<<x1+i*Dbin<<"  "<<Sum[i]<<endl;
  
  
}


int main(){
  const char *ArchIn   ="LOC_SISMOS_PROF";
  const char *ArchOut  ="HIST";
  int bin;
  double x1,x2;
  bin =7;
  x1=0;
  x2=70;
  Histograma(ArchIn,ArchOut,x1,x2,bin);

  return 0;
}
