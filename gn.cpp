#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "INV3.h"

using namespace std;

void ModeloVelocidad(const char *ArchVelIn,const char *ArchVelOut){
  int i,*NumBlo,NumBloT;
  double *L,*L1,*Vel;
  string *LS,*NumBloS,*VelS;
  ifstream FILEIn (ArchVelIn);
  ofstream FILEOut(ArchVelOut);
  
  LS     =new string[3*2];NumBloS=new string [3];
  L      =new double[3*2];
  NumBlo =new int [3];
  
  LS=AdquieraDatos2(ArchVelIn,3,2,0);
  L=ConviertaDouble(LS,3*2);
  
  NumBloS=AdquieraDatos2(ArchVelIn,1,3,3);
  NumBlo =ConviertaInt(NumBloS,3);
  NumBloT=NumBlo[0]*NumBlo[1]*NumBlo[2];
  
  Vel    =new double[NumBloT];
  VelS   =new string[NumBloT];
  VelS=AdquieraDatos2(ArchVelIn,1,NumBloT,4);
  Vel =ConviertaDouble(VelS,NumBloT);

  
}


int main(){
  const char *ArchVelIn   ="MOD_VEL";
  const char *ArchVelOut  ="MOD_VELGN";
  ModeloVelocidad(ArchVelIn,ArchVelOut);

  return 0;
}
