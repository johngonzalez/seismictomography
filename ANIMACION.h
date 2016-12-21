#include <iostream>
#include <cmath>
#include <stdio.h>

void AnimacionInicial(double *loc,double *pos,double *tobs,int NEST,int NumPicP,string *NomEst,double V){
  int i,j;
  double D[NEST];
  
  for(i=NumPicP;i<NEST;i++){
    D[i]=V/0.7*tobs[i];
    if(i==NumPicP)
      cout<<"plot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
    else
      cout<<"replot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
  }
  
  cout<<"pause 0.1"<<endl; 
}

void Animacion(double *loc,double *pos,double *tobs,int NEST,int NumPicP,double error,string *NomEst,double V){
  int i;
  double D[NEST];
  
  cout<<"set title ' e="<<int(error)<<"  x="<<int(loc[0])<<"  y="<<int(loc[1])<<"  z="<<int(loc[2])<<"  t="<<int(loc[3])<<"  v="<<V<<"'"<<endl;
  for(i=0;i<NEST;i++){
    if(i<NumPicP)
      D[i]=V*(tobs[i]-loc[3]);
    else
      D[i]=V/0.7*tobs[i];
    if(i==0)
      cout<<"plot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
    else
      cout<<"replot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
  }
  
  cout<<"replot "<<2<<"*cos(t)+"<<loc[0]<<","<<2<<"*sin(t)+"<<loc[1]<<" lc 9 t 'Localizacion' "<<endl;
  cout<<"pause 2"<<endl; 
  
}
void AnimacionFinal(double *loc,double *pos,double *tobs,int NEST,int NumPicP,double error,string *NomEst,double V){
  int i;
  double D[NEST];
  cout<<"set title ' e="<<int(error)<<"  x="<<int(loc[0])<<"  y="<<int(loc[1])<<"  z="<<int(loc[2])<<"  t="<<int(loc[3])<<"  v="<<V<<"'"<<endl;
  
  for(i=0;i<NEST;i++){
    if(i<NumPicP)
      D[i]=V*(tobs[i]-loc[3]);
    else
      D[i]=V/0.7*tobs[i];
    if(i==0)
      cout<<"plot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
    else
      cout<<"replot "<<D[i]<<"*cos(t)+"<<pos[i+0*NEST]<<","<<D[i]<<"*sin(t)+"<<pos[i+1*NEST]<<" t '"<<NomEst[i]<<"'"<<endl;
  }
  
  //cout<<"replot "<<2<<"*cos(t)+"<<loc[0]<<","<<2<<"*sin(t)+"<<loc[1]<<" ps 1.5 t 'Localizacion' "<<endl;
  cout<<"replot '<echo "<<loc[0]<<"  "<<loc[1]<<"' pt 7 ps 3 lc 0 title 'Localizacion'"<<endl
    //<<"set term postscript"<<endl
    //<<"set output 'MapaCirculos2.png'"<<endl
      <<"replot"<<endl
      <<"pause 10"<<endl; 
  
}
void AnimacionEstaciones(const char *ArchMap,const char *ArchMap2,double *Lim,double *PosT,int NEST,string *NomEst,int NumSis,double *Mod,int NumEstC,string *Fase,double * Pos){

  int i,k;
  ofstream FILE(ArchMap);
  ofstream FILE2(ArchMap2);
  //_________________________________________________________________________________________________________________________________

  FILE<<"set size square"<<endl
      <<"set xrange["<<Lim[0]<<":"<<Lim[3]<<"]"<<endl
      <<"set yrange["<<Lim[1]<<":"<<Lim[4]<<"]"<<endl
      <<"set key out vert"<<endl
      <<"set xlabel 'Posicion x (Km)'"<<endl
      <<"set ylabel 'Posicion y (Km)'"<<endl
      <<"unset key"<<endl;
  
  for(i=0;i<NEST;i++)
    FILE<<"set label "<<i+1<<" '"<<NomEst[i]<<"' "<<"at "<<2+PosT[i+0*NEST]<<", "<<2+PosT[i+1*NEST]<<", "<<0<<" nopoint tc def"<<endl;
  
  for(i=0;i<NEST;i++)
    if (i==0)
      FILE<<"plot '<echo "<<PosT[i+0*NEST]<<"  "<<PosT[i+1*NEST]<<"  "<<0<<"' pt 9 ps 1 lc 0"<<endl;
    else
      FILE<<"replot '<echo "<<PosT[i+0*NEST]<<"  "<<PosT[i+1*NEST]<<"' pt 9 ps 1 lc 0"<<endl;
  
  for(k=0;k<NumSis;k++)
    FILE<<"replot '<echo "<<Mod[0+4*k]<<"  "<<Mod[1+4*k]<<"' pt 7 ps 2 lc 1"<<endl;      
  
  FILE<<"pause 60"<<endl;
  FILE.close();
  //________________________________________________________________________________________________________________________________
  
  FILE2<<"set size square"<<endl
       <<"set xrange["<<Lim[0]<<":"<<Lim[3]<<"]"<<endl
       <<"set yrange["<<Lim[1]<<":"<<Lim[4]<<"]"<<endl
       <<"set key out vert"<<endl
       <<"set xlabel 'Posicion x (Km)'"<<endl
       <<"set ylabel 'Posicion y (Km)'"<<endl
       <<"unset key"<<endl
       <<"set style arrow 7 nohead ls 1"<<endl;
  for(i=0;i<NEST;i++)
    FILE2<<"set label "<<i+1<<" '"<<NomEst[i]<<"' "<<"at "<<2+PosT[i+0*NEST]<<", "<<2+PosT[i+1*NEST]<<", "<<0<<" nopoint tc def"<<endl;
  
  for(k=0;k<NumSis;k++)
    for(i=0;i<NumEstC;i++)
      if(Fase[NumEstC*k+i]!="N")
	FILE2<<"set arrow from "<<Mod[0+4*k]<<","<<Mod[1+4*k]<<" to "<<Pos[NumEstC*0+i]<<","<<Pos[NumEstC*1+i]<<" as 7"<<endl;
  
  for(i=0;i<NEST;i++)
    if (i==0)
      FILE2<<"plot '<echo "<<PosT[i+0*NEST]<<"  "<<PosT[i+1*NEST]<<"  "<<0<<"' pt 9 ps 1 lc 0"<<endl;
    else
      FILE2<<"replot '<echo "<<PosT[i+0*NEST]<<"  "<<PosT[i+1*NEST]<<"' pt 9 ps 1 lc 0"<<endl;
  
  for(k=0;k<NumSis;k++)
    FILE2<<"replot '<echo "<<Mod[0+4*k]<<"  "<<Mod[1+4*k]<<"' pt 7 ps 2 lc 1"<<endl;              
  FILE2<<"pause 60"<<endl;
  FILE2.close();

}
  
void Impresion(const char *ArchSal2,int iteracion,double e,double lim,double e1,double De,int NumSis,int NumBlo,double *Mod1){
  int i,j,k;
  ofstream FILE(ArchSal2);
  FILE<<e<<" "<<lim<<endl<<endl
      <<iteracion<<"  "<<e1<<" "<<abs(De)<<endl<<endl;
  
  for(k=0;k<NumSis;k++){
    for(i=0;i<4;i++) FILE<<Mod1[4*k+i]<<"  ";
    FILE<<endl;
  }
  for(j=0;j<NumBlo;j++) FILE<<Mod1[4*NumSis+j]<<"  ";
  FILE<<endl<<endl;
  FILE.close();
}
void Impresion2(int iteracion,double e1,double De,int NumSis,int NumEstC,int *Fecha,double *Mod1,string *NomEstC,string *Fase,double *TObs,double *TTeo){
  int i,k;
  
  cout<<endl<<iteracion<<"  "<<e1<<" "<<abs(De)<<endl<<endl;
  for(k=0;k<NumSis;k++){
    for(i=0;i<3      ;i++) cout<<Fecha[NumSis*i+k]<<" "; cout<<endl;
    for(i=0;i<4      ;i++) cout<<Mod1[4*k+i]<<"  "; cout<<endl<<endl;
    for(i=0;i<NumEstC;i++) cout<<NomEstC [i]<<" "<<Fase[NumEstC*k+i]<<" "<<TObs[NumEstC*k+i]<<" "<<TTeo[NumEstC*k+i]<<" "<<TObs[NumEstC*k+i]-TTeo[NumEstC*k+i]<<endl;
    cout<<endl;
  }
}
void ImpresionTTeo(const char *ArchSal3,int NumEstC,int NumSis,double *TTeo,string *Fase,string *NomEstC,int *Fecha){
  int i,k;
  ofstream FILE(ArchSal3);
  
  for(k=0;k<NumSis;k++){
    FILE<<"SISMO "<<Fecha[NumSis*0+k]<<" "<<Fecha[NumSis*1+k]<<" "<<Fecha[NumSis*2+k]<<endl;
    for(i=0;i<NumEstC;i++){
      if (Fase[i+k*NumSis]!="N")
	FILE<<NomEstC[i]<<" "<<Fase[i+k*NumEstC]<<" "<<TTeo[i+k*NumEstC]<<endl;
    }
  }
  FILE.close();
}

