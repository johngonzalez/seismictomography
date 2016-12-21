/*
  PROGRAMA PARA INVERTIR HIPOCENTROS SIMULTANEAMENTE Y MODELO DE VELOCIDADES A PARTIR DE PICADA DE PRIMEROS ARRIBOS P y SP.

  ARCHIVOS DE ENTRADA: Picada de tiempos con fase P y SP.
                       Archivo con informacion de la posicion de las estaciones
                       Modelo inicial. Hipocentros y velocidades
  ARCHIVO DE SALIDA:   Informacion de cada iteracion. Modelo final. Hipocentro y velocidades
*/

#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include "SL.h"        //Libreria para resolver sistema de ecuaciones
#include "INV3.h"      //Libreria con herramientas de inversion
#include "MD5.h"       //Libreria para resolver modelo directo (Diferencias Finitas)
#include "ANIMACION.h" //Libreria para impresion y animacion

using namespace std;

int main() {
  
  //____________________________________________INICIO.ARCHIVOS DE LECTURA_________________________________________________________//

  const char *ArchPar ="DIR_CAL";                             //Archivo de parametros
  const char *ArchPic ="DIR_TTEO";                             //Archivo de Datos
  const char *ArchEst ="DIR_PosEst";                               //Archivo de Estaciones
  const char *ArchMod ="DIR_MODI4";                               //Archivo con Modelo Inicial
  
  const char *ArchSal  ="INV_SIN_ITER";                             //Archivo Salida. Iteraciones
  const char *ArchSal2 ="INV_SIN_MODF";                             //Archivo Salida. Modelo Final
  const char *ArchMap ="INV_SIN_MAP";                              //Archivo Salida. Mapa de localizacion    
  const char *ArchMap2="INV_SIN_MAPTR";                              //Archivo Salida. Mapa de trazado de rayos    
  //_____________________________________________ DEFINICION DE VARIABLES________________________________________________________//
  
  int      i,j,k,*NumBlo,NumBloT,NumSis,NumEst,NumFil,NumCol,NumPicT,NumEstC,iteracion,*Fecha;
  int      NumEstCP,NumEstCSP,NumPicTP,NumPicTSP,*NumPicP,*NumPicSP;
  Sistema  A;
  double   e1,e2,De,D2e,V,TC,e,lim;
  double   *Mod1,*Mod2,*ModP,*TObs,*Pos,*PosT,*L;
  double   *tres,*TTeo,*TTeoP,*TGG,*TGtres,*DMod,*dato,*G,*TG,*I, *E,*TGG2;
  ofstream FILE(ArchSal);
  string   *NomEstP,*NomEstSP,*NomEstT,*NomEstC,*NomEstCP,*NomEstCSP,*Fase;
  
  //_________________________________________________ADQUISION DE VALORES_______________________________________________________//

  L =new double [3*2];
  VariablesControl(ArchPar,e,lim,L);                                          //Adquiere Variables de Control
  
  NumEst   =Numero(ArchEst,1);                                                //Adquiere Numero Estaciones
  PosT     =new double [NumEst*3];                    
  NomEstT  =NombreEstacionesTodas  (ArchEst);                                 //Adquiere Nombre de todas las Estaciones
  PosT     =PosicionEstacionesTodas(ArchEst);                                 //Adquiere Posicion de todas las Estaciones
  //InicieAnimacion(L[0],L[3],L[1],L[4]);
  //AnimacionEstaciones(PosT,NumEst,NomEstT);
  NumeroSismosPicadas(ArchPic,NumSis,NumPicT);                                //Adquiere Numero de Sismos y numero de picadas en total
  NumeroPicadasFase(ArchPic,NumPicT,NumPicTP,NumPicTSP);                      //Adquiere Numero de picadas por cada fase
  
  NumEstCP =NumeroEstacionesConsideradasFase(ArchPic,NumPicT,"P");            //Adquiere Numero de Estaciones consideradas con fase P
  NumEstCSP=NumeroEstacionesConsideradasFase(ArchPic,NumPicT,"SP");           //Adquiere Numero de Estaciones consideradas con fase SP
  NumEstC= NumEstCP+NumEstCSP;
  
  NumBlo    =new int [3];
  NumBlo    =NumeroBloques(ArchMod);                                               //Adquiere numero de bloques
  NumBloT   =NumBlo[0]*NumBlo[1]*NumBlo[2];
  NumFil    =NumEstC*NumSis;                                                  //Calcula numero de filas (ecuaciones) de la matriz sensibilidad
  NumCol    =4*NumSis+NumBloT;                                                 //Calcula numero de columnas (incognitas) de la matriz sensibilidad
  cout<<NumBloT<<endl;
  //_____________________________________________GUARDADO EN MEMORIA DE VARIABLES_________________________________________________//
  
  NomEstP   =new string [NumPicTP];
  NomEstSP  =new string [NumPicTSP];
  NomEstCP  =new string [NumEstCP];
  NomEstCSP =new string [NumEstCSP];
  NomEstC   =new string [NumEstC];
  
  TObs      =new double [NumFil];
  Fase      =new string [NumFil];
  Pos       =new double [NumEstC*3];
  Mod1      =new double [NumCol];
  Mod2      =new double [NumCol];
  TTeo      =new double [NumFil];
  NumPicP   =new int    [NumSis];
  NumPicSP  =new int    [NumSis];
  Fecha     =new int    [NumSis*3];
  
  G         =new double [NumFil*NumCol];
  TG        =new double [NumCol*NumFil];
  TGG       =new double [NumCol*NumCol];
  I         =new double [NumCol*NumCol];  
  E         =new double [NumCol*NumCol];  
  TGG2      =new double [NumCol*NumCol];
  TGtres    =new double [NumCol];
  DMod      =new double [NumCol];
  tres      =new double [NumFil];
  
  //___________________________________________________ADQUISICION DE VALORES_____________________________________________________//
  
  NumPicP   =NumeroEstacionesSismoFase(ArchPic,NumSis,"P");                        //Adquiere numero de picadas P  por cada sismo
  NumPicSP  =NumeroEstacionesSismoFase(ArchPic,NumSis,"SP");                       //Adquiere numero de picadas SP por cada sismo
  
  Fecha     =NumeroFecha(ArchPic,NumSis);                                          //Adquiere Fecha (ano,mes,dia) de cada sismo
  NomEstCP  =EstacionesConsideradasFase(ArchPic,NumEstC,"P");                      //Adquiere nombre de estaciones consideradas P
  NomEstCSP =EstacionesConsideradasFase(ArchPic,NumEstC,"SP");                     //Adquiere nombre de estaciones consideradas SP
  NomEstC   =EstacionesConsideradas(NumEstCP,NumEstCSP,NomEstCP,NomEstCSP);        //Adquiere nombre de estaciones consideradas
  
  NomEstP  =NombreEstacionesSismoFase(ArchPic,NumPicT,"P");                        //Adquiere Nombre de las estaciones con fase P
  NomEstSP  =NombreEstacionesSismoFase(ArchPic,NumPicT,"SP");                      //Adquiere Nombre de las estaciones con fase SP
  
  Observacion(ArchPic,NumSis,NumPicTP,NumPicTSP,NumEstCP,NumEstCSP,NumPicP,NumPicSP,NomEstCP,NomEstCSP,NomEstP,NomEstSP,TObs,Fase); 
                                                                                   //Adquiere Tiempo y la fase de cada picada 
  Pos    =PosicionEstaciones(ArchEst,NomEstC,NumEst,NumEstC);                      //Adquiere la posicion de las estaciones
  Mod1   =ModeloInicial(ArchMod,NumSis,NumBloT);                                           //Adquiere hipocentros y modelo de velocidad inicial
  //AnimacionInicial(Mod1,Pos,Tobs,NumPic,NumPicP,NomEst,V);
    
  De=1e30; e1=1e30; iteracion=0; TC=0;
  
  FILE<<e<<" "<<lim<<endl<<endl;
  for(k=0;k<NumSis;k++){
    for(i=0;i<4;i++) FILE<<Mod1[4*k+i]<<"  ";
    FILE<<endl;
  }
  for(j=0;j<NumBloT;j++) FILE<<Mod1[4*NumSis+j]<<"  ";
  FILE<<endl<<endl;
  //________________________________________________________ITERACIONES____________________________________________________________//
  
  while(abs(De)>lim){
    //while(iteracion<1e30){
    
    //______________________IMPRESION______________________________//
    
    //Impresion(iteracion,e1,De,NumSis,NumBlo,Mod1);
    FILE<<iteracion<<"  "<<e1<<" "<<abs(De)<<endl;
    
    //_____________ARMAR MATRIZ Y RESOLVER SISTEMA________________//

    e2=e1;
    Mod2=Mod1;
    
    Mod1=ModeloDentroLimites(NumSis,Mod1,L);                          //Controla que el modelo de entrada no se salga de limites
    TTeo=CalculeTT(NumEstC,NumBlo,NumSis,Pos,Mod1,Fase,L,TC);         //Calcula tiempo teorico
    for(i=0;i<NumFil;i++) tres[i]=TObs[i]-TTeo[i];                    //Calcula tiempo residual
    G=CalculeG2(NumEstC,NumBlo,NumSis,Pos,Mod1,Fase,L,TC);            //Calcula matriz de sensibilidad
    TG=Transpuesta(G,NumFil,NumCol);                                  //Calcula matriz transpuesta
    TGG=Multiplique(TG,NumCol,NumFil,G,NumFil,NumCol);                //Multiplica matriz transpuesta por matriz (queda cuadrada)
    TGtres=Multiplique(TG,NumCol,NumFil,tres,NumFil,1);               //Multiplica matriz transpuesta por tiempo residual
    
    I=Identidad(NumCol);                                              //Matriz identidad
    E=MulEscalar(e,I,NumCol,NumCol);                                  //Matriz whitening
    TGG2=Sume(TGG,E,NumCol,NumCol);                                   //Matriz principal
    
    A.Inicie(TGG2,TGtres,NumCol);                                     //Inicia sistema Ax=b
    A.LU();                                                           //Resuelve sistema por LU
    DMod=A.Muestre();                                                 //Resultado (cambio de modelo, hipocentros y velocidad)
    
    for(i=0;i<NumCol;i++) Mod1[i]=Mod1[i]+DMod[i];                    //Calcula nuevo modelo
    e1=norma2(tres,NumFil);                                           //Calcula error cuadratico
    De=e2-e1;                                                         //Calcula diferencias de error entre modelo anterior y actual
    iteracion++;                                                      //Calcula iteracion
    cout<<iteracion<<"  "<<e1<<" "<<abs(De)<<"      "<<Mod1[4*NumSis]<<endl<<endl;
    AnimacionEstaciones(ArchMap,ArchMap2,L,PosT,NumEst,NomEstT,NumSis,Mod1,NumEstC,Fase,Pos);                 //Imprime archivo gnuplot con el mapa de estaciones y localizacion
    Impresion(ArchSal2,iteracion,e,lim,e1,De,NumSis,NumBloT,Mod1);           //Imprime informacion
  }
  
  //___________________________________________________IMPRESION FINAL__________________________________________________________//
  
  //Impresion2(iteracion-1,e1,De,NumSis,NumEstC,Fecha,Mod1,NomEstC,Fase,TObs,TTeo);
  FILE<<endl<<endl;
  if(e2<e1){
    FILE<<e2<<endl<<endl;
    for(k=0;k<NumSis;k++){
      for(i=0;i<4;i++)
	FILE<<Mod2[4*k+i]<<"  ";
      FILE<<endl;
    }
    for(j=0;j<NumBloT;j++) FILE<<Mod2[j]<<"  ";
    FILE<<endl;
  }
  else{
    FILE<<e1<<endl<<endl;
    for(k=0;k<NumSis;k++){
      for(i=0;i<4;i++)
	FILE<<Mod2[4*k+i]<<"  ";
      FILE<<endl;
    }
    for(j=0;j<NumBloT;j++) FILE<<Mod2[j]<<"  ";
    FILE<<endl;
  }
  
  
  //__________________________________________________BORRADO DE MEMORIA_________________________________________________________//
  
  delete [] L;
  delete [] PosT;
  delete [] NomEstCP;
  delete [] NomEstCSP;
  delete [] NomEstC;
  delete [] NumBlo;
  
  delete [] TObs;
  delete [] Fase;
  delete [] Pos;
  //delete [] Mod1;
  delete [] Mod2;
  delete [] TTeo;
  delete [] NumPicP;
  delete [] NumPicSP;
  delete [] Fecha;
 
  delete [] G;
  delete [] TG;
  delete [] TGG;
  delete [] I;
  delete [] E;
  delete [] TGG2;
  delete [] TGtres;
  delete [] DMod;
  delete [] tres;

  FILE.close();

  //___________________________________________________________FIN________________________________________________________________//

    return 0;  
}
