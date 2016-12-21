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
  const char *ArchPic ="DIR_DAT2";                             //Archivo de Datos
  const char *ArchEst ="DIR_PosEst";                               //Archivo de Estaciones
  const char *ArchMod ="DIR_MODI2";                               //Archivo con Modelo Inicial
  
  //const char *ArchSal  ="DIR_ITER";                             //Archivo Salida. Iteraciones
  //const char *ArchSal2 ="DIR_MODF";                             //Archivo Salida. Modelo Final
  const char *ArchSal3 ="DIR_TTEO";                             //Archivo Salida. Tiempos Teoricos
  const char *ArchMap ="DIR_MAP";                              //Archivo Salida. Mapa de localizacion    
  const char *ArchMap2="DIR_MAPTR";                              //Archivo Salida. Mapa de trazado de rayos    
  //_____________________________________________ DEFINICION DE VARIABLES________________________________________________________//
  
  int      i,j,k,*NumBlo,NumBloT,NumSis,NumEst,NumFil,NumCol,NumPicT,NumEstC,iteracion,*Fecha;
  int      NumEstCP,NumEstCSP,NumPicTP,NumPicTSP,*NumPicP,*NumPicSP;
  Sistema  A;
  double   e1,e2,De,D2e,V,TC,e,lim;
  double   *Mod1,*Mod2,*ModP,*TObs,*Pos,*PosT,*L;
  double   *tres,*TTeo,*TTeoP,*TGG,*TGtres,*DMod,*dato,*G,*TG,*I, *E,*TGG2;
  //ofstream FILE(ArchSal);
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
 
  
  De=1e30; e1=1e30;
  iteracion=0;
  TC=0;
  //_______________________________________________CALCULA TIEMPO E IMPRESION____________________________________________________//
  TTeo=CalculeTT(NumEstC,NumBlo,NumSis,Pos,Mod1,Fase,L,TC);         //Calcula tiempo teorico
  ImpresionTTeo(ArchSal3,NumEstC,NumSis,TTeo,Fase,NomEstC,Fecha);
  AnimacionEstaciones(ArchMap,ArchMap2,L,PosT,NumEst,NomEstT,NumSis,Mod1,NumEstC,Fase,Pos);                 //Imprime archivo gnuplot con el mapa de estaciones y localizacion
  
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

  //FILE.close();

  //___________________________________________________________FIN________________________________________________________________//

    return 0;  
}
