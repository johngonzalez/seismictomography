#include <iostream>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <sstream>

using namespace std;

double Z(double x, double y){
  return x*y;
}

int main(){
  int i,j;
  double x,y;
  ofstream FILE("d.dat");
  
  for(i=0;i<100;i++){
    x=i*0.1;
    for(j=0;j<100;j++){
      y=j*0.1;
      cout<<x<<" "<<y<<" "<<Z(x,y)<<endl;
    }
    cout<<endl;
  }
  FILE.close();
  
  return 0;
}
