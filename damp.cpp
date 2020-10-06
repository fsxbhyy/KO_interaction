#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>


void Gap::DampIteration(const int& counter){ 
        int i,j;
	int mm=10;
  if(counter==1){
    mm=5;
  }
  // else if(counter==2){
  //   mm=1;
  // }
	i=0;
	std::vector<double> kH1(delta_H.size(),0);
	kH1=delta_H;
	for(i=0;i<mm;i++){
		MOperate(3);
		for(j=0;j<delta_H.size();j++){
			delta_H[j]=delta_H[j]+inhomoH[j];
                	kH1[j]=(kH1[j]*(i+1)+delta_H[j])/(i+2);
		}
		delta_H=kH1;
	}
}


