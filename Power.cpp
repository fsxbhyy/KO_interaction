#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>

/*double Gap::Power(){
        int i,j;
	int mm=10;
	int size_f_l=delta_L.size();
	shift=0;
	double lamu=0;
	double modulus=0;
	std::vector<double> initial(size_f_l,0);
	initial=delta_L;
	std::vector<double> kL(size_f_l,0);
	
        for(i=0;i<10;i++){
		modulus=sqrt(Inner(delta_L,delta_L));
		for(j=0;j<size_f_l;j++){
			delta_L[j]=delta_L[j]/modulus;
		}
	//	for(j=0;j<alpha.NNL;j++){
          //                      printf("%f 1\n",alpha.inhomoL[j]);
            //    }
		 kL=delta_L;
                MOperate(0);
	//	for(j=0;j<alpha.NNL;j++){
          //                      printf("%f 2\n",alpha.Low[j]);
            //    }

		for(j=0;j<size_f_l;j++){
			delta_L[j]=delta_L[j]+inhomoL[j];
		}
                lamu=Inner(delta_L,kL);
       // 	printf("lamu=%f\n",lamu);
	}
	

	if(lamu<0){
        	shift=-lamu;
		delta_L.swap(initial);
		for(i=0;i<mm;i++){
			modulus=sqrt(Inner(delta_L,delta_L));
			for(j=0;j<size_f_l;j++){
				delta_L[j]=delta_L[j]/modulus;
			}
			kL=delta_L;
			MOperate(0);
			for(j=0;j<size_f_l;j++){
				delta_L[j]=delta_L[j]+shift*kL[j]+inhomoL[j];
			}
			lamu=Inner(delta_L,kL);
			 printf("lamushift=%f\n",lamu-shift);
		}
	}
	
	modulus=sqrt(Inner(delta_L,delta_L));
        for(j=0;j<size_f_l;j++){
                delta_L[j]=delta_L[j]/modulus;
        }


	return lamu-shift;
}*/


double Gap::Power(){
        int i,j;
	int mm=100;
	int size_f_l=delta_L.size();
	double lamu=0;
	double modulus=0;
	std::vector<double> kL(size_f_l,0);
if(shift==0){	
        std::vector<double> initial(size_f_l,0);
	initial=delta_L;
	for(i=0;i<mm;i++){
		modulus=sqrt(Inner(delta_L,delta_L));
		for(j=0;j<size_f_l;j++){
			delta_L[j]=delta_L[j]/modulus;
		}
	//	for(j=0;j<alpha.NNL;j++){
          //                      printf("%f 1\n",alpha.inhomoL[j]);
            //    }
		 kL=delta_L;
                MOperate(0);
	//	for(j=0;j<alpha.NNL;j++){
          //                      printf("%f 2\n",alpha.Low[j]);
            //    }

		for(j=0;j<size_f_l;j++){
			delta_L[j]=delta_L[j]+inhomoL[j];
		}
                lamu=Inner(delta_L,kL);
        	printf("lamu=%f\n",lamu);
	}
/*	if(lamu<0){
        	shift=shift-lamu;
		delta_L.swap(initial);
	}*/
}
	
	if(shift!=0){
		
		for(i=0;i<mm;i++){
			modulus=sqrt(Inner(delta_L,delta_L));
			for(j=0;j<size_f_l;j++){
				delta_L[j]=delta_L[j]/modulus;
			}
			kL=delta_L;
			MOperate(0);
			for(j=0;j<size_f_l;j++){
				delta_L[j]=delta_L[j]+shift*kL[j]+inhomoL[j];
			}
			lamu=Inner(delta_L,kL);
			 printf("lamushift=%f\n",lamu-shift);
		
		}
	//	if(lamu<0) shift=shift-lamu;
		
	}
	
	modulus=sqrt(Inner(delta_L,delta_L));
        for(j=0;j<size_f_l;j++){
                delta_L[j]=delta_L[j]/modulus;
        }


	return lamu-shift;
}

