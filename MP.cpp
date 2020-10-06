#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>

void Gap::MOperate(const int& switF){
        int i,j,n,m;//dif is the label of w_n-w_m	
	int alpha;
	int size_q=Grid_q->size();
//	int size_f=Grid_f_l->size()+Grid_f_h->size();
	int size_f_l=Grid_f_l->size();
	int size_f_h=Grid_f_h->size();
	int size_f_t=size_f_l+size_f_h;
	double sum;
	double dif,dif1;
	double difq;
	double difqsum=0;
	double q_dum;
	if(switF==0){ //K11
		std::vector<double> kL;
	      		for(i=0;i<size_f_l;i++){
				for(j=0;j<size_q;j++){
					sum=0;
					for(n=0;n<size_f_l;n++){
						alpha=0;
						difqsum=0;
						for(m=0;m<size_q;m++){	
							while(alpha<(*W->Grid_int_q)[j].size()&&(*W->Grid_int_q)[j][alpha]<=(*Grid_q)[m]){
								q_dum=(*W->Grid_int_q)[j][alpha];
									
								if(alpha==0) difq=q_dum;	
								else difq=q_dum-(*W->Grid_int_q)[j][alpha-1];			
									
								if(m==0) dif1=delta_L[n*size_q+m];
								else dif1=delta_L[n*size_q+m-1]+(delta_L[n*size_q+m]-delta_L[n*size_q+m-1])/((*Grid_q)[m]-(*Grid_q)[m-1])*(q_dum-(*Grid_q)[m-1]);

								difqsum+=difq;
								sum-=(*GG->cgfactor_l)[n]*difq/((*Grid_f_l)[n]*(*Grid_f_l)[n]+(q_dum*q_dum-miu)*(q_dum*q_dum-miu))*Interpolate(i,n,j,alpha)*dif1;
								alpha++;
							}
						}
		 
					}
					sum=sum*Temp/4/M_PI/M_PI;
					kL.push_back(sum);
//					printf("%f %f %f\n",(*Grid_f_l)[i],(*Grid_q)[j],sum);
				}	      
			}
			delta_L.swap(kL);
	}

	else if(switF==1){ //K12
        		for(i=0;i<size_f_l;i++){
				for(j=0;j<size_q;j++){
					sum=0;
					for(n=0;n<size_f_h;n++){
						alpha=0;
						for(m=0;m<size_q;m++){	
							while(alpha<(*W->Grid_int_q)[j].size()&&(*W->Grid_int_q)[j][alpha]<=(*Grid_q)[m]){
								q_dum=(*W->Grid_int_q)[j][alpha];
									
								if(alpha==0) difq=q_dum;	
								else difq=q_dum-(*W->Grid_int_q)[j][alpha-1];			
									
								if(m==0) dif1=delta_H[n*size_q+m];
								else dif1=delta_H[n*size_q+m-1]+(delta_H[n*size_q+m]-delta_H[n*size_q+m-1])/((*Grid_q)[m]-(*Grid_q)[m-1])*(q_dum-(*Grid_q)[m-1]);

								sum-=(*GG->cgfactor_h)[n]*difq/((*Grid_f_h)[n]*(*Grid_f_h)[n]+(q_dum*q_dum-miu)*(q_dum*q_dum-miu))*Interpolate(i,n+size_f_l,j,alpha)*dif1;
								alpha++;
							}
						}
						 
					}
				sum=sum*Temp/4/M_PI/M_PI;
				inhomoL[i*size_q+j]=sum;
				}	      
			}
	}
	else if(switF==2){ //K21
 
		        for(i=0;i<size_f_h;i++){
				for(j=0;j<size_q;j++){
					sum=0;
					for(n=0;n<size_f_l;n++){
						alpha=0;
						for(m=0;m<size_q;m++){	
							while(alpha<(*W->Grid_int_q)[j].size()&&(*W->Grid_int_q)[j][alpha]<=(*Grid_q)[m]){
								q_dum=(*W->Grid_int_q)[j][alpha];
									
								if(alpha==0) difq=q_dum;	
								else difq=q_dum-(*W->Grid_int_q)[j][alpha-1];			
									
								if(m==0) dif1=delta_L[n*size_q+m];
								else dif1=delta_L[n*size_q+m-1]+(delta_L[n*size_q+m]-delta_L[n*size_q+m-1])/((*Grid_q)[m]-(*Grid_q)[m-1])*(q_dum-(*Grid_q)[m-1]);

								sum-=(*GG->cgfactor_l)[n]*difq/((*Grid_f_l)[n]*(*Grid_f_l)[n]+(q_dum*q_dum-miu)*(q_dum*q_dum-miu))*Interpolate(i+size_f_l,n,j,alpha)*dif1;
								alpha++;
							}
						}
												
					}
				sum=sum*Temp/4/M_PI/M_PI;
				inhomoH[i*size_q+j]=sum;
				}	      
			}

        }
	else {	//K22
		std::vector<double> kH;
	        	for(i=0;i<size_f_h;i++){
				for(j=0;j<size_q;j++){
					sum=0;
					for(n=0;n<size_f_h;n++){
						alpha=0;	
						for(m=0;m<size_q;m++){	
							while(alpha<(*W->Grid_int_q)[j].size()&&(*W->Grid_int_q)[j][alpha]<=(*Grid_q)[m]){
								q_dum=(*W->Grid_int_q)[j][alpha];
									
								if(alpha==0) difq=q_dum;	
								else difq=q_dum-(*W->Grid_int_q)[j][alpha-1];			
									
								if(m==0) dif1=delta_H[n*size_q+m];
								else dif1=delta_H[n*size_q+m-1]+(delta_H[n*size_q+m]-delta_H[n*size_q+m-1])/((*Grid_q)[m]-(*Grid_q)[m-1])*(q_dum-(*Grid_q)[m-1]);

								sum-=(*GG->cgfactor_h)[n]*difq/((*Grid_f_h)[n]*(*Grid_f_h)[n]+(q_dum*q_dum-miu)*(q_dum*q_dum-miu))*Interpolate(i+size_f_l,n+size_f_l,j,alpha)*dif1;
								alpha++;
							}
						}
						
					}
				sum=sum*Temp/4/M_PI/M_PI;	
				kH.push_back(sum);
				}	      
			}
		delta_H.swap(kH);
	}
}


