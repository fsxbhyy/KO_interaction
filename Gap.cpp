#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>

Gap::Gap(std::vector<double> &delta_L_swap, std::vector<double> &delta_H_swap,std::vector<double> *griq,std::vector<double> *gridfl,std::vector<double> *gridfh,std::vector<int> *map_w,W_function* w,G_function* gg,const double& T):Grid_q(griq),Grid_f_l(gridfl),Grid_f_h(gridfh),Map(map_w),W(w),GG(gg),Temp(T),miu(1.0){
		int i;
    this->H_accum=0.0;
		this->shift=5.0;
		this->delta_L.swap(delta_L_swap);
		this->delta_H.swap(delta_H_swap);
		
		for(i=0;i<Grid_q->size()*Grid_f_l->size();i++){
//			delta_L.push_back(1.0);
			inhomoL.push_back(0);
		}
		for(i=0;i<Grid_q->size()*Grid_f_h->size();i++){
//			delta_H.push_back(1.0);
			inhomoH.push_back(0);
		}
}

double Gap::Interpolate(const int& k1,const int& k2,const int& n,const int&m){ //k1,k2 range from 0 to size_l+size_h-1
        double result=0;
        double size_q=Grid_q->size();
        double size_f_l=Grid_f_l->size();
        double size_f=2*(Grid_f_l->size()+Grid_f_h->size());
	double size_f_w=W->Grid_f->size();
        double w,w1,w2;
	int k;
	int k0=2*k2;
        int i; 
	if(k1<size_f_l) w1=(*Grid_f_l)[k1];
        else w1=(*Grid_f_h)[k1-size_f_l];
	if(k2<size_f_l) w2=(*Grid_f_l)[k2];
        else w2=(*Grid_f_h)[k2-size_f_l];
	
	double q1=(*W->Grid_q)[n];
	double q2=(*W->Grid_int_q)[n][m];

for(i=0;i<2;i++){//positive and negative w2
	k=(*Map)[k1*size_f+k0];
	w=fabs(w1-w2);
        if(k==W->Grid_f->size()-1) result+=W->value[n][m*size_f_w+k];
        else result+=W->value[n][m*size_f_w+k]+(W->value[n][m*size_f_w+k+1]-W->value[n][m*size_f_w+k])*(w-(*W->Grid_f)[k])/((*W->Grid_f)[k+1]-(*W->Grid_f)[k]);
//      if(k1==1&&k2==0) printf("%f %e %e %e %e\n",w,(*W->Grid_q)[n]-(*W->Grid_int_q)[n][m],W->value[n][m*size_f_w+k+1],W->value[n][m*size_f_w+k],result);
//	if(k1==k2&&n==0&&m==0) printf("%e %e %e\n",W->value[n][m*size_f_w+k+1],W->value[n][m*size_f_w+k],result);
	w2=-w2;	
	k0++;

//	if(fabs(w)>1e-8) result+=H1_analy0(W->raw->equ_cut,q1,q2,W->raw->g_pi,W->chan);
//	if(i==0) result=0;
}
	return result;
}

