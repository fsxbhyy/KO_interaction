#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>
W_function::W_function(PI_function* c,std::vector<double>* a,std::vector<double>* b,std::vector<std::vector<double> >* e,const int& gsize,const int& channel):Grid_q(a),Grid_f(b),Grid_int_q(e),raw(c),grid_size(gsize),chan(channel){	
	int i,j,k,i1;
	int size_f=Grid_f->size();
	double q1=0;
	double q2=0;
	double H=0;
	double H1=0;
	double H2=0;
	double HF=0;
	double w;
	double g_ph=1.0;
	double omega_ph=0.05;
	for(j=0;j<Grid_q->size();j++){
		std::vector<double> dum;
		for(k=0;k<(*Grid_int_q)[j].size();k++){
			q1=(*Grid_q)[j];
			q2=(*Grid_int_q)[j][k];
			for(i=0;i<size_f;i++){
				H=raw->Integrate2(i,q1,q2,channel);
				//H=H+2*g_ph*(*Grid_f)[i]*(*Grid_f)[i]/((*Grid_f)[i]*(*Grid_f)[i]+omega_ph*omega_ph);
				H=H*q1*q2;
				if(i>0) H+=H1_analy0(raw->equ_cut,q1,q2,raw->g_pi,channel);
				dum.push_back(H);
				

			/*	H=raw->Integrate(raw->value,i,q1,q2);
                                if(channel==0) {
                                        HF=H/q1/q2;
                                }
                                else if(channel==1) {
                                        H1=raw->Integrate(raw->value1,i,q1,q2);
                                        HF=(H*(q1*q1+q2*q2)-H1)/2/q1/q1/q2/q2;
                //                      HF=H*(q1*q1+q2*q2)/q1/q1/q2/q2/2;
                        //              HF=-H1/2/q1/q1/q2/q2;
                                }
                                else if(channel==2) {
                                        H1=raw->Integrate(raw->value1,i,q1,q2);
                                        H2=raw->Integrate(raw->value2,i,q1,q2);
                                        HF=(H*(3*(q1*q1+q2*q2)*(q1*q1+q2*q2)-4*q1*q1*q2*q2)-6*H1*(q1*q1+q2*q2)+3*H2)/8/q1/q1/q1/q2/q2/q2;
                                }
                                dum.push_back(HF);*/
		//		else{
		//			q=value[i*size_q*size_q+k*size_q+j]; //inverse j and k, since symmetric matrix;
	//	this->value.push_back(q);
	//			}
			//	H=H1_analy((*Grid_f)[i],(*Grid_q)[j],(*Grid_q)[k]);
			//	this->value.push_back(H);
			}
		}
		this->value.push_back(dum);
	}
	
}




long double H1_analy(const double& w1,const double& w2,const double& k1,const double& q1){
//	double Omega=0.2*M_PI;
	int i;
	double w;
	double H=0;
	double dum1,dum2,dum3;
	double g=3.0;
	
	dum2=k1-q1;        
        if(fabs(dum2)==0) dum2=1e-8;
	dum1=4*k1*q1/dum2/dum2;
	w=w1-w2;
for(i=0;i<1;i++){
        dum3=4*k1*q1/(dum2*dum2+w*w+g);
	H+=2*M_PI*g*(log1p(dum1)-g/(w*w+g)*(log1p(dum1)-log1p(dum3)));
	w=w1+w2;
}
	
//	H=1.0*w1*w1/(w1*w1+Omega*Omega);
	return H;
}


