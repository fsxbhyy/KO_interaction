#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>
PI_function::PI_function(std::vector<double> &a,std::vector<double> *b,std::vector<double> *c,const double& cut,const double& g):Grid_dq(b),Grid_f(c),equ_cut(cut),g_pi(g){
	this->value.swap(a);
};

double PI_function::Integrate(std::vector<double> &val,const int& h,const double& q1,const double& q2){//h is the label for frequency
        int i;
	int j;
	int size=Grid_dq->size();
        double sum=0;
	double sum1=0;
	double sum_total=0;
	double difq=fabs(q1-q2);
//	if(qdif_cutoff<(*Grid_dq)[0]) {
//		throw "Need to have a larger qdif_cutoff!";
//	}
	if(difq==0) difq=equ_cut;
	for(i=0;(*Grid_dq)[i]<difq;i++){
		sum+=(val[h*size+i]+val[h*size+i+1])/2.0*((*Grid_dq)[i+1]-(*Grid_dq)[i]);
	}
	//Integrate 0 to |q1-q2|
	if(i==0) sum=-val[h*size+i]*((*Grid_dq)[i]-difq);
	else sum=sum-(val[h*size+i]+(val[h*size+i-1]+(val[h*size+i]-val[h*size+i-1])*(difq-(*Grid_dq)[i-1])/((*Grid_dq)[i]-(*Grid_dq)[i-1])))/2.0*((*Grid_dq)[i]-difq);
	
	sum1=sum;
	j=i;

	//Integrate 0 to q1+q2
	sum=0;
	for(i=0;(*Grid_dq)[i]<(q1+q2)&&i<size-1;i++){
		sum+=(val[h*size+i]+val[h*size+i+1])/2.0*((*Grid_dq)[i+1]-(*Grid_dq)[i]);
	}
	if(i==0) sum=val[h*size+i]*((*Grid_dq)[i]-difq);
	else sum=sum-(val[h*size+i]+(val[h*size+i-1]+(val[h*size+i]-val[h*size+i-1])*(q1+q2-(*Grid_dq)[i-1])/((*Grid_dq)[i]-(*Grid_dq)[i-1])))/2.0*((*Grid_dq)[i]-q1-q2);
//	if((sum-sum1)/q1/q2>100) printf("q1=%f q2=%f i1=%d i2=%d sum=%f sum1=%f final=%f\n",q1,q2,i,j,sum,sum1,(sum-sum1)/q1/q2);

	sum_total=(sum-sum1);		
	return sum_total;
}


double PI_function::Integrate2(const int& h,const double& q1,const double& q2,const int& channel){//h is the label for frequency
        int i;
	int j;
	double chi0=-1.0;
        double qdum=fabs(q1-q2);	
	std::vector<double> chi_grid;
	double grid_size=100.0;
	double cut=0.0005;
	if(qdum<cut){
		double step=1e-8;
		std::vector<double> chi_grid_dum;
		i=0;
		chi0=1.0;
//		if(fabs(q1-q2)<5e-6&&q1+q2<3e-5) step=1e-9;
		while(qdum<q1+q2&&qdum<cut){
			if(chi0<=1.0) chi_grid_dum.push_back(chi0);
			qdum=qdum+step;
			i++;
			if(i%300==0) step=step*2.0;
			chi0=(q1*q1+q2*q2-qdum*qdum)/2/q1/q2;
		}
		if(q1+q2>cut){
		//	int chi_num=chi_grid_dum.size()-1;
		//	chi0=chi_grid_dum[chi_num];
			chi0=chi_grid_dum[chi_grid_dum.size()-1];
			step=(chi0+1.0)/grid_size;
			chi0=chi0-step;
			while(chi0>-1.0){
				chi_grid_dum.push_back(chi0);
		                chi0=chi0-step;
			}
		}
		if(fabs(chi_grid_dum[chi_grid_dum.size()-1]+1)>1e-7){
			chi_grid.push_back(-1.0);
		}
		for(i=chi_grid_dum.size()-1;i>-1;i--){
			chi_grid.push_back(chi_grid_dum[i]);
		}
	/*	if(q1+q2>cut){
			for(i=1;i<chi_grid.size();i++){
				if(chi_grid[i]-chi_grid[i-1]<0){
					printf("%f %f\n",q1,q2); 
					printf("%f %f\n",chi_grid[i],chi_grid[i]-chi_grid[i-1]);
				}
			}
		//	printf("*");
		}*/
	}
	else{
		for(i=0;i<grid_size;i++){
			chi_grid.push_back(chi0);
			chi0=chi0+2.0/grid_size;	
		}
	
		chi_grid.push_back(1.0);
	}
	int size=Grid_dq->size();
        double sum=0;
	double difq=fabs(q1-q2);
	double sumq=q1+q2;
	double delta=chi_grid[1]-chi_grid[0];
	
	double middle=0;
	double middle_chi=0;
	double chi_dum=0;
	double chi_old=0;
	qdum=0;
//	if(qdif_cutoff<(*Grid_dq)[0]) {
//		throw "Need to have a larger qdif_cutoff!";
//	}
	if(difq==0) difq=equ_cut;
	
	
	for(i=size-1;(*Grid_dq)[i]>sumq;i--){
	}
	
	

/*	for(;(*Grid_dq)[i]>difq;i--){
		chi=(q1*q1+q2*q2-(*Grid_dq)[i]*(*Grid_dq)[i])/2/q1/q2;
	//	printf("%f\n",chi);
		if((chi-old_chi)<0) throw 1;
		if(swit==1){ 
			middle=value[h*size+i]+(value[h*size+i+1]-value[h*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])*(sumq-(*Grid_dq)[i])/2.0;
			middle_chi=((*Grid_dq)[i]+sumq)/2.0;
			middle_chi=(q1*q1+q2*q2-middle_chi*middle_chi)/2.0/q1/q2;
			swit=0;
		}
		else {
			middle=(value[h*size+i]+value[h*size+i+1])/2.0;
			middle_chi=(chi+old_chi)/2.0;
		}
//		middle=value[h*size+i];
		if(i==0) printf("%f\n",middle);
		delta=delta+(chi-old_chi);
		if(channel==0) sum+=(chi-old_chi)*middle;	
		else if(channel==1) sum+=(chi-old_chi)*middle*middle_chi;
		else if(channel==2) sum+=(chi-old_chi)*middle*(3*middle_chi*middle_chi-1)/2.0;
		old_chi=chi;
	}
	
	chi=1.0;
	middle=value[h*size+i]+(value[h*size+i+1]-value[h*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])*(((*Grid_dq)[i+1]+difq)/2-(*Grid_dq)[i]);
       	middle_chi=(difq+(*Grid_dq)[i+1])/2.0;	
   	if(channel==0) sum+=(chi-old_chi)*middle;
       	else if(channel==1) sum+=(chi-old_chi)*middle*middle_chi;
	else if(channel==2) sum+=(chi-old_chi)*middle*(3*middle_chi*middle_chi-1)/2.0;*/
	int test=0;
	double sum_delta=0;

	for(j=0;j<chi_grid.size()-1;j++){
		chi_dum=(q1*q1+q2*q2-(*Grid_dq)[i]*(*Grid_dq)[i])/2/q1/q2;
		if(chi_dum<chi_grid[j+1]){
			int swit=1;
			chi_old=chi_grid[j];
			while(swit==1){
				if(chi_dum>=chi_grid[j+1]){
				       	chi_dum=chi_grid[j+1];
					swit=0;
				}
				delta=chi_dum-chi_old;
				sum_delta+=delta;
				if(delta<0) printf("%f\n",delta);
				middle_chi=(chi_dum+chi_old)/2.0;
//				if(fabs(middle_chi)>1.0) std::cout<<middle_chi<<std::endl;
				qdum=sqrt(fabs(q1*q1+q2*q2-2*q1*q2*middle_chi));
				if(qdum==0) qdum=equ_cut;
				middle=value[h*size+i]+(value[h*size+i+1]-value[h*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])*(qdum-(*Grid_dq)[i]);
				if(h>0) middle=middle*4*M_PI*g_pi/qdum/qdum;
				if(channel==0) sum+=delta*middle;
        		        else if(channel==1) sum+=delta*middle*middle_chi;
		                else if(channel==2) sum+=delta*middle*(3*middle_chi*middle_chi-1)/2.0;
			
				//update chi_dum and chi_old
				chi_old=chi_dum;
				if(swit==1) i--;
				chi_dum=(q1*q1+q2*q2-(*Grid_dq)[i]*(*Grid_dq)[i])/2/q1/q2; 
				if(i<0) {
					swit=0;
					i=0;
				}
			}
		}
		
		else{
			qdum=sqrt(q1*q1+q2*q2-2*q1*q2*(chi_grid[j]+chi_grid[j+1])/2.0);
			if(qdum==0) qdum=equ_cut;
			middle=value[h*size+i]+(value[h*size+i+1]-value[h*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])*(qdum-(*Grid_dq)[i]);
			if(h>0) middle=middle*4*M_PI*g_pi/qdum/qdum;
			middle_chi=(chi_grid[j]+chi_grid[j+1])/2.0;
			delta=chi_grid[j+1]-chi_grid[j];
			sum_delta+=delta;
			if(channel==0) sum+=delta*middle;
                	else if(channel==1) sum+=delta*middle*middle_chi;
                	else if(channel==2) sum+=delta*middle*(3*middle_chi*middle_chi-1)/2.0;	
		}
		
	}	
//	printf("%f\n",sum_delta);
	
//	for(i=0;i<N;i++){
			
//	}	
//	printf("%f\n",delta);
	return sum;
}
double PI_function::Analy_limit(const int& w,const double& k1,const double& q1,const int& channel){
	double x=sqrt(k1*k1+q1*q1);
	double y=k1*q1;
	double result=0;
	double W_f;
	double slope;
	int size=Grid_dq->size();
	int i=Search(x,(*Grid_dq));
	if(i<size) W_f=value[w*size+i]+(value[w*size+i+1]-value[w*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])*(x-(*Grid_dq)[i]);
	else W_f=value[w*size+i];
	if(channel==0) {
		W_f=W_f*4*M_PI*g_pi/x/x;
		result=2.0*y*W_f;
	}
	else if(channel==1){
		if(i>0&&i<size) slope=((value[w*size+i+1]-value[w*size+i])/((*Grid_dq)[i+1]-(*Grid_dq)[i])+(value[w*size+i]-value[w*size+i-1])/((*Grid_dq)[i]-(*Grid_dq)[i-1]))/2.0;
		else slope=0;
		result=-(2.0/3.0)*y*y*4*M_PI*g_pi/x/x/x*(slope-2.0*W_f/x);	
	}
	
	return result;
}

double H1_analy0(const double& equ_cut,const double& k1,const double& q1,const double& g,const int& channel){
//      double Omega=0.2*M_PI;
        int i;
	double result;
        long double dum1,dum2,dum3;
        dum2=k1-q1;
        if(fabs(dum2)==0) dum2=equ_cut;
        dum1=4*k1*q1/dum2/dum2;
	dum3=(k1+q1)*(k1+q1)/(dum2*dum2);
	if(channel==0) {
		if(dum1<1e-5) result=2*M_PI*g*log1p(dum1); 
		else result=2*M_PI*g*log(dum3);
	}
	else if(channel==1) result=4*M_PI*g*k1*q1;
	else if(channel==2) result=4*M_PI*g*(2*k1*q1*(k1*k1+q1*q1));
//    	result=result*k1*q1;
    	return result;
}


