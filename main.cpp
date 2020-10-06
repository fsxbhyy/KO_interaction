#include <iostream>
#include "head.h"
#include <functional>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <signal.h>
#include <string>
#include "H5Cpp.h"
#include <fenv.h>



/* run this program using the console pauser or add your own getch, system("pause") or input loop */
/*#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

typedef boost::mt19937 RNGType;
RNGType rng(1);
boost::uniform_01<> zero_to_one;
boost::variate_generator< RNGType, boost::uniform_01<> >
             dice(rng, zero_to_one);

double RNG(){

        return dice();
}
*/

//#include "boost/random.hpp"

/*#include <random>
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<double> dis(0, 1);

double RNG(){
        return dis(gen);
}*/


volatile int flag=0;

bool is_file_exist(const char *fileName)
{
    std::ifstream infile(fileName);
    return infile.good();
}


void signal_handler(int sig_code){
        if(sig_code==SIGINT)
                flag=1;
}

double inf_sum(const double& q,const int& n){
  double a=q*q;
  double sum=1.0;
  int i=0;
  double j=1.0;
  double k=2.0;
  for(i=0;i<n;i++){
    sum=sum+a/j/k;
    a=a*q*q;
    j=j*(i+2.0);
    k=k*(i+3.0);
  } 
  return 1.0/sum/sum;
}


double GG_cal(const double& q,const double& w){
	double result=q*q/(w*w+(q*q-1)*(q*q-1));
	return result; 
}

void Initialize(std::vector<std::vector<double> >& grid_int_q,const std::vector<double>& grid_q,const double& qc,const double& cut){
	int i,j;
	double w;
        double step;
	int grid_size=4;
	double slow_down=1.0;
	int i0,j0;
	for(i=0;i<grid_q.size();i++){
		j=0;
		step=cut;
		w=grid_q[i]-step;
		std::vector<double> dum1;
		std::vector<double> dum2;
		std::vector<double> dum3;
		while(w>0){
			dum1.push_back(w);
			if(j%grid_size==0) {
				if(step<slow_down) step=step*2;	
				else step=slow_down;
				/*	else{
					if(swit==0){
						swit=1;
						step0=step;
					}	
					step=step+step0;
				}*/
			}
			w=w-step;
			j++;
		}
		for(j=dum1.size()-1;j>=0;j--){
			dum2.push_back(dum1[j]);
		}
	
		j=0;
                step=cut;
                w=grid_q[i]+step;
		while(w<qc){
			dum2.push_back(w);
			if(j%grid_size==0) {
				if(step<slow_down) step=step*2;	
				else step=slow_down;
				/*else{
					if(swit==0){
						swit=1;
						step0=step;
					}	
					step=step+step0;
				}*/
			}
                        w=w+step;
                        j++;
		}
		dum2.push_back(w);

		i0=0;
		j0=0;
		double mark;
		while(i0<grid_q.size()&&j0<dum2.size()){
			if(grid_q[i0]<dum2[j0]) {
				if(i0!=i){
					dum3.push_back(grid_q[i0]);
				}
				i0++;
			}	
			else if(grid_q[i0]==dum2[j0]) {
				if(i0!=i) {
					dum3.push_back(grid_q[i0]);
				}
				i0++;
				j0++;
			}
			else if(grid_q[i0]>dum2[j0]){
				dum3.push_back(dum2[j0]);
				j0++;
			}
		
		}
		while(i0<grid_q.size()){
			dum3.push_back(grid_q[i0]);
			i0++;
		}
		while(j0<dum2.size()){
                        dum3.push_back(dum2[j0]);
                        j0++;
                }


		grid_int_q.push_back(dum3);
	//	grid_int_q.push_back(grid_q);
	}
}
	
        



int main(int argc,char** argv) {
	flag=0;
	int channel=0;
	double g=0;
	feenableexcept(FE_INVALID | FE_OVERFLOW);
	signal(SIGINT,signal_handler);
	std::string slabel;
        std::string name1="w0_";
	std::string name2="Gapfunction_";
	std::string name3="w1_";
	std::string name4="MC/pi_";
	std::string name5="delta_";
	std::string name6="MC/configure_";

	
	std::ofstream file;
        std::ofstream file2;
	std::ofstream file3;
	std::ifstream file4;
        if(argc==3){
               	slabel=argv[1];
		channel=atoi(argv[2]);
		g=3.126370;
        }        
	else {
		std::cin>>slabel>>channel>>g;	
	}
	name1=name1+slabel+".txt";
        name2=name2+slabel+".txt";
	name3=name3+slabel+".txt";
   	name4=name4+slabel+".h5";
	name5=name5+slabel+".h5";
	name6=name6+slabel+".txt";
	
//	dice.engine().seed(27*atoi(slabel.c_str())+time(0));
        file.open(name1.c_str());
        file2.open(name2.c_str());
	file3.open(name3.c_str());
	
	H5::H5File file_data;
	H5::DataSet dataset;
	std::vector<double> cgf_l;
	std::vector<double> cgf_h;
	std::vector<double> grid_q;
	std::vector<double> grid_f_l;
	std::vector<double> grid_f_h;
	std::vector<int> map_wdif_to_w;
	std::vector<std::vector<double> > grid_int_q;
	int i,j;
	int grid_size1=4;
	int grid_size2=8;
	double w,w2;
	double fc=15.0;
	double qc=10.0;
        double qf=1.0;
       	double Temp=0.001;
 	double gg=0;
	double omega=0.1;
	if(is_file_exist(name6.c_str())){
		file4.open(name6.c_str());
		printf("read configure");
		file4>>fc>>qf>>qc>>gg>>Temp;
		file4.close();	
	}
	else printf("no configure\n");


  fc=50.0;
	qc=10.0;
	omega=0.25;
  double r_s=1.0;
	g=r_s*1.042;
  Temp=0.0125;
	printf("fc=%f qf=%f qc=%f g=%f Temp=%f\n",fc,qf,qc,g,Temp);
	printf("omega=%f\n",omega);
	i=0;
	w=0;
	double equal_cut=1e-6;

	double A1,A2,B1,B2,C1,C2,D;
	double r_s_dl=sqrt(4*0.521*r_s/M_PI);
	C1=1-r_s_dl*r_s_dl/4.0*(1+0.07671*r_s_dl*r_s_dl*((1+12.05*r_s_dl)*(1+12.05*r_s_dl)+4.0*4.254/3.0*r_s_dl*r_s_dl*(1+7.0/8.0*12.05*r_s_dl)+1.5*1.363*r_s_dl*r_s_dl*r_s_dl*(1+8.0/9.0*12.05*r_s_dl))/
				(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl)/(1+12.05*r_s_dl+4.254*r_s_dl*r_s_dl+1.363*r_s_dl*r_s_dl*r_s_dl));
	C2=1-r_s_dl*r_s_dl/4.0*(1+r_s_dl*r_s_dl/8.0*(log(r_s_dl*r_s_dl/(r_s_dl*r_s_dl+0.990))-(1.122+1.222*r_s_dl*r_s_dl)/(1+0.533*r_s_dl*r_s_dl+0.184*r_s_dl*r_s_dl*r_s_dl*r_s_dl)));
	D=inf_sum(r_s_dl,100);
	A1=(2.0-C1-C2)/4.0/g*M_PI;
	A2=(C2-C1)/4.0/g*M_PI;
	B1=6*A1/(D+1.0);
	B2=2*A2/(1.0-D);
//read frequency-difference gridf
	file_data.openFile(name4.c_str(),H5F_ACC_RDWR);
	hsize_t dims_out[1];
	dataset=file_data.openDataSet("/pi/w");
	dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
	std::vector<double> grid_f(dims_out[0],0);
	dataset.read(&(grid_f[0]),H5::PredType::IEEE_F64LE);
	
//read momentum_dif grid_dq	
	dataset=file_data.openDataSet("/pi/q");
	dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
	std::vector<double> grid_dq(dims_out[0],0);
	dataset.read(&(grid_dq[0]),H5::PredType::IEEE_F64LE);
//read value of PI	
	dataset=file_data.openDataSet("/pi/pi");
	dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
	std::vector<double> PI_value(dims_out[0],0);
//	std::vector<double> PI_value1(dims_out[0],0);
//	std::vector<double> PI_value2(dims_out[0],0);
	dataset.read(&(PI_value[0]),H5::PredType::IEEE_F64LE);
	for(j=0;j<grid_dq.size();j++){
	  double G_s=A1*grid_dq[j]*grid_dq[j]/(1.0+B1*grid_dq[j]*grid_dq[j])+A2*grid_dq[j]*grid_dq[j]/(1.0+B2*grid_dq[j]*grid_dq[j]);
	  double G_a=A1*grid_dq[j]*grid_dq[j]/(1.0+B1*grid_dq[j]*grid_dq[j])-A2*grid_dq[j]*grid_dq[j]/(1.0+B2*grid_dq[j]*grid_dq[j]);
	  for(i=0;i<grid_f.size();i++){
	
//		printf("q=%f w=%f PI=%f\n",grid_dq[j],grid_f[i],PI_value[i*grid_dq.size()+j]);
	//		PI_value[i*grid_dq.size()+j]=-PI_value[i*grid_dq.size()+j]/(grid_dq[j]*grid_dq[j]/4/M_PI/g+PI_value[i*grid_dq.size()+j]);
		if(i>0){
		  //PI_value[i*grid_dq.size()+j]=-PI_value[i*grid_dq.size()+j]/((grid_dq[j]*grid_dq[j]+1.0)/4/M_PI/g+PI_value[i*grid_dq.size()+j]);
		  PI_value[i*grid_dq.size()+j]=-PI_value[i*grid_dq.size()+j]*(1-G_s)*(1-G_s)/(grid_dq[j]*grid_dq[j]/4/M_PI/g+PI_value[i*grid_dq.size()+j]*(1-G_s))+3*PI_value[i*grid_dq.size()+j]*G_a*G_a/(grid_dq[j]*grid_dq[j]/4/M_PI/g-PI_value[i*grid_dq.size()+j]*G_a);
		}
		else{
		  //PI_value[i*grid_dq.size()+j]=4*M_PI*g/(grid_dq[j]*grid_dq[j]+1.0+4*M_PI*g*PI_value[i*grid_dq.size()+j]);
		  PI_value[i*grid_dq.size()+j]=4*M_PI*g/grid_dq[j]/grid_dq[j]*G_s+4*M_PI*g*(1-G_s)/(grid_dq[j]*grid_dq[j]+4*M_PI*g*PI_value[i*grid_dq.size()+j]*(1-G_s))+3*PI_value[i*grid_dq.size()+j]*(4*M_PI*g/grid_dq[j]/grid_dq[j]*G_a)*(4*M_PI*g/grid_dq[j]/grid_dq[j]*G_a)/(1.0-PI_value[i*grid_dq.size()+j]*G_a*4*M_PI*g/grid_dq[j]/grid_dq[j]);
		}
	//		PI_value[i*grid_dq.size()+j]=4*M_PI*g/grid_dq[j]/grid_dq[j]/(1+g/(grid_dq[j]*grid_dq[j]+grid_f[i]*grid_f[i]));
	//		if(i==0) printf("%f %f %f\n",PI_value[i*grid_dq.size()+j],grid_dq[j],grid_f[i]);
	//		PI_value1[i*grid_dq.size()+j]=PI_value[i*grid_dq.size()+j]*grid_dq[j]*grid_dq[j];
	//		PI_value2[i*grid_dq.size()+j]=PI_value1[i*grid_dq.size()+j]*grid_dq[j]*grid_dq[j];

		}
	}
	file_data.close();
        std::vector<double> delta_L_swap;
        std::vector<double> delta_H_swap;

//Manually add grid_f_dif and grid_q_dif
/*	i=0;
	w=0;
	std::vector<double> grid_f;
	while(w<2*fc){ //0->wf
		grid_f.push_back(w);
		w=w+2*M_PI*Temp*pow(2.0,i/grid_size2);
		i++;
	}
	grid_f.push_back(w);

	int grid_size_dq=20;
	i=0;
	w=1e-8;
	double step_dq=w;
	std::vector<double> grid_dq;
	while(w<2*qc){ 
		grid_dq.push_back(w);
	//	w=w+1e-6;
		w=w+step_dq*pow(2.0,i/grid_size_dq);
		i++;
	}
	printf("lastq=%f\n",w);
	grid_dq.push_back(w);

	std::vector<double> PI_value;
	for(i=0;i<grid_f.size();i++){
		for(j=0;j<grid_dq.size();j++){
			PI_value.push_back(4*M_PI*g/grid_dq[j]/(1+g/(grid_dq[j]*grid_dq[j]+grid_f[i]*grid_f[i])));
		//	PI_value.push_back(grid_dq[j]);
		}
	}*/

//read frequency-difference gridf
if(is_file_exist(name5.c_str())){	
	printf("read old figure\n");
	file_data.openFile(name5.c_str(),H5F_ACC_RDWR);
	dataset=file_data.openDataSet("/delta/w");
        dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
        std::vector<double> grid_ff(dims_out[0],0);
        dataset.read(&(grid_ff[0]),H5::PredType::IEEE_F64LE);
	for(i=0;i<grid_ff.size();i++){
		if(grid_ff[i]<omega){
			grid_f_l.push_back(grid_ff[i]);
			if(i==grid_ff.size()-1) cgf_l.push_back(round((grid_ff[i]-grid_ff[i-1])/2/M_PI/Temp));
			else cgf_l.push_back(round((grid_ff[i+1]-grid_ff[i])/2/M_PI/Temp));
		}
		else{
			grid_f_h.push_back(grid_ff[i]);
			if(i==grid_ff.size()-1) cgf_h.push_back(round((grid_ff[i]-grid_ff[i-1])/2/M_PI/Temp));
                        else cgf_h.push_back(round((grid_ff[i+1]-grid_ff[i])/2/M_PI/Temp));

		}
	}
	
//read momentum_dif grid_q     
        dataset=file_data.openDataSet("/delta/q");
        dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
        for(i=0;i<dims_out[0];i++){
                grid_q.push_back(0);
        }
        dataset.read(&(grid_q[0]),H5::PredType::IEEE_F64LE);
	
//read delta L and delta H
	dataset=file_data.openDataSet("/delta/delta_L");
        dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
        for(i=0;i<dims_out[0];i++){
                delta_L_swap.push_back(0);
        }

	dataset.read(&(delta_L_swap[0]),H5::PredType::IEEE_F64LE);
/*	for(i=0;i<grid_q.size();i++){
        	for(j=0;j<grid_f_l.size();j++){
			if(grid_q[i]<1.0) delta_L_swap[j*grid_q.size()+i]/=grid_q[i]*grid_q[i];
        	}
	}*/

//read delta L and delta H
	dataset=file_data.openDataSet("/delta/delta_H");
        dataset.getSpace().getSimpleExtentDims(dims_out,NULL);
        for(i=0;i<dims_out[0];i++){
                delta_H_swap.push_back(0);
        }
//        dataset.read(&(delta_H_swap[0]),H5::PredType::IEEE_F64LE);
/*	for(i=0;i<grid_q.size();i++){
                for(j=0;grid_f_h.size();j++){
                        delta_H_swap[j*grid_q.size()+i]=delta_H_swap[j*grid_q.size()+i]/10000.0;
                }
        }*/

	file_data.close();
}

	//Establish the frequency grids f_l and f_h;
else{
	i=0;
        w=M_PI*Temp;
        double step_f=1.0;
        double f_init=1.0;
        double f_middle=10.0;
        int N_highf1=40;
        int N_highf2=20;
        while(w<fc){ //0->wf
                if(w<omega){
                        grid_f_l.push_back(w);
                        cgf_l.push_back(step_f);
                }
                else{
                        grid_f_h.push_back(w);
                        cgf_h.push_back(step_f);
                }
                w=w+2*M_PI*Temp*step_f;
                if(w<f_init) {
                        i++;
                        step_f=pow(2,i/grid_size2);

                }
                else if(w<f_middle) step_f=floor((f_middle-f_init)/2/M_PI/Temp/N_highf1);
               // else if(w<fc) step_f=floor((fc-f_middle)/2/M_PI/Temp/N_highf2);
		else if(w<fc) {
			i++;
			if(i%16==0) step_f=step_f*2;
		}

        }





// Establish the momentum grid q;
	i=0;
        w=1e-5;
        double step_q=w;
	int counter;
        while(w<0.3){ //0->qf
                grid_q.push_back(w);
                w=w+step_q;
                if(i%grid_size1==0) { 
			step_q=step_q*3;	
		}
		i++;
        }
	counter=grid_q.size()-1;
	while(w<0.7){
		grid_q.push_back(w);
		w=w+0.02;
	}
        for(i=counter;i>-1;i--){
                w=qf-grid_q[i];
                grid_q.push_back(w);
        }

        i=0;
	w=1.0;
	step_q=1e-5;
        while(w<qc){ //qf->qc
                grid_q.push_back(w);
                w=w+step_q;
		if(i%grid_size1==0) step_q=step_q*3;
		if(step_q>0.5) step_q=0.5;
        	i++;
	}
	for(i=0;i<grid_q.size()*grid_f_l.size();i++){
                        delta_L_swap.push_back(1.0);
               }
	for(i=0;i<grid_q.size()*grid_f_h.size();i++){
                        delta_H_swap.push_back(-0.1);
               }


	
}

	for(i=0;i<grid_q.size();i++){
		if(i==0) file<<"q\n";
		else file<<grid_q[i]<<"\t"<<grid_q[i]-grid_q[i-1]<<"\n";
	}
	for(i=0;i<grid_dq.size();i++){
                if(i==0) file<<"qdif\n";
                file<<grid_dq[i]<<"\n";
        }
	for(i=0;i<grid_f.size();i++){
                if(i==0) file<<"wdif\n";
                file<<grid_f[i]<<"\n";
        }
	for(i=0;i<grid_f_l.size();i++){
                if(i==0) file<<"wlow\n";
		else file<<grid_f_l[i]<<"\t"<<(grid_f_l[i]-grid_f_l[i-1])/2/M_PI/Temp<<"\t"<<cgf_l[i]<<"\n";
        }
	for(i=0;i<grid_f_h.size();i++){
                if(i==0) file<<"whigh\n";
		else file<<grid_f_h[i]<<"\t"<<(grid_f_h[i]-grid_f_h[i-1])/2/M_PI/Temp<<"\t"<<cgf_h[i]<<"\n";
        }
	double qcc=qc;
	Initialize(grid_int_q,grid_q,qcc,equal_cut);

	Map_difw_to_w2(grid_f_l,grid_f_h,grid_f,map_wdif_to_w);
	PI_function PI(PI_value,&grid_dq,&grid_f,equal_cut,g);
	W_function W(&PI,&grid_q,&grid_f,&grid_int_q,grid_size2,channel);
	
	G_function GG(&grid_q,&grid_f_l,&grid_f_h,&cgf_l,&cgf_h);
	Gap delta(delta_L_swap,delta_H_swap,&grid_q,&grid_f_l,&grid_f_h,&map_wdif_to_w,&W,&GG,Temp);	
/*
for(NF=0;NF<grid_f_l.size()+grid_f_h.size();NF++){
for(wf=0;wf<grid_f_l.size()+grid_f_h.size();wf++){
	if(wf<grid_f_l.size()) freq=grid_f_l[wf];
	else freq=grid_f_h[wf-grid_f_l.size()];
	for(i=0;i<grid_q.size();i++){
                sum=0;
		step=1e-12;
		distance=step;              
	      	for(j=0;distance<equal_cut;j++){
			aa1=GG_cal(grid_q[i]+distance,freq)*delta.Interpolate_singular(PI,NF,wf,grid_q[i],distance);
			aa2=GG_cal(grid_q[i]-distance,freq)*delta.Interpolate_singular(PI,NF,wf,grid_q[i],-distance);	
			//	sum=sum+step*(GG_cal(grid_q[i]+distance,freq)*delta.Interpolate_singular(PI,NF,wf,grid_q[i],distance)+GG_cal(grid_q[i]-distance,freq)*delta.Interpolate_singular(PI,NF,wf,grid_q[i],-distance));
                  	sum=sum+step*(aa1+aa2);
		  	if(j%grid_size0==0) step=step*2;
                        distance+=step;
                }
                delta.sing.push_back(sum);
//	printf("q=%f w=%f aa1=%f aa2=%f\n",grid_q[i],grid_f[NF],aa1,aa2);
        }
}
}*/
	
	printf("size_q=%d size_f_l=%d size_f_h=%d\n",int(grid_q.size()),int(grid_f_l.size()),int(grid_f_h.size()));
 	
	int i1;
/*//	i=Search(46.2254,grid_f);
//	i=0;
	double w_test;
	int max=2;
	printf("%f\n",grid_f[i]);
	for(i=0;i<max;i=i+max-1){	
		for(j=0;j<grid_q.size();j++){
			for(i1=0;i1<grid_int_q[j].size();i1++){
				w=W.value[j][i1*grid_f.size()+i];
				w_test=PI.Analy_limit(i,grid_q[j],grid_int_q[j][i1],channel);
//				if(w_test!=0) w=fabs(w-w_test)/w;
				//if(j==0&&grid_int_q[j][i1]<1e-4) 
			if(i==1)	file<<grid_q[j]<<"\t"<<grid_int_q[j][i1]<<"\t"<<w<<"\t"<<w_test<<"\n";
//			else if(i==max-1)	file3<<grid_q[j]<<"\t"<<grid_int_q[j][i1]<<"\t"<<w<<"\t"<<w_test<<"\n";
			//	if(fabs((w-w_test)/w_test)<0.1) file3<<grid_q[j]<<"\t"<<grid_int_q[j][i1]<<"\t"<<w<<"\n";
			}
		
		file<<"\n";		
		}
	}*/
/*	for(j=0;j<grid_dq.size();j++){
               	for(i1=0;i1<grid_f.size();i1++){
                                w=PI.value[i1*grid_dq.size()+j];
                              //  if(w!=0) w=log(w);
                                file3<<grid_dq[j]<<"\t"<<grid_f[i1]<<"\t"<<w<<"\n";
                        }
                        file3<<"\n";    
        }*/


//Check the map from f_dif to f;	
/*	for(i=0;i<grid_f_l.size()+grid_f_h.size();i++){
		for(j=0;j<grid_f_l.size()+grid_f_h.size();j++){
			if(i<grid_f_l.size()) w=grid_f_l[i];
			else w=grid_f_h[i-grid_f_l.size()];
			if(j<grid_f_l.size()) w2=grid_f_l[j];
			else w2=grid_f_h[j-grid_f_l.size()];
			printf("w1-w2=%f w1+w2=%f map1=%f map2=%f\n",w-w2,w+w2,grid_f[map_wdif_to_w[2*i*(grid_f_l.size()+grid_f_h.size())+2*j]],grid_f[map_wdif_to_w[2*i*(grid_f_l.size()+grid_f_h.size())+2*j+1]]);		
		}
	}*/

/*	int n,m;
	double H,H0;
	int size_00=grid_f_l.size()+grid_f_h.size()-1;
	for(i=0;i<grid_f_l.size()+grid_f_h.size();i++){
		for(j=0;j<grid_f_l.size()+grid_f_h.size();j++){                                      
                    	if(i<grid_f_l.size()) w=grid_f_l[i];
			else w=grid_f_h[i-grid_f_l.size()];
			if(j<grid_f_l.size()) w2=grid_f_l[j];
                        else w2=grid_f_h[j-grid_f_l.size()];
		     	for(n=0;n<grid_q.size();n++){
                               for(m=0;m<grid_int_q[n].size();m++){
				        
					H0=delta.Interpolate(i,j,n,m);
		//			if(i==grid_f_l.size()+grid_f_h.size()-1&&j==grid_f_l.size()+grid_f_h.size()-1) file3<<w-w2<<"\t"<<w+w2<<"\t"<<grid_q[n]<<"\t"<<grid_int_q[n][m]<<"\t"<<H0<<"\t"<<H<<"\n";
					if(i==0&&j==1) file3<<w-w2<<"\t"<<w+w2<<"\t"<<grid_q[n]<<"\t"<<grid_int_q[n][m]<<"\t"<<H0<<"\n";
				
				}
			}

		}
	}*/

/*	for(i=0;i<grid_f.size();i++){
		     	for(n=0;n<grid_q.size();n++){
                               for(m=0;m<grid_q.size();m++){
				        H=H1_analy(grid_f[i],grid_q[n],grid_q[m]);
					if(n==m) file3<<grid_f[i]<<"\t"<<grid_q[n]<<"\t"<<grid_q[m]<<"\t"<<W.value[i*grid_q.size()*grid_q.size()+n*grid_q.size()+m]<<"\t"<<H<<"\n";
				}
			}
		}*/
	double lamu_dum=1;						
	double lamu=0;
	i=0;
  j=0;
	do{
		if(flag==0){
			lamu_dum=lamu;
		//	printf("%f %f\n",delta.inhomoL[0],delta.delta_L[0]);
//		printf("%f %f %f %f\n",delta.inhomoL[0],delta.inhomoH[0],delta.delta_L[0],delta.delta_H[0]);
			delta.MOperate(2);	//calculate K21*Gap1
		//		printf("%f %f\n",delta.inhomoL[0],delta.delta_L[0]);
//						printf("%f %f %f %f\n",delta.inhomoL[0],delta.inhomoH[0],delta.delta_L[0],delta.delta_H[0]);
			delta.DampIteration(j);
		//	printf("%f %f\n",delta.inhomoL[0],delta.delta_L[0]);
//						printf("%f %f %f %f\n",delta.inhomoL[0],delta.inhomoH[0],delta.delta_L[0],delta.delta_H[0]);
			delta.MOperate(1);	
			//	printf("%f %f\n",delta.inhomoL[0],delta.delta_L[0]);
//						printf("%f %f %f %f\n",delta.inhomoL[0],delta.inhomoH[0],delta.delta_L[0],delta.delta_H[0]);
			lamu=delta.Power();	//calculate K12*Gap2
			//	printf("%f %f\n",delta.inhomoL[0],delta.delta_L[0]);
//						printf("%f %f %f %f\n",delta.inhomoL[0],delta.inhomoH[0],delta.delta_L[0],delta.delta_H[0]);
			printf("%f\n",lamu);
			i++;
      if(fabs(lamu-lamu_dum)<1e-4){
        j=1;
      }
      // if(fabs(lamu-lamu_dum)<1e-4){
      //   j=2;
      // }
		}
		else i=1000;
	}while(fabs(lamu-lamu_dum)>1e-5&&i<100);
	
	for(i1=0;i1<grid_q.size();i1++){
		for(j=0;j<grid_f_l.size();j++){
			file2<<grid_q[i1]<<"\t"<<grid_f_l[j]<<"\t"<<delta.delta_L[j*grid_q.size()+i1]<<"\n";
		}
		for(j=0;j<grid_f_h.size();j++){
                	file2<<grid_q[i1]<<"\t"<<grid_f_h[j]<<"\t"<<delta.delta_H[j*grid_q.size()+i1]<<"\n";
		}
	}
	file<<lamu;
/*	if(is_file_exist(name5.c_str())){ 
        if(remove( name5.c_str() ) != 0)
         	perror( "Error deleting delta file" );
         	else 
        	puts( "delta file successfully deleted" );
	}


	H5::H5File file_Gap;

	file_Gap=H5::H5File(name5.c_str(),H5F_ACC_TRUNC);
	H5::Group g2(file_Gap.createGroup("/delta"));
 	hsize_t dims_in[1];
    //store freq grid
    	std::vector<double> combine_l_h;
	for(i=0;i<delta.Grid_f_l->size();i++){
		combine_l_h.push_back((*delta.Grid_f_l)[i]);
	}
	for(i=0;i<delta.Grid_f_h->size();i++){
                combine_l_h.push_back((*delta.Grid_f_h)[i]);
        }
	dims_in[0]=combine_l_h.size();
    	H5::DataSpace dataspace=H5::DataSpace(1,dims_in);
    	dataset=H5::DataSet(file_Gap.createDataSet("/delta/w",H5::PredType::IEEE_F64LE,dataspace));
	dataset.write(&(combine_l_h[0]),H5::PredType::IEEE_F64LE);
	dataspace.close();
 	dataset.close();


    //store val grid
   	dims_in[0]=delta.Grid_q->size();
	dataspace=H5::DataSpace(1,dims_in);
    	dataset=H5::DataSet(file_Gap.createDataSet("/delta/q",H5::PredType::IEEE_F64LE,dataspace));
    	dataset.write(&((*delta.Grid_q)[0]),H5::PredType::IEEE_F64LE);
    	dataspace.close();
    	dataset.close();
	
	dims_in[0]=delta.delta_L.size();
        dataspace=H5::DataSpace(1,dims_in);
        dataset=H5::DataSet(file_Gap.createDataSet("/delta/delta_L",H5::PredType::IEEE_F64LE,dataspace));
        dataset.write(&(delta.delta_L[0]),H5::PredType::IEEE_F64LE);
        dataspace.close();
        dataset.close();
	
	dims_in[0]=delta.delta_H.size();
        dataspace=H5::DataSpace(1,dims_in);
        dataset=H5::DataSet(file_Gap.createDataSet("/delta/delta_H",H5::PredType::IEEE_F64LE,dataspace));
        dataset.write(&(delta.delta_H[0]),H5::PredType::IEEE_F64LE);
        dataspace.close();
        dataset.close();
*/
		
	file.close();
	file2.close();
	file3.close();
	
}





