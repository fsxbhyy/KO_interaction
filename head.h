#ifndef CIRCLE_H
#define CIRCLE_H
#include <vector>
#define M_PI 3.1415926535897932384

class PI_function{
	public:
		std::vector<double> value;
		std::vector<double> value1;
		std::vector<double> value2;
		std::vector<double> *Grid_dq;
		std::vector<double> *Grid_f;
		double equ_cut;
		double g_pi;
		PI_function(std::vector<double>&,std::vector<double> *,std::vector<double> *,const double&,const double&);			
		double Integrate(std::vector<double>&,const int&,const double&,const double&);
		double Integrate2(const int& h,const double& q1,const double& q2,const int& channel);
		double Analy_limit(const int& w,const double& k1,const double& q1,const int& channel);
};

class W_function{
	public:
		std::vector<std::vector<double> > value;
		std::vector<double> *Grid_q;
		std::vector<double> *Grid_f;		
		std::vector<std::vector<double> > *Grid_int_q;
		W_function(PI_function*,std::vector<double>*,std::vector<double>*,std::vector<std::vector<double> >*,const int&,const int&);	
		PI_function* raw;
                int grid_size;
		int chan;
};

class G_function{
	public:
		std::vector<double> value_L;
		std::vector<double> value_H;	
		std::vector<double> *Grid_q;
		std::vector<double> *Grid_f_l;
		std::vector<double> *Grid_f_h;
		std::vector<double> *cgfactor_l;
                std::vector<double> *cgfactor_h;

		G_function(std::vector<double>*,std::vector<double>*,std::vector<double>*,std::vector<double>*,std::vector<double>*);
};



class Gap{
	public:
		std::vector<double> inhomoH;
		std::vector<double> inhomoL;
		std::vector<double> *Grid_q;
		std::vector<double> *Grid_f_l;
		std::vector<double> *Grid_f_h;
		std::vector<int> *Map; 
//		std::vector<double> sing;
		W_function* W;
	        G_function* GG;	
		
		double Temp;		
		double shift;	
		
		std::vector<double> delta_L;
		std::vector<double> delta_H;
		double miu;
		double H_accum;	
		double Interpolate(const int&,const int&,const int&,const int&);		
//		double Interpolate_singular(PI_function&,const int&,const int&,const double&,const double&);		
		Gap(std::vector<double> &,std::vector<double> &,std::vector<double> *,std::vector<double> *gridfl,std::vector<double> *gridfh,std::vector<int>*,W_function* w,G_function* gg,const double& T);		
		void DampIteration(const int&);
		double Power();
		void MOperate(const int&);
	
};
double Inner(const std::vector<double>&,const std::vector<double>&);
double RNG();
long double H1_analy(const double& w1,const double& w2,const double& k1,const double& q1);
double H1_analy0(const double& w,const double& k1,const double& q1,const double&, const int&);
void Map_difw_to_w2(const std::vector<double>&,const std::vector<double>&,const std::vector<double>&,std::vector<int>&);
int Search(const double&,const std::vector<double>&);
void signal_handler(int sig_code);
#endif
