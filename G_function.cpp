#include "head.h"
#include <functional>
#include <stdio.h>
#include <math.h>
#include <iostream>

G_function::G_function(std::vector<double>* a,std::vector<double>* b,std::vector<double>* c,std::vector<double>* d,std::vector<double>* e):Grid_q(a),Grid_f_l(b),Grid_f_h(c),cgfactor_l(d),cgfactor_h(e){
	int i,j;
	double w,q;
	for(i=0;i<Grid_f_l->size();i++){
		for(j=0;j<Grid_q->size();j++){
			w=(*Grid_f_l)[i];
			q=(*Grid_q)[j];
			q=q*q-1;
			value_L.push_back(1.0/(w*w+q*q));
	//		value_L.push_back(1.0/fabs(w));
		}
	}
	for(i=0;i<Grid_f_h->size();i++){
                for(j=0;j<Grid_q->size();j++){
                        w=(*Grid_f_h)[i];
                        q=(*Grid_q)[j];
                        q=q*q-1;
                        value_H.push_back(1.0/(w*w+q*q));
		//	value_H.push_back(1.0/fabs(w));
                }
        }

}
