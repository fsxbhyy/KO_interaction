#include "head.h" 
#include <iostream>
#include <math.h>
double Inner(const std::vector<double> &a,const std::vector<double> &b){
        double k=0;
        int i;
        for(i=0;i<a.size();i++){
                k=k+a[i]*b[i];
        }
        return k;
}



void Map_difw_to_w2(const std::vector<double>& wl,const std::vector<double>& wh, const std::vector<double>& wdif, std::vector<int>& list){
	int i,j;
	double w1,w2;
	for(i=0;i<wl.size()+wh.size();i++){
		for(j=0;j<wl.size()+wh.size();j++){
			if(i<wl.size()) w1=wl[i];
			else w1=wh[i-wl.size()];
			if(j<wl.size()) w2=wl[j];
			else w2=wh[j-wl.size()];
			list.push_back(Search(fabs(w1-w2),wdif));
			w2=-w2;
			list.push_back(Search(fabs(w1-w2),wdif));
		}
	}
}

int Search(const double& w,const std::vector<double>& wdif){ //binary search
	int result;
	int front,middle,end;
	front=0;
	end=wdif.size()-1;
	if(w>=wdif[end]) result=end;
	else {
		do{
			middle=(front+end)/2;
			if(w<wdif[middle]) end=middle;
			else front=middle;	
		}while((w-wdif[middle])*(w-wdif[middle+1])>0);
		if(w==wdif[middle+1]) result=middle+1;
		else result=middle;
	}
//	printf("w=%f wdif=%f\n",w,wdif[result]);
	return result;
}

