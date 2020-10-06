import os

size=8
T=0.05
for i in range (0,size):
##	if(i<size):
##		 j=0
##	else:
##		j=1
#	myCmd='cp ./MC/pi_1.h5 ./MC/pi_{0}.h5 '.format(str(i),str(j),str(i))
#	os.system(myCmd)
#	myCmd='nohup python config.py {0} {1}'.format(str(i),str(T))
#	os.system(myCmd)
	myCmd='nohup ./calc_pi.exe {0} >out{1}.txt 2>&1 &'.format(str(i),str(i))
	os.system(myCmd)
#	T=T/2**0.5
