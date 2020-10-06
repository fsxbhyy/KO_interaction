import os
size=8

for i in range (0,size):
	if(i<size):
		 j=0
	else:
		j=1
	myCmd='nohup ./main {0} {1} >out{2}.txt 2>&1 &'.format(str(i),str(j),str(i))
	os.system(myCmd)

