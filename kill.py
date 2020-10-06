import os

size=8
j=19925
for i in range (0,size):
	myCmd='kill -9 {0}'.format(str(j))
	j=j+2
	os.system(myCmd)


