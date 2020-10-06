# this is for generating configurations of calc_pi
# w_exp and calc_delta

import numpy as np
import sys
import math
T=1.0
mu=1.0
m=0.5
#rs=1.0/1.0421235224
rs=3
label=0
Kc=15.0
wc=40.0


kc=0.000000001

if len(sys.argv)==3:
    label=sys.argv[1]
    T=float(sys.argv[2])
qc=0.000000001;
# calculate grids
w=[0,] # freq grid of pi and w
count=0
multi=1
multi_2=1
NN=40
grid_size=16
step1=math.floor((wc*2-5*mu)/(2*np.pi*T)/NN)

while w[-1]<2*wc:    
    w.append(w[-1]+2*np.pi*T)


# while w[-1]<2*wc:
#     if w[-1]<1*mu:
#         if count==grid_size-1:
# 		    multi=multi*2
#         w.append(w[-1]+2*np.pi*T*multi)
#     elif w[-1]<5*mu:
#         w.append(w[-1]+((0.1*mu)//(2*np.pi*T)+1)*2*np.pi*T)
#     else:
#         if count==grid_size-1:
# 		    multi_2=multi_2*2
#         w.append(w[-1]+multi_2*((0.5*mu)//(2*np.pi*T)+1)*2*np.pi*T)
#     count+=1
#     if count==grid_size:
#         count=0

v=[np.pi*T,] #freq grid of delta
count=0
multi=1
while v[-1]<wc:
    if v[-1]<0.1*mu:
        v.append(v[-1]+2*np.pi*T*multi)
    elif v[-1]<0.2*mu: 
        v.append(v[-1]+((0.05*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    elif v[-1]<0.4*mu:
        v.append(v[-1]+((0.02*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    elif v[-1]<2*mu:
        v.append(v[-1]+((0.1*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    else:
        v.append(v[-1]+((0.8*mu)//(2*np.pi*T)+1)*2*np.pi*T)
    count+=1
    if count==5:
        count=0
        multi*=2
# while v[-1]<wc:
#     v.append(v[-1]+2*np.pi*T)


EF=2.0 
q=[qc,]
step=1e-6
flat=0.00029

q.append(flat)
i=0
while q[-1]+step<2.0*mu:
    q.append(q[-1]+step)
    i=i+1
    if(i%16==0):
        if(step<0.005):
	   step=step*2
	else:
	   step=0.005
	



i=0
while q[-1]<2*Kc:
    q.append(q[-1]+step)
    i=i+1
    if(i%4==0):
	step=step*2
    

k=[kc,]
Nlog=20
while k[-1]<Kc:
    if k[-1]<0.1*mu:
        step=k[-1]
        for i in range(7):
            k.append(k[-1]+step)
    elif k[-1]<0.7*mu:
        k.append(k[-1]+0.1*mu)
    elif (mu-k[-1])>kc:
        n=7
        step=(mu-k[-1])/n
        for j in range(n-1):
            k.append(k[-1]+step)
    elif k[-1]<mu:
        k.append(mu)
        k.append(mu+kc)
    elif k[-1]<1.3*mu:
        step=k[-1]-mu
        for i in range(7):
            k.append(k[-1]+step)
    elif k[-1]<4.0*mu:
        k.append(k[-1]+0.2*mu)
    else:
        k.append(k[-1]+1.0*mu)

# write pi.in, containing configs for calc_pi

name1="pi_"+label+".in"
name2="w_"+label+".in"
name3="delta"+label+".in"
file=open(name1,"w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set freq=grid\n")
for i in w:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.write("set mmt=grid\n")
for i in q:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()

# write w.in
file=open(name2,"w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set rs=%.16f\n" %rs)
file.write("set mmt=grid\n")
for i in k:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()

# write delta.in
file=open(name3,"w")
file.write("set T=%.16f\n" %T)
file.write("set mu=%.16f\n" %mu)
file.write("set m=%.16f\n" %m)
file.write("set rs=%.16f\n" %rs)
file.write("set v=grid\n")
for i in v:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.write("set mmt=grid\n")
for i in k:
    file.write("%.16f\n" %i)
file.write("end grid\n")
file.close()
