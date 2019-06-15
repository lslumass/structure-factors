#!/usr/bin/python
import math,sys,os,subprocess
import string               
import numpy as np
import multiprocessing

##calculate the strucutre factor in alt-bulk system

##definite calculate function for each frame
def sq_cal(frame):
    Sq_frame = [0]*bin_num
    q_bin_frame, Nq_frame = [0]*bin_num, [0]*bin_num
    for nx in range(nxmin,nxmax):
        for ny in range(nymin,nymax):
            for nz in range(nzmax):
                if nx == 0 and ny == 0 and nz == 0:
                   continue   
                qx, qy, qz = dqx*nx, dqy*ny, dqz*nz
                q = math.sqrt(pow(qx,2) + pow(qy,2) + pow(qz,2))
                if q >= qmax:
                   continue
                bin_i = int(q/bin_size)
                Nq_frame[bin_i] += 1
                
                cossum, sinsum = 0.0, 0.0
                for p in points[frame]:
                    cossum = cossum + math.cos(qx*p[1]+qy*p[2]+qz*p[3])
                    sinsum = sinsum + math.sin(qx*p[1]+qy*p[2]+qz*p[3])
                Sq_frame[bin_i] = Sq_frame[bin_i] + pow(cossum,2)+pow(sinsum,2)

    for bin_i in range(0,bin_num):
        if Nq_frame[bin_i] == 0:
           c = 0
        else:
           c = 1.0/float(Nq_frame[bin_i])
        Sq_frame[bin_i] = Sq_frame[bin_i]*c/part_num[frame]
    return Sq_frame

##read files and definite varibles
nframe = 10                
f = open('strfact.xyz','r')          #input file strfact.xyz
fout = open('sq.dat','w')            #output file sq.dat

#atom = str(raw_input("the atom type: "))
atom = 'O'            #以体系中O粒子来计算Sq。也可以指定别的类型。
bin_factor = float(raw_input("the bin_factor: "))         #bin_factor用来控制q的bin大小。

fgro = open('md.gro','r')        #输入md.gro读取盒子大小。
lines = fgro.readlines()
N_tot = int(lines[1])
line = lines[-1]
box = line.split()
Lx, Ly, Lz = float(box[0])*10, float(box[1])*10, float(box[2])*10
dqx, dqy, dqz = float(2*math.pi/Lx), float(2*math.pi/Ly), float(2*math.pi/Lz)
dq = min(dqx,dqy,dqz) 
nxmax, nymax, nzmax = 30,30,30
nxmin, nymin = (-1)*nxmax, (-1)*nymax
qxmax, qymax, qzmax = nxmax*dqx, nymax*dqy, nzmax*dqz
qmax = min(qxmax,qymax,qzmax)
q_cut = 1.5
bin_size = bin_factor * dq
bin_num = int(math.ceil(qmax/bin_size))
fgro.close()

## read all the points
points = []
part_num = [0]*nframe
for i in range(nframe):
    m = (N_tot+2)*i+2
    n = (N_tot+2)*(i+1)
    f = open('strfact.xyz','r')
    point = []
    for lines in f.readlines()[m:n]:
        line = lines.split()
        line[0] = str(line[0])
        line[1] = float(line[1])
        line[2] = float(line[2])
        line[3] = float(line[3])
        if line[0] == atom:
           point.append(line)
           part_num[i] += 1
    points.append(point)
    f.close()

## calculate q_bin and Nq
q_bin, Nq = [0]*bin_num, [0]*bin_num
for nx in range(nxmin,nxmax):
    for ny in range(nymin,nymax):
        for nz in range(nzmax):
            if nx == 0 and ny == 0 and nz == 0:
               continue   
            qx, qy, qz = dqx*nx, dqy*ny, dqz*nz
            q = math.sqrt(pow(qx,2) + pow(qy,2) + pow(qz,2))
            if q >= qmax:
               continue
            bin_i = int(q/bin_size)
            q_bin[bin_i] += q
            Nq[bin_i] += 1
for i in range(bin_num):
    if Nq[i] == 0:
       continue
    else:
       c = 1.0/Nq[i]
    q_bin[i] = q_bin[i]*c
    
##multiprocess to calculate Sq for each frame
p = multiprocessing.Pool(processes=10)
frames = [0,1,2,3,4,5,6,7,8,9]
Sq_all = p.map(sq_cal, frames)
p.close()
p.join()

Sq = []
for bin_i in range(bin_num):
    Sq_i = []
    for i in range(nframe):
        Sq_i.append(Sq_all[i][bin_i])
    Sq.append(Sq_i)

Sq_max, q_max = 0,0
for bin_i in range(bin_num):
    if Nq[bin_i] == 0:
       continue
    Sq_aver = np.mean(Sq[bin_i])
    Sq_std = np.std(Sq[bin_i], ddof=1)
    print >> fout, q_bin[bin_i], Sq_aver, Sq_std
    if Sq_max < Sq_aver:
       Sq_max = Sq_aver
       q_max = q_bin[bin_i]

print "q_max = ", q_max, "AND d = ", float(2*math.pi/q_max)
