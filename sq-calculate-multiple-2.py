#!/usr/bin/python
import math,sys,os,subprocess
import string               
import numpy as np
import multiprocessing

##calculate the strucutre factor in alt-bulk system

##definite calculate function
def sq_cal(frame):
    Sq_frame = []
    Sq_q = {}
    Nq_q = {}
    for nx in range(nxmin,nxmax):
        for ny in range(nymin,nymax):
            for nz in range(nzmax):
                if nx == 0 and ny == 0 and nz == 0:
                   continue   
                qx, qy, qz = dqx*nx, dqy*ny, dqz*nz
                q = math.sqrt(pow(qx,2) + pow(qy,2) + pow(qz,2))
                if q >= qmax:
                   continue
                cossum, sinsum = 0.0, 0.0
                for p in points[frame]:
                    cossum = cossum + math.cos(qx*p[1]+qy*p[2]+qz*p[3])
                    sinsum = sinsum + math.sin(qx*p[1]+qy*p[2]+qz*p[3])
                Sq = pow(cossum,2) + pow(sinsum,2)

                if Sq_q.has_key(q) == True:
                    Sq_q[q] = Sq_q[q]+Sq
                    Nq_q[q] += 1
                else:
                    Sq_q[q] = Sq
                    Nq_q[q] = 1
    c = 1.0/float(part_num[frame])
    for key in sorted(Sq_q):
        Sq_frame.append(float(Sq_q[key]*c/Nq_q[key]))
    return Sq_frame

##read files and definite varibles
nframe = 10                
f = open('strfact.xyz','r')
fout = open('sq.dat','w')

#atom = str(raw_input("the atom type: "))
atom = 'O'
bin_factor = float(raw_input("the bin_factor: "))

fgro = open('md.gro','r')
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
    
## calculate the value of q
q_all = []
for nx in range(nxmin,nxmax):
    for ny in range(nymin,nymax):
        for nz in range(nzmax):
            if nx == 0 and ny == 0 and nz == 0:
                continue   
            qx, qy, qz = dqx*nx, dqy*ny, dqz*nz
            q = math.sqrt(pow(qx,2) + pow(qy,2) + pow(qz,2))
            if q >= qmax:
                continue
            if (q in q_all) == True:
                continue
            else:
                q_all.append(q)
q_all.sort()

##multiprocess to calculate Sq for each frame
p = multiprocessing.Pool(processes=10)
frames = [0,1,2,3,4,5,6,7,8,9]
Sq_all = p.map(sq_cal, frames)
p.close()
p.join()

##calculate q_max and average Sq and std
Sq_max, q_max = 0,0
Sq_bin, q_bin, std_bin = [0]*bin_num, [0]*bin_num, [0]*bin_num
for i in range(len(q_all)):
    Sq_i = []
    for j in range(nframe):
        Sq_i.append(Sq_all[j][i])
    Sq_aver = np.mean(Sq_i)
    if Sq_max < Sq_aver:
        Sq_max = Sq_aver
        q_max = q_all[i]
    bin_i = int(q_all[i]/bin_size)
    if Sq_bin[bin_i] < Sq_aver:
        Sq_bin[bin_i] = Sq_aver
        q_bin[bin_i] = q_all[i]
        std_bin[bin_i] = np.std(Sq_i, ddof=1)

for bin_i in range(bin_num):
    if q_bin[bin_i] != 0:
        print >> fout, q_bin[bin_i], Sq_bin[bin_i], std_bin[bin_i]

print "q_max = ", q_max, "     AND d = ", float(2*math.pi/q_max)
