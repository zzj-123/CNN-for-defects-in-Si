# -*- coding: utf-8 -*-
"""
Created on Mon Sep  7 20:24:39 2020

@author: 34595
"""

import numpy as np
import os
def discal(c,t):
    dis=round(np.linalg.norm(c-t),4)
    O=[]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                l=np.array([i,j,k])
                if (l==np.array([1,1,1])).all()==False:
                    o=l-np.array([1,1,1])
                    O.append(o)
    O=np.array(O)
    O=O*16.29
    for i in range(26):
        nt=t+O[i]
        if round(np.linalg.norm(nt-c),4)<dis:
            dis=round(np.linalg.norm(nt-c),4)
    return dis


def rdfcal(atomlist,r_cutoff,dr):
    n=atomlist.shape[0]
    grl=np.linspace(0,r_cutoff,round(r_cutoff/dr)+1)
    gr=np.zeros(round(r_cutoff/dr)+1)
    for i in range(n):
        for j in range(i+1,n):
            dis=discal(atomlist[i],atomlist[j])
            if dis<r_cutoff:
                gr[int(dis/dr)]+=1
    for i in range(1,gr.shape[0]):
        gr[i]=gr[i]/n
        
        gr[i]=gr[i]/(4*np.pi*dr*((dr*i)**2))
    return grl,gr
for kk in range(1,10):
    f1=open("../6.071/"+str(kk)+"/dump2.atom")
    data=f1.readlines()
    n=int(data[3])
    La=round(float(data[5].split()[1]))
    Lb=round(float(data[6].split()[1]))
    Lc=round(float(data[7].split()[1]))
    for i in range(0,6):
        coor=[]
        for j in range(0,n):
            # print(data[i*(9+n)+9+j])
            t1=[float(data[i*(9+n)+9+j].split()[2]),float(data[i*(9+n)+9+j].split()[3]),float(data[i*(9+n)+9+j].split()[4])]
            t=np.array(t1)
            coor.append(t)
        atomlist=np.array(coor)
        atomlist=atomlist*16.29
        grl,gr=rdfcal(atomlist,8.2,0.02)
        grho=216/(16.29**3)
        gr=gr/grho
        np.savetxt("../DNV6/rdf/6.071/"+str((kk-1)*5+i)+".txt",gr)
        np.save("../DNV6/rdf/6.071/"+str((kk-1)*5+i),gr)
        break
    break
    f1.close()