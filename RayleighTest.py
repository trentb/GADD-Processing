# -*- coding: utf-8 -*-
"""
Created on Tue Feb 03 17:22:40 2015

@author: Trent
"""

import diel
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import os
filename="7-31-14_An_XI-D1_100Hz-1MHz_2V.diel"

thickness=1
diameter=600
minpts=10
maxpts=32
maxi=8
freq=1000

s,m=diel.getdata("7-31-14_An_XI-D1_100Hz-1MHz_2V.diel")
E=[x*10/thickness for x in s["C AC_set"]]
E=np.asarray(E)
data=m["CMPLX CD0"]
data=[x.real*thickness*10**-6/((8.854*10**-12)*((diameter/2*10**-6)**2*np.pi))for x in data]
k=np.array(data)
k=k.reshape(len(s["A FREQUENCY_set"]),len(s["B COUNTER"]),len(s["C AC_set"]))

freqindex=s["A FREQUENCY_set"].index(freq)
#Placeholder variables:
bestslopediff=1000
bestslope=1
bestrsq=0

for i in range(0,maxi):
    for j in range(minpts,maxi+maxpts):
        if j-i < maxpts and j-i > minpts:
            slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(E[i:j],k[freqindex,1,i:j])
            slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(E[i:j],k[freqindex,2,i:j])
            if bestrsq < (r_value2**2+r_value3**2)/2:
                bestslopediff=abs(slope2-slope3)
                besti=i
                bestj=j
                bestslope=(slope2+slope3)/2
                bestintercept=(intercept2+intercept3)/2
                bestrsq=(r_value2**2+r_value3**2)/2
                bestpvalue=(p_value2+p_value3)/2
                beststderror=(std_err2+std_err3)/2

averagek=(k[freqindex,1,:]+k[freqindex,2,:])/2
Efit=np.linspace(min(E),max(E),100)
Kfit=bestintercept+bestslope*Efit

plt.plot(E,averagek, 'o', linestyle='None', linewidth=0)
plt.plot(Efit,Kfit)
plt.ticklabel_format(style='plain')
plt.tick_params(labelsize=14)
plt.ylabel('Dielectric Constant', fontsize=18)
plt.xlabel('Electric Field (kV/cm)', fontsize=18)
plt.title('Rayleigh Data', fontsize=24)

Ekdata=np.vstack((E,averagek))
savefile=os.path.splitext(filename)[0]
savefile=savefile+"_average.csv"
np.savetxt(savefile, Ekdata, delimiter=",",)
