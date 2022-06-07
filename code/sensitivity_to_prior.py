# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:52:00 2020

@author: ewanc
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

import pandas as pd


delta=1

### Lognormal incubation period distribution
# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean,dispersion=5.5,1.52
mu,sigma=np.log(mean),np.log(dispersion)
testable_period=30

R=[]
for i in range(testable_period):
    # daily probability is the difference of cumulative distributions at consecutive time points
    x=i+1
    u=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    x=i
    if x>0:
        l=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    else:
        l=0
    p=u-l
    R.append(p)

#### Probability of positive test #################
# the value in position i is th probability of testing positive i days after infection
#from https://www.medrxiv.org/content/10.1101/2020.11.24.20229948v1.full.pdf
df=pd.read_csv('../raw_data/PCR_curve_summary.csv')
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]

df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]

#S=S+[0 for i in range(30,testable_period)]

g1=0.85
g2=1.15

colour={g1:'r',1:'g',g2:'b'}
PCR_mean={}
LFD_mean={}

for g in colour:
    print('growth rate '+str(g)+',' )
    # make the probability function as a list
    P_pcr=[]
    for j in range(delta,testable_period):
        numerator=R[j-delta]*S_pcr[j]*(g**(-j))
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([R[i-delta]*S_pcr[i]*(g**(-i)) for i in range(testable_period)])
        P_pcr.append(numerator/denominator) 
    
        
    # make the probability function as a list
    P_lfd=[]
    for j in range(testable_period):
        # use -j (and -i) since we are lookig backwards in time
        numerator=R[j]*S_lfd[j]*(g**(-j))
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([R[i]*S_lfd[i]*(g**(-i)) for i in range(testable_period)])
        P_lfd.append(numerator/denominator) 
    
    

        lab1='$P_{PCR}(\\tau)$, with daily growth rate ='+str(g)
        lab2='$P_{LFD}(\\tau)$, with daily growth rate ='+str(g)
    PCR_mean[g]=sum(i*P_pcr[i] for i in range(len(P_pcr)))
    LFD_mean[g]=sum(i*P_lfd[i] for i in range(len(P_lfd)))
    print('PCR mean='+str(PCR_mean[g]))
    print('LFD mean='+str(LFD_mean[g]))
    plt.plot(P_pcr,c=colour[g], label=lab1) 
    plt.plot(P_lfd,c=colour[g],linestyle=':',label=lab2) 

plt.legend()
plt.xlabel('Days since exposure, $\\tau$')

print()
print('PCR')
print('g='+str(g1)+', change of '+str(round(PCR_mean[1]-PCR_mean[g1],2))+' days')
print('g='+str(g2)+', change of '+str(round(PCR_mean[1]-PCR_mean[g2],2))+' days')
print()
print('LFD')
print('g='+str(g1)+', change of '+str(round(LFD_mean[1]-LFD_mean[g1],2))+' days')
print('g='+str(g2)+', change of '+str(round(LFD_mean[1]-LFD_mean[g2],2))+' days')
  