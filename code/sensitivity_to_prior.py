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
    x=i+1
    u=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    x=i
    l=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
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
g=1.3

# make the probability function as a list
P_pcr=[]
for j in range(delta,testable_period):
    numerator=R[j-delta]*S_pcr[j]*(g**j)
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i-delta]*S_pcr[i]*(g**i) for i in range(testable_period)])
    P_pcr.append(numerator/denominator)   

# make the probability function as a list
P_lfd=[]
for j in range(testable_period):
    numerator=R[j]*S_lfd[j]*(g**j)
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i]*S_lfd[i]*(g**j) for i in range(testable_period)])
    P_lfd.append(numerator/denominator) 


plt.plot(P_pcr) 
   
    

    