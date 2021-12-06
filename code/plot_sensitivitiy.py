# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 16:52:00 2020

@author: ewanc
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import erf

import pandas as pd

#fig=plt.figure(figsize=(10,10))
fig = plt.figure(figsize=(10,3))
fs=12

ax = fig.add_subplot(1,2,1)


df=pd.read_csv('../raw_data/PCR_curve_summary.csv')

print(len(df))

S=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]
S=S+[0 for i in range(30)]
#[0,0.1,0.20,0.75,0.80,0.80,0.75,0.70,0.65,0.60,0.55,0.50,0.45,0.40,0.35,0.30,0.25,0.20,0.16,0.13,0.10,0.08,0.06,0.05]
#S=S+[0 for i in range(30,testable_period)]
df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]



#sensitivity=sensitivity+[0,0,0,0]

plt.plot(S, c='k',linewidth=3, label='$q=PCR$')
plt.plot(S_lfd,':', c='k',linewidth=3, label='$q=LFD$')
plt.xlabel('Time since exposure (days), $\\tau$',size=fs)
plt.ylabel('Test sensitivity, $S_{q}(\\tau)$',size=fs)
plt.yticks([0,0.4,0.8],[0,0.4,0.8],size=fs)
plt.xticks([0,10,20,30,40,50],size=fs)
plt.ylim([0,1])
plt.xlim([0,35])
plt.text(-10,1.04,'A',size=20)


#ax = fig.add_subplot(1,3,2)

#variant_sensitivity=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,120,10)]+[np.mean(df['median'].tolist()[120:130]) for i in range(120,300,10)]
#testable_period=45
#peak=0
#s_max=0.7
#variant_sensitivity=[i*s_max/peak for i in range(peak)]+[s_max*(testable_period-i)/(testable_period-peak) for i in range(peak,testable_period)]
s_max=0.4
sigma=2
pi=4
testable_period=60#int(30*stretch)
#interval=int(10/stretch)
#print(interval)
'''
curves={(0.5,1.5):'--c',(0.9,2):':m'}

for (s_max,sigma) in curves.keys():
    #Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,40,10)]
    #Sv=Sv+[np.mean(df['median'].tolist()[i:i+interval]) for i in range(40,300,interval)]
    #Sv=Sv+[0 for i in range(len(Sv),testable_period)]
    #Sv=[s*s_max/0.785 for s in Sv]
    
    Sv=S[0:pi]
    for tau in range(pi,testable_period):
        f=(tau-pi)/sigma
        y=S[pi+int(f)]+(f-int(f))*(S[pi+int(f)+1]-S[pi+int(f)])    
        Sv.append(y)
    
    plt.plot(Sv,curves[(s_max,sigma)],linewidth=3,label='$s_{max}='+str(round(100*s_max))+'$%, $\phi='+str(sigma)+'$')
'''
plt.legend(prop={'size': fs})


ax = fig.add_subplot(1,2,2)

# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean=5.5
dispersion=1.52

sigma=np.log(dispersion)
mu=np.log(mean)    

incubation_probability=[]
for i in range(len(S)):
    x=i+1
    u=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    x=i
    l=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    p=u-l
    incubation_probability.append(p) 

plt.plot(incubation_probability, c='k',linewidth=3)
plt.xlabel('Incubation period (days), $\\tau$',size=fs)
plt.ylabel('Probability, $R(\\tau)$',size=fs)
plt.yticks([0,0.1,0.2],size=fs)
plt.xticks([0,5,10,15,20,25],size=fs)
plt.ylim([0,0.21])
plt.xlim([0,25])

plt.text(-4,0.22,'B',size=20)


plt.subplots_adjust(wspace=0.4)
plt.savefig('../figures/figure1.pdf',format='pdf',dpi=300,bbox_inches='tight')


    
    

    