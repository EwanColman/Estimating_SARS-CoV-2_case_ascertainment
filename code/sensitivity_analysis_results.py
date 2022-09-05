# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 17:49:58 2021

@author: ewanc
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df=pd.read_csv('../output/time_independent_rates.csv')

#df2=pd.read_csv('../output/time_independent_rates_(higher_sensitivity).csv')
#df2=pd.read_csv('../output/time_independent_rates_(longer_delay).csv')
df2=pd.read_csv('../output/time_independent_rates_(all_LFD).csv')



changes=[]
for variant in ['Wild','Alpha','Delta','BA.1','BA.2']:
    group=df['Group'].tolist()
    original=df[variant+'_Rate'].tolist()
    adjusted=df2[variant+'_Rate'].tolist()
    
    for i in range(len(group)):
        print(group[i],round(adjusted[i]-original[i],2))
        changes.append(adjusted[i]-original[i])
print()   
print('Mean:',np.mean(changes))
print('SD:',np.std(changes))

#print(df.columns)

