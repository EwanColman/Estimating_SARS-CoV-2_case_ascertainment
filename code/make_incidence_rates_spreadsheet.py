# -*- coding: utf-8 -*-
"""
Created on Thu May 13 23:05:03 2021

@author: ewanc
"""
import pickle as pk
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from datetime import datetime



dic={}
# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
print(df.head())
variant_proportion=df['NV_proportion'].tolist()
dic['Date']=df['Date'].tolist()
dic['England_alpha_proportion']=variant_proportion



#I=pk.load(open('../pickles/England_incidence.p','rb'))
#dic['incidence_England']=[int(i) for i in I]


######## AGES ##########################################################
population_of={'02_10':6264662,
               '11_15':2664513,
               '16_24':5860347,
               '25_34':7567108,
               '35_49':10864046,
               '50_69':13688509,
               '70+':8114862
               }

#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(age).p','rb'))

ages=['02_10','11_15','16_24','25_34','35_49','50_69','70+']

i=0
for age in ages:
      
    I=data['incidence_'+age]
   
    dic[age+'_incidence']=[int(i*population_of[age]/100) for i in I]
    
############# REGIONS ###################################

population_of={'NorthWest':7341196,
               'YorkshireandTheHumber':5502967,
               'NorthEast':2669941,
               'WestMidlands':5934037,
               'EastMidlands':4835928,
               'EastofEngland':6236072,
               'SouthWest':5624696,
               'London':8961989,
               'SouthEast':9180135,
               'Wales':3152879,
               'Scotland':5463300,
               'NorthernIreland':1893667
               }

#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(region).p','rb'))

regions=['SouthWest',
        'London',
        'SouthEast',
        'WestMidlands',
        'EastMidlands',
        'EastofEngland',
        'NorthWest',
        'YorkshireandTheHumber',
        'NorthEast']

nations=['Wales',
         'Scotland',
         'NorthernIreland']

for region in regions+nations:
    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    cases=df['cases'].tolist()
    time_in_days=df['Days_since_March1'].tolist()
    variant_proportion=df['NV_proportion'].tolist()
    
    I=data['incidence_'+region]
   
    dic[region+'_incidence']=[int(i*population_of[region]/100) for i in I]
    dic[region+'_alpha_proportion']=variant_proportion

    
    


length=min([len(d) for d in dic.values()])

dic['England_incidence']=[sum([dic[region+'_incidence'][t] for region in regions])for t in range(length)]



for d in dic:
    while len(dic[d])>length:
        dic[d].pop()

pd.DataFrame(dic).to_csv('../output/Incidence_estimates_and_alpha_proportions.csv',index=False)
