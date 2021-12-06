# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 10:49:03 2021

@author: ewanc
"""

import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
age_categories=['02_10','11_15','16_24','25_34','35_49','50_69']#,'70+']

groups=['00_04','05_09','10_14','15_19','20_24',
'25_29',
'30_34',
'35_39',
'40_44',
'45_49',
'50_54',
'55_59',
'60_64',
'65_69']

last_groups=['70_74',
'75_79',
'80_84',
'85_89','90+']

population_of={'02_10':6264662,#7527576,
               '11_15':2664513,
               '16_24':5860347,
               '25_34':7567108,
               '35_49':10864046,
               '50_69':13688509,
               '70+':8114862
               }
###### Vaccinations #############

## vax by age ###

url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=vaccinationsAgeDemographics&format=csv'

df=pd.read_csv(url)
# add a new column with header days since Mar1
time_in_days=[]
date=df['date'].tolist()
while date:
    d=date.pop(0)
    #print(d)
    day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
    time_in_days.append(day_numerical)

df['days_since_march1']=time_in_days
df=df.sort_values('days_since_march1',ascending=True)
df=df[df['days_since_march1']>=0]


dates=list(set(df['date'].tolist()))

vaccinated={}
for age in age_categories:
    vaccinated[age]=[0 for i in range(len(dates))]
vaccinated['70+']=[0 for i in range(len(dates))]


ages_of={'02_10':[],
         '11_15':[],
         '16_24':['16_17','18_24'],
         '25_34':['25_29','30_34'],
         '35_49':['35_39','40_44','45_49'],
         '50_69':['50_54','55_59','60_64','65_69'],
         '70+':['70_74','75_79','80_84','85_89','90+']}

plt.figure()

for group in vaccinated:
    for age in ages_of[group]:
        age_df=df[df['age']==age]    
        new_vaccinations=age_df['newPeopleVaccinatedFirstDoseByVaccinationDate'].tolist()
        
        vaccinated[group]=[vaccinated[group][i]+new_vaccinations[i] for i in range(len(dates))]


    # add 0s from march 1
    vaxed=[0 for i in range(min(df['days_since_march1']))]
    vaxed=vaxed+vaccinated[group]
    # calculate as a proportion of pop
    vax_pop=[sum(vaxed[:i])/population_of[group] for i in range(len(vaxed))]
    
    plt.plot(vax_pop,label=group)
    

plt.legend()

plt.figure()
###### Vaccinations #############
for age in population_of:
     #### vaccination ####   
    df=pd.read_csv('../raw_data/Vaccination/'+age+'_vaccinations.csv')
    # add a new column with header days since Mar1
    time_in_days=[]
    date=df['date'].tolist()
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        time_in_days.append(day_numerical)
    
    df['days_since_march1']=time_in_days
    df=df.sort_values('days_since_march1',ascending=True)
    df=df[df['days_since_march1']>=0]
    
    # add 0s from march 1
    vaxed=[0 for i in range(min(df['days_since_march1']))]
    vaxed=vaxed+df['vaccinated'].tolist()
    # calculate as a proportion of pop
    vax_pop=[sum(vaxed[:i])/population_of[age] for i in range(len(vaxed))]
    print(vax_pop)
    plt.plot(vax_pop,label=group)