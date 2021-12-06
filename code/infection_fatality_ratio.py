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


max_day=610
fs=12
start_date='1 August 2020'
end_date='10 October 2021'

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
start_date=(datetime.strptime(start_date, '%d %B %Y')-time_zero).days
end_date=(datetime.strptime(end_date, '%d %B %Y')-time_zero).days
###########################################

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,10)]+['01/'+str(i)+'/2021' for i in range(10,11)]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
  

fig = plt.figure(figsize=(12,8))
gs = fig.add_gridspec(2,2)

#fig=plt.figure(figsize=(8,4))
#plt.subplots_adjust(hspace=0,wspace=0)


######## AGES ##########################################################


#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(age).p','rb'))

name_of={#'02_10':'2 to 10',
         #'11_15':'11 to 15',
         #'16_24':'16 to 24',
         '25_34':'25 to 34',
         '35_49':'35 to 49',
         '50_69':'50 to 69',
         '70+':'70+'}

population_of={'02_10':6264662,#7527576,
               '11_15':2664513,
               '16_24':5860347,
               '25_34':7567108,
               '35_49':10864046,
               '50_69':13688509,
               '70+':8114862
               }

initials_of={'02_10':'2-10',
         '11_15':'11-15',
         '16_24':'16-24',
         '25_34':'25-34',
         '35_49':'35-49',
         '50_69':'50-69',
         '70+':'70+'}

date=data['date']
start=0
while date[start]<start_date:
    start=start+1

end=len(date)-1
while date[end]>end_date:
    end=end-1
    
# adust the date range
date=date[start:end]
    
top={ '25_34':0.05,
         '35_49':0.1,
         '50_69':1,
         '70+':7}

i=0
for age in name_of:
    
#    n=i % 3
#    m=1+int(n/3)
    
    #ax=fig.add_subplot(3,3,i) 
    ax = fig.add_subplot(gs[int(i/2),i%2])
    i=i+1
    ax.set_xlim([start_date,end_date+10])
    #ax.axvline((datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,linewidth=1)
    
    ax.set_ylim([0,top[age]])
    ax.text(start_date+10,0.9*top[age],'Age '+name_of[age],size=fs)
    


    if i>2:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i in [1,3]: 
        plt.ylabel('IFR (%)',size=fs)
        #plt.text(100,200,'A',size=20)
          
    ##### DEATHS ##########
    df=pd.read_csv('../raw_data/Death/'+age+'_deaths.csv')
    
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
    deaths=df['deaths'].tolist()
    death_time_in_days=df['days_since_march1'].tolist()
    dates=df['date'].tolist()


    I=data['incidence_'+age][0:end_date]
    IFR=[100*100*sum(deaths[i:i+14])/(14*max(I[i]*population_of[age],1)) for i in range(end_date)]
    plt.plot(IFR)
  
    leg2='Infection fatality ratio'
        
  

# fake plot
#ax.fill_between([0,0],[0,0],[0,0],color='m',linewidth=0,edgecolor='k',label=leg2)
ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))



plt.savefig('../figures/IFR_fig.png',format='png',dpi=300,bbox_inches='tight')

