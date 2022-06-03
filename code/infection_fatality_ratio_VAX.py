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
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

lag=21
name_of={#'02_10':'2 to 10',
         #'11_15':'11 to 15',
         #'16_24':'16 to 24',
         '25_34':'25 to 34',
         '35_49':'35 to 49',
         '50_69':'50 to 69',
         '70+':'70+'}
max_day=800
fs=12
start_date='1 August 2020'
end_date='15 April 2022'

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
start_date=(datetime.strptime(start_date, '%d %B %Y')-time_zero).days
end_date=(datetime.strptime(end_date, '%d %B %Y')-time_zero).days
###########################################
months=26

# get the range of days corresponding to each month
day_range={}

first_day_of_month=0
for month in range(1,months+1):
    # find what month it is
    m=1+(2+month)%12
    if m>9:
         date='01/'+str(m)+'/202'+str(int((month+2)/12))
    else:
        date='01/0'+str(m)+'/202'+str(int((month+2)/12))
    print(month,date) 
    # find what the first day of the month is in days since march 1st 2020
    days_since_March1=(datetime.strptime(str(date), '%d/%m/%Y')-time_zero).days
    #print(days_since_March1)

    day_range[month-1]=(first_day_of_month-lag,days_since_March1-lag)
    first_day_of_month=days_since_March1

for month in range(months-1):
    print(month,day_range[month])



# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]\
    +['01/0'+str(i)+'/2021' for i in range(1,10)]+['01/'+str(i)+'/2021' for i in range(10,13)]\
    +['01/0'+str(i)+'/2022' for i in range(1,4)]    
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
  


fig = plt.figure(figsize=(12,8))
gs = fig.add_gridspec(2,2)

#fig=plt.figure(figsize=(8,4))
#plt.subplots_adjust(hspace=0,wspace=0)


######## AGES ##########################################################


data=pk.load(open('../pickles/reporting_rates_(age).p','rb'))
vax_data=pk.load(open('../pickles/reporting_rates_(age)_VAX.p','rb'))



population_of={'02_10':6254603,
               '11_15':3370248,
               '16_24':5950637,
               '25_34':7596145,
               '35_49':10853151,
               '50_69':13618246,
               '70+':7679719
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
         '35_49':0.15,
         '50_69':1,
         '70+':11}

i=0
for age in name_of:
    
#    n=i % 3
#    m=1+int(n/3)
    
    #ax=fig.add_subplot(3,3,i) 
    ax = fig.add_subplot(gs[int(i/2),i%2])
    i=i+1
    
    #ax.axvline((datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,linewidth=1)
    ax.set_xlim([start_date,end_date])
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
    ####### add dashboard deaths ############
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
    
    monthly_deaths2=[None]+[sum(deaths[day_range[m][0]+lag:day_range[m][1]+lag]) for m in range(1,months)]
    #####################
    
        
    if age=='35_49':
        leg1='Original IFR'
        leg1_vax='With vaccine effect'
        leg2='Plausible region'
    else:
        leg1=None
        leg1_vax=None
        leg2=None

    I_dic=data['incidence_'+age]
    print(I_dic)
    I=I_dic['Rate']
    
    # aggregate to months 
    monthly_infections=[population_of[age]*sum(I[day_range[month][0]:day_range[month][1]])/100 for month in range(months-1)]
    #ONS
    #IFR=[100*monthly_deaths[age][i]/monthly_infections[i] for i in range(5,months-1)]
    #28days
    IFR=[100*monthly_deaths2[i]/monthly_infections[i] for i in range(5,months-1)]
    
    
    
    
    plt.plot([day_range[m][1]+lag for m in range(5,months-1)],IFR,marker='o',markersize=2,color='k',linewidth=1,zorder=1,label=leg1)
    #plt.fill_between([day_range[m][1]+lag for m in range(5,months-1)],IFR['Lower'],IFR['Upper'],color='c',alpha=0.2,linewidth=1)
   
    I_dic=vax_data['incidence_'+age]
    print(I_dic)
    I=I_dic
    
    # aggregate to months 
    monthly_infections=[population_of[age]*sum(I[day_range[month][0]:day_range[month][1]])/100 for month in range(months-1)]
    #ONS
    #IFR_vax=[100*monthly_deaths[age][i]/monthly_infections[i] for i in range(5,months-1)]
    #28days
    IFR_vax=[100*monthly_deaths2[i]/monthly_infections[i] for i in range(5,months-1)]
    
    
    plt.plot([day_range[m][1]+lag for m in range(5,months-1)],IFR_vax,color='k',linewidth=0.5,zorder=2,label=leg1_vax)
    
    ax.fill_between([day_range[m][1]+lag for m in range(5,months-1)],IFR,IFR_vax,color='c',alpha=0.2,linewidth=1,label=leg2)
    
    if age=='35_49':
        plt.legend(loc=1,prop={'size':fs})#,frameon=False,bbox_to_anchor=(0.5, 0.5))



plt.savefig('../figures/Supplementary_figure_IFR.pdf',format='pdf',dpi=300,bbox_inches='tight')

