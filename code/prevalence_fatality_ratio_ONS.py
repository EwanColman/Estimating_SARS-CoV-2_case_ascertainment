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

lag=20
name_of={#'02_10':'2 to 10',
         #'11_15':'11 to 15',
         #'16_24':'16 to 24',
         '25_34':'25 to 34',
         '35_49':'35 to 49',
         '50_69':'50 to 69',
         '70+':'70+'}

##### Process the ONS data #################
df=pd.read_excel('../raw_data/Death/ONS_Deaths.xlsx',
                 sheet_name='Table 1',
                 skiprows=4)



print(df.head())

df=df[df['Cause']=='Due to COVID-19']
# count the number of months

months=len(df)
print(df.columns)

monthly_deaths={}
for age in name_of:
    
    if age=='70+':
        column_names=[i for i in range(70,105)]+['105 and over']
    else:
        column_names=[i for i in range(int(age[0:2]),1+int(age[3:5]))]
    
    # store the list 0 index corresponds to march 2020
    # index= months since march 2020
    monthly_deaths[age]=df[column_names].sum(axis=1).tolist()
    print(age)
    print(monthly_deaths[age][-1])
    print()
##############################################

# get the range of days corresponding to each month
day_range={}

first_day_of_month=0
for month in range(1,months):
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



population_of={'02_10':6254603
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
         '35_49':0.1,
         '50_69':1,
         '70+':8}

i=0
for age in name_of:
    
#    n=i % 3
#    m=1+int(n/3)
    
    #ax=fig.add_subplot(3,3,i) 
    ax = fig.add_subplot(gs[int(i/2),i%2])
    i=i+1
    
    #ax.axvline((datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,linewidth=1)
    ax.set_xlim([start_date,end_date+10])
    #ax.set_ylim([0,top[age]])
    #ax.text(start_date+10,0.9*top[age],'Age '+name_of[age],size=fs)
    


    if i>2:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i in [1,3]: 
        plt.ylabel('IFR (%)',size=fs)
        #plt.text(100,200,'A',size=20)
          
   #### SURVEILLANCE ########
    df=pd.read_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',sep=',')
   
    
    # add a new column with header days since Mar1
    time_in_days=[]
    date=df['Date'].tolist()
    while date:
        d=date.pop(0)
        #day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
        time_in_days.append(day_numerical)
    
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1',ascending=True)
    #df.drop(df.tail(1).index,inplace=True)
    df=df[df['days_since_march1']>100] 
    # convert data frame to lists for plotting
    rate=df.reset_index()
    ONS_day=df['days_since_march1'].tolist()
    
    
    # aggregate to months 
    monthly_prevalence=[population_of[age]*np.mean([rate['Rate'][i] for i in range(len(rate)) if ONS_day[i]>=day_range[month][0] and ONS_day[i]<day_range[month][1]]) for month in range(months-1)]
        
    PFR=[monthly_deaths[age][i]/monthly_prevalence[i] for i in range(5,months-1)]
    plt.plot([day_range[m][1]+lag for m in range(5,months-1)],PFR,'-o')
  
    leg2='Infection fatality ratio'
        
  

# fake plot
#ax.fill_between([0,0],[0,0],[0,0],color='m',linewidth=0,edgecolor='k',label=leg2)
ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))



plt.savefig('../figures/PFR.png',format='png',dpi=300,bbox_inches='tight')

