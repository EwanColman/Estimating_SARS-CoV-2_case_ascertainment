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


max_day=825
fs=12
start_date='1 August 2020'
end_date='1 June 2022'

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
start_date=(datetime.strptime(start_date, '%d %B %Y')-time_zero).days
end_date=(datetime.strptime(end_date, '%d %B %Y')-time_zero).days
###########################################

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(11,10,2)]+['01/'+str(i)+'/2020' for i in range(11,13,2)]
dates=dates+['01/0'+str(i)+'/2021' for i in range(1,10,2)]+['01/'+str(i)+'/2021' for i in range(11,13,2)]
dates=dates+['01/0'+str(i)+'/2022' for i in range(1,6,2)]

dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
  

fig = plt.figure(figsize=(12,16))
gs = fig.add_gridspec(8, 3, height_ratios=[2,2,2,1,2,2,2,2])

#fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)


######## AGES ##########################################################


#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(age)_200.p','rb'))




name_of={'02_10':'2 to 10',
         '11_15':'11 to 15',
         '16_24':'16 to 24',
         '25_34':'25 to 34',
         '35_49':'35 to 49',
         '50_69':'50 to 69',
         '70+':'70+'}

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
    
i=0
for age in name_of:
    
#    n=i % 3
#    m=1+int(n/3)
    
    #ax=fig.add_subplot(3,3,i) 
    ax = fig.add_subplot(gs[int(i/3),i%3])
    i=i+1
    ax.text(start_date+10,85,'Age '+name_of[age],size=fs)
    ax.set_xlim([start_date,end_date+10])
    #ax.axvline((datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,linewidth=1)
    ax.set_ylim([0,100])
    if i%3==1:
        ax.set_yticks([0,20,40,60,80])
        ax.set_yticklabels([0,20,40,60,80],size=fs)
    else:
        ax.set_yticks([])   
    if i>4:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i==4: 
        plt.ylabel('Percentage of infections',size=fs)
        plt.text(40,200,'A',size=20)
          

    # medians and confidence intervals
    reporting_multiplier_list=data['reporting_multiplier_'+age]
    reporting_multiplier={'Lower':[],'Rate':[],'Upper':[]}
    for t in range(len(reporting_multiplier_list[0])):
        rate_list=sorted([reporting_multiplier_list[j][t] for j in range(200)])
        
        reporting_multiplier['Lower'].append(rate_list[4])
        reporting_multiplier['Rate'].append(rate_list[99])
        reporting_multiplier['Upper'].append(rate_list[194])
    
    incidence_list=data['incidence_'+age]
    incidence={'Lower':[],'Rate':[],'Upper':[]}
    for t in range(len(incidence_list[0])):  
        inci_list=sorted([incidence_list[j][t] for j in range(200)])
        incidence['Lower'].append(inci_list[4])
        incidence['Rate'].append(inci_list[99])
        incidence['Upper'].append(inci_list[194])
    ##########
        
    #if i==1:
    leg1='Ascertainment rate'
    leg2='Ascertainment rate 95% confidence interval'

    ax.plot(date,reporting_multiplier['Rate'][start:end],linestyle=':',marker='o',markersize=2,color='k',linewidth=1,zorder=1,label=leg1)
    ax.fill_between(date,reporting_multiplier['Lower'][start:end],reporting_multiplier['Upper'][start:end],color='m',alpha=0.1,linewidth=1,label=leg2)
    
    
    ax2 = ax.twinx()
    ax2.set_xlim([start_date,end_date+10])
    ax2.set_ylim([0,2])
    if i in [3,6,7]:
        ax2.set_yticks([0,0.5,1,1.5])
        ax2.set_yticklabels([0,0.5,1,1.5],size=fs)
    else:
        ax2.set_yticks([])
    if i==6: 
        ax2.set_ylabel('Percentage of population',size=fs)
        
    #I=data['incidence_'+age]
    #plt.fill_between(range(Sep1,len(I)),[0 for i in range(Sep1,len(I))],I[Sep1:],color='k',linewidth=0,alpha=0.1)

    I=incidence['Rate'][0:end_date]
    ax2.plot(range(date[0],len(I)),I[date[0]:],color='g',linewidth=0.75,label='Daily new infections')


# fake plot
#ax.fill_between([0,0],[0,0],[0,0],color='m',linewidth=0,edgecolor='k',label=leg2)
ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))
ax2.legend(loc=2,ncol=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.2))


############# REGIONS ###################################
data=pk.load(open('../pickles/reporting_rates_(region)_200.p','rb'))

name_of={'London':'London',
        'SouthEast':'South East',
        'EastofEngland':'East of England',
        'SouthWest':'South West',
        'EastMidlands':'East Midlands',
        'WestMidlands':'West Midlands',
        'YorkshireandTheHumber':'Yorkshire and the Humber',
        'NorthEast':'North East',
        'NorthWest':'North West',
        'Wales':'Wales',
        'Scotland':'Scotland',
        'NorthernIreland':'Northern Ireland'
        }

initials_of={'London':'LDN',
        'SouthEast':'SE',
        'EastofEngland':'EoE',
        'SouthWest':'SW',
        'EastMidlands':'EM',
        'WestMidlands':'WM',
        'YorkshireandTheHumber':'Y&H',
        'NorthEast':'NE',
        'NorthWest':'NW'
        }

regions=['SouthWest',
        'London',
        'SouthEast',
        'WestMidlands',
        'EastMidlands',
        'EastofEngland',
        'NorthWest',
        'YorkshireandTheHumber',
        'NorthEast',
        'Wales',
        'Scotland',
        'NorthernIreland']

m=0
i=0
print()
for region in regions:
    
    ax = fig.add_subplot(gs[7-int(i/3),i%3])
    i=i+1
    ax.text(start_date+10,85,name_of[region],size=fs)
    ax.set_xlim([start_date,end_date+10])
    ax.set_ylim([0,100])
    
    if i%3==1:
        ax.set_yticks([0,20,40,60,80])
        ax.set_yticklabels([0,20,40,60,80],size=fs)
    
    else:
        ax.set_yticks([])
        
    if i<4:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i==7: 
        plt.ylabel('Percentage of infections',size=fs)
        plt.text(40,200,'B',size=20)

    # get index of sep 1
    date=data['date_'+region]
    start=0
    while date[start]<start_date:
        start=start+1

    end=len(date)-1
    while date[end]>end_date:
        end=end-1
    # adust the date range
    date=date[start:end]
    
    
    # medians and confidence intervals
    reporting_multiplier_list=data['reporting_multiplier_'+region]
    reporting_multiplier={'Lower':[],'Rate':[],'Upper':[]}
    for t in range(len(reporting_multiplier_list[0])):
        rate_list=sorted([reporting_multiplier_list[j][t] for j in range(200)])
        
        reporting_multiplier['Lower'].append(rate_list[4])
        reporting_multiplier['Rate'].append(rate_list[99])
        reporting_multiplier['Upper'].append(rate_list[194])
    
    incidence_list=data['incidence_'+region]
    incidence={'Lower':[],'Rate':[],'Upper':[]}
    for t in range(len(incidence_list[0])):  
        inci_list=sorted([incidence_list[j][t] for j in range(200)])
        incidence['Lower'].append(inci_list[4])
        incidence['Rate'].append(inci_list[99])
        incidence['Upper'].append(inci_list[194])
    ##########
    

    ax.plot(date,reporting_multiplier['Rate'][start:end],linestyle=':',marker='o',markersize=2,color='k',linewidth=1,zorder=1)
    ax.fill_between(date,reporting_multiplier['Lower'][start:end],reporting_multiplier['Upper'][start:end],color='m',alpha=0.1,linewidth=1)
    

    ax2 = ax.twinx()
    ax2.set_xlim([start_date,end_date+10])
    ax2.set_ylim([0,2])
    if i%3==0:
        ax2.set_yticks([0,0.5,1,1.5])
        ax2.set_yticklabels([0,0.5,1,1.5],size=fs)
    else:
        ax2.set_yticks([])
    if i==9: 
        ax2.set_ylabel('Percentage of population',size=fs)
      
    I=incidence['Rate'][0:end_date]
    ax2.plot(range(date[0],len(I)),I[date[0]:],color='g',linewidth=0.75)


plt.savefig('../figures/figure2.pdf',format='pdf',dpi=300,bbox_inches='tight')

