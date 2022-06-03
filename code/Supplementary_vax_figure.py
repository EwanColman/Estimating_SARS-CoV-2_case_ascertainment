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
dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]
dates=dates+['01/0'+str(i)+'/2021' for i in range(1,10)]+['01/'+str(i)+'/2021' for i in range(10,13)]
dates=dates+['01/0'+str(i)+'/2022' for i in range(1,5)]

dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
  

fig = plt.figure(figsize=(12,16))
gs = fig.add_gridspec(8, 3, height_ratios=[2,2,2,1,2,2,2,2])

#fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)


######## AGES ##########################################################


#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(age).p','rb'))
vax_data=pk.load(open('../pickles/reporting_rates_(age)_VAX.p','rb'))


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
        plt.text(100,200,'A',size=20)
          

    reporting_multiplier=data['reporting_multiplier_'+age]
    reporting_multiplier_vax=vax_data['reporting_multiplier_'+age]
    #testing_multiplier=data['testing_multiplier_'+age]
    new_variant_proportion=data['sgtf_proportion_'+age]
    sgtf_proportion=data['sgtf_proportion_'+age]
    #plt.scatter(date,reporting_multiplier,s=10,facecolors='#006666')
    
    #if i==1:
    leg1='Original ascertainment rate, no vaccination effect'
    leg1_vax='Ascertainment rate assuming strong vaccination effect'
    leg2='Plausible region of improved estimate'

    ax.plot(date,reporting_multiplier['Rate'][start:end],linestyle=':',marker='o',markersize=2,color='k',linewidth=1,zorder=1,label=leg1)
    ax.plot(date,reporting_multiplier_vax['Rate'][start:end],color='k',linewidth=0.5,zorder=1,label=leg1_vax)
    
    ax.fill_between(date,reporting_multiplier['Rate'][start:end],reporting_multiplier_vax['Rate'][start:end],color='c',alpha=0.2,linewidth=1,label=leg2)
    

ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))


############# REGIONS ###################################
data=pk.load(open('../pickles/reporting_rates_(region).p','rb'))
vax_data=pk.load(open('../pickles/reporting_rates_(region)_VAX.p','rb'))


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
        plt.text(100,200,'B',size=20)

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
    
    reporting_multiplier=data['reporting_multiplier_'+region]
    reporting_multiplier_vax=vax_data['reporting_multiplier_'+region]
    #testing_multiplier=data['testing_multiplier_'+region]
    #new_variant_proportion=data['variant_proportion_'+region]
    sgtf_proportion=data['sgtf_proportion_'+region]
    
    ax.plot(date,reporting_multiplier['Rate'][start:end],linestyle=':',marker='o',markersize=2,color='k',linewidth=1,zorder=1,label=leg1)
    ax.plot(date,reporting_multiplier_vax['Rate'][start:end],color='k',linewidth=0.5,zorder=1,label=leg1)
    
    ax.fill_between(date,reporting_multiplier['Rate'][start:end],reporting_multiplier_vax['Rate'][start:end],color='c',alpha=0.2,linewidth=1,label=leg2)
    

plt.savefig('../figures/Supplementary_figure_vax.pdf',format='pdf',dpi=300,bbox_inches='tight')

