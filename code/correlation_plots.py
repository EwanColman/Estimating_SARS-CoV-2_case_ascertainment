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
Sep1='1 August 2020'
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days

###########################################
#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,8)]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
  

fig = plt.figure(figsize=(12,16))
gs = fig.add_gridspec(3, 1, height_ratios=[6,1,8])

#fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)


######## AGES ##########################################################


#folder='../raw_data/Regions'
data=pk.load(open('../pickles/reporting_rates_(age).p','rb'))

ages=['02_10','11_15','16_24','25_34','35_49','50_69','70+']

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
while date[start]<Sep1:
    start=start+1

date=date[start:]
    


print()
ax = fig.add_subplot(gs[0,0])
plt.text(-1.3,5.5,'A',size=20)

M=[]
#for region1 in regions:
for i in range(1,7):
    

    multiplier1=data['reporting_multiplier_'+ages[i]]['Rate']
    m=[]
    for j in range(i):
        multiplier2=data['reporting_multiplier_'+ages[j]]['Rate']
        print(len(multiplier1))
        print(len(multiplier2))
        
        pearson, p=stats.pearsonr(multiplier1,multiplier2)
        print(name_of[ages[i]],'vs',name_of[ages[j]],'p=',round(p,2),round(pearson,2)) 
        if p<0.05 and pearson>0:
            m.append(pearson)
            plt.text(j-0.25,i-1.07,str(round(pearson,2)),size=fs)
        else:
            m.append(0)
            #plt.text(8-i,8-j,'p<0.05')
        
        
    for j in range(i,6):
        m.append(0)

    M.append(m)
    

Values=np.array(M)

plt.xticks(range(6),[initials_of[r] for r in ages[0:6]],size=fs,rotation=0)
plt.yticks(range(6),[initials_of[r] for r in ages[1:7]],size=fs)
ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.imshow(Values,cmap='Greys',vmax=1.5,vmin=0,origin='lower') # displays in color


############# REGIONS ###################################
data=pk.load(open('../pickles/reporting_rates_(region).p','rb'))


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
        'NorthWest':'NW',
        'Wales':'WA',
        'Scotland':'SCO',
        'NorthernIreland':'NI'
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


# data={}
# for region in regions:
#     df=pd.read_csv(folder+'/'+region+'_surveillance_1e.csv',sep=',') 
#     data[region]=df['Rate'].tolist()

ax = fig.add_subplot(gs[2,0])
plt.text(-1.5,10.8,'B',size=20)


s=15
#ax.set_title('Randomized edges',size=fs,pad=10)
#plt.ylabel('VoC reporting rate',size=fs)
#plt.xlabel('Peak sensitivity',size=fs)
regions=['SouthWest',
        'London',
        'SouthEast',
        'EastofEngland',
        'EastMidlands',
        'WestMidlands',
        'NorthWest',
        'YorkshireandTheHumber',
        'NorthEast',
        'Wales',
        'Scotland',
        'NorthernIreland']


M=[]
#for region1 in regions:
for i in range(1,12):
    print()
    print(regions[i])

    # regions have one more week than nations so remove final point
    multiplier1=data['reporting_multiplier_'+regions[i]]['Rate'][:115]
    

    
    m=[]
    for j in range(i):
        # regions have one more week than nations so removefinal point (117)
        multiplier2=data['reporting_multiplier_'+regions[j]]['Rate'][:115]
        


        
        
        
        pearson, p=stats.pearsonr(multiplier1,multiplier2)
        print(name_of[regions[i]],'vs',name_of[regions[j]],'p=',round(p,2),round(pearson,2)) 
        if p<0.05 and pearson>0:
            m.append(pearson)
            plt.text(j-0.3,i-1.1,str(round(pearson,2)),size=fs)
        else:
            m.append(0)
            #plt.text(8-i,8-j,'p<0.05')
        
        
    for j in range(i,11):
        m.append(0)

    M.append(m)
    

Values=np.array(M)
plt.xticks(range(11),[initials_of[r] for r in regions[0:11]],size=fs,rotation=0)
plt.yticks(range(11),[name_of[r] for r in regions[1:12]],size=fs)
ax.xaxis.tick_top()
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

plt.imshow(Values,cmap='Greys',vmax=1.5,vmin=0,origin='lower') # displays in color


plt.savefig('../figures/figureS3.pdf',format='pdf',dpi=300,bbox_inches='tight')

