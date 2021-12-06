
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

folder='../raw_data/Regions'

max_day=506

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)

name_of={'NorthWest':'North West',
         'YorkshireandTheHumber':'Yorkshire and The Humber',
        'NorthEast':'North East',
        'WestMidlands':'West Midlands',
        'EastMidlands':'East Midlands',
        'EastofEngland':'East of England',
        'SouthWest':'South West',
        'London':'London',
        'SouthEast':'South East'}


# #################

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,5)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]


### plot template ####
fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,4,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,5.2])
    if i%3==1:
        ax.set_yticks([0,1,2,3,4])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Percentage of people')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####


    # get surveillance data for the region    
    df=pd.read_csv(folder+'/'+region+'_surveillance_1e.csv',sep=',')    
    # add a new column with header days since Mar1
    time_in_days=[]
    date=df['Date'].tolist()
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        # minus 4 to make it a mid-week estimate
        time_in_days.append(day_numerical)
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    
    # convert data frame to lists for plotting
    rate=df['Rate'].tolist()
    upper=df['Upper'].tolist()
    lower=df['Lower'].tolist()
    ONS_day=df['days_since_march1'].tolist()
    date=df['Date'].tolist()
    
    plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
    plt.scatter(ONS_day,rate,s=10,color='k',marker='^')  

plt.savefig('../raw_data_figures/ONS_number_test-positive_regions.png',format='png', bbox_inches='tight',dpi=256)



# make the plot of the Case data
### plot template ####


fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,12000,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,18000])
    if i%3==1:
        ax.set_yticks([i*5000 for i in range(4)])
        ax.set_yticklabels([i*5 for i in range(4)])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Thousands of cases')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####

    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    cases=df['cases'].tolist()
    time_in_days=df['Days_since_March1'].tolist()

    plt.bar(time_in_days,cases,width=1,color='k')
    
plt.savefig('../raw_data_figures/Pillar_cases_regions.png',format='png', bbox_inches='tight',dpi=256)

'''
fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,400000,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,500000])
    if i%3==1:
        ax.set_yticks([i*100000 for i in range(5)])
        ax.set_yticklabels([round(i*0.1,1) for i in range(5)])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Millions of tests')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####

    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    PCR_tests=df['PCR_tests'].tolist()
    time_in_days=df['Days_since_March1'].tolist()

    plt.bar(time_in_days,PCR_tests,width=1,color='k')
    
plt.savefig('../raw_data_figures/PCR_tests_regions.png',format='png', bbox_inches='tight',dpi=256)

### LFD tests ###########
fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,120000,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,150000])
    if i%3==1:
        ax.set_yticks([i*40000 for i in range(4)])
        ax.set_yticklabels([i*40 for i in range(4)])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Thousands of tests')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####

    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    LFD_tests=df['LFD_tests'].tolist()
    time_in_days=df['Days_since_March1'].tolist()

    plt.bar(time_in_days,LFD_tests,width=1,color='k')
    
plt.savefig('../raw_data_figures/LFD_tests_regions.png',format='png', bbox_inches='tight',dpi=256)
'''

### New variant proportion ###########
fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,85,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,102])
    if i%3==1:
        ax.set_yticks([i*20 for i in range(1,6)])
        #ax.set_yticklabels([i*20 for i in range(4)])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Percentage')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####

    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    variant_proportion=df['NV_proportion'].tolist()
    time_in_days=df['Days_since_March1'].tolist()
    
    x_values=[time_in_days[i] for i in range(0,len(variant_proportion)-1,7)]
    proportion=[100*variant_proportion[i] for i in range(0,len(variant_proportion)-1,7)]
          
    plt.plot(x_values,proportion,color='k')

plt.savefig('../raw_data_figures/new_variant_proportion_regions.png',format='png', bbox_inches='tight',dpi=256)


### Cycle threshold ###########
fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)
m=0
i=0
for region in name_of:    
    i=i+1 
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,35,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([20,40])
    if i%3==1:
        ax.set_yticks([20+i*5 for i in range(4)])
        #ax.set_yticklabels([i*20 for i in range(4)])
    else:
        ax.set_yticks([])        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])   
    if i==4: 
        plt.ylabel('Cycle threshold')
    # if i==8:
    #     plt.xlabel('Date')
#### Plot template end ####

    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    CT_mean=df['mean_CT'].tolist()
    time_in_days=df['Days_since_March1'].tolist()
    
    x_values=[time_in_days[i] for i in range(205,len(CT_mean)-1,7)]
    ct_value=[CT_mean[i] for i in range(205,len(CT_mean)-1,7)]
          
    plt.plot(x_values,ct_value,color='k')

plt.savefig('../raw_data_figures/cycle_threshold_regions.png',format='png', bbox_inches='tight',dpi=256)


