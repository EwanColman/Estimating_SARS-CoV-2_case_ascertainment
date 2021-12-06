
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt


#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

population_of_England=56287000

df=pd.read_csv('../raw_data/ONS_incidence.csv')
#df.drop(df.tail(1).index,inplace=True)


# add a new column with header days since Mar1
time_in_days=[]
date=df['Date'].tolist()
while date:
    d=date.pop(0)
    #print(d)
    day_numerical=(datetime.strptime(str(d), '%d-%b-%y')-time_zero).days
    # minus 4 to make it a mid-week estimate
    time_in_days.append(day_numerical-3)
df['days_since_march1']=time_in_days

rate=df['Rate'].tolist()
upper=df['Upper'].tolist()
lower=df['Lower'].tolist()
ONS_day_eng=df['days_since_march1'].tolist()
date=df['Date'].tolist()


total_total={}
for i in range(len(ONS_day_eng)):
    day=ONS_day_eng[i]
    total_total[day]=rate[i]*population_of_England/100




df=pd.read_csv('../processed_data/England_daily_data.csv')
total_cases_total=df['all_cases'].tolist()


dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,7)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]
dates_words2=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
variant_proportion=df['NV_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()


######## AGES ####################################
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

folder='../raw_data/'

population_of={'02_10':6264662,
               '11_15':2664513,
               '16_24':5860347,
               '25_34':7567108,
               '35_49':10864046,
               '50_69':13688509,
               '70+':8114862
               }

print('Eng pops:',population_of_England)
print('age pops:',sum(population_of.values()))


total_ages={}

total_cases_age=[0 for i in range(1000)]


for age in population_of:
   
    df=pd.read_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    
    date=df['Date'].tolist()
    
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
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
    #date=df['Date'].tolist()
    
    for i in range(len(ONS_day)):
        day=ONS_day[i]
        if not day in total_ages:
            total_ages[day]=0
        total_ages[day]=total_ages[day]+rate[i]*population_of[age]/100

    ####### CASES ###############
    print(age)
    df=pd.read_csv('../raw_data/Diagnostic/'+age+'_cases.csv')
    #print(df.head())
    #print()
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
    cases=df['cases'].tolist()
    case_time_in_days=df['days_since_march1'].tolist()
    dates=df['date'].tolist()

    total_cases_age=[total_cases_age[t]+cases[t] for t in range(len(cases))]


######## REGIONS   ##############################


# use different axis labels


# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
population_of={'SouthWest':5624696,
               'London':8961989,
               'SouthEast':9180135,
               'WestMidlands':5934037,
               'EastMidlands':4835928,
               'EastofEngland':6236072,
               'NorthWest':7341196,
               'YorkshireandTheHumber':5502967,
               'NorthEast':2669941,
               }

print('region pops:',sum(population_of.values()))

total_regions={}
total_cases_region=[0 for i in range(1000)]


for region in population_of:
    
    df=pd.read_csv('../raw_data/Surveillance/'+region+'_surveillance_1k.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    
    date=df['Date'].tolist()
    
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
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
    
    for i in range(len(ONS_day)):
        day=ONS_day[i]
        if not day in total_regions:
            total_regions[day]=0
        total_regions[day]=total_regions[day]+rate[i]*population_of[region]/100

    ###### CASES ##########
    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    #print(df.head())
    #print(len(df))
    cases=df['cases'].tolist()
    time_in_days=df['Days_since_March1'].tolist()
    variant_proportion=df['NV_proportion'].tolist()
    
    total_cases_region=[total_cases_region[t]+cases[t] for t in range(len(cases))]

for i in range(len(total_cases_age)):
    print(i,total_cases_age[i],total_cases_region[i],total_cases_total[i])


plt.figure()
plt.plot(total_cases_age,label='age')
plt.plot(total_cases_region,label='region')
plt.plot(total_cases_total,label='all')

plt.figure()
plt.plot(ONS_day,[total_ages[day] for day in ONS_day],label='age')
plt.plot(ONS_day,[total_regions[day] for day in ONS_day],label='region')
plt.plot(ONS_day_eng,[total_total[day] for day in ONS_day_eng],label='all')



#for day in total_ages:
#    print(day,total_ages[day],total_regions[day],total_total[day])



