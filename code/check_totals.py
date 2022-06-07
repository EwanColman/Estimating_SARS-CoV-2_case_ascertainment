
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt


#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

population_of_England=56287000


df=pd.read_csv('../processed_data/England_daily_data.csv')
total_cases_total=df['all_cases'].tolist()


dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,7)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]
dates_words2=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')

LFD_proportion=df['LFD_proportion'].tolist()


######## AGES ####################################
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

folder='../raw_data/'

population_of={'02_10':6254603,
                '11_15':3370248,
                '16_24':5950637,
                '25_34':7596145,
                '35_49':10853151,
                '50_69':13618246,
                '70+':7679719
               }

print('Eng pops:',population_of_England)
print('age pops:',sum(population_of.values()))


total_ages={}

total_cases_age=[0 for i in range(1000)]


for age in population_of:
    print(age)
    df=pd.read_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
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
    
    df=pd.read_csv('../raw_data/Diagnostic/'+age+'_cases.csv')
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['date'].tolist()]
    print(time_in_days)
    df['days_since_march1']=time_in_days
    df=df.sort_values('days_since_march1',ascending=True)
    df=df[df['days_since_march1']>=0]
    cases=df['cases'].tolist()
    case_time_in_days=df['days_since_march1'].tolist()
    dates=df['date'].tolist()

    total_cases_age=[total_cases_age[t]+cases[t] for t in range(len(cases))]


######## REGIONS   ##############################

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
              # 'Wales':3152879,
              # 'Scotland':5463300,
              # 'NorthernIreland':1893667
               }

print('region pops:',sum(population_of.values()))

total_regions={}
total_cases_region=[0 for i in range(1000)]


for region in population_of:
    
    df=pd.read_csv('../raw_data/Surveillance/'+region+'_surveillance_1k.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
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
    
    total_cases_region=[total_cases_region[t]+cases[t] for t in range(len(cases))]

for i in range(len(total_cases_region)):
    print(i,total_cases_age[i],total_cases_region[i],total_cases_total[i])


plt.figure()
plt.plot(total_cases_age,label='age')
plt.plot(total_cases_region,label='region')
plt.plot(total_cases_total,label='all')

plt.figure()
plt.plot(ONS_day,[total_ages[day] for day in ONS_day],label='age')
plt.plot(ONS_day,[total_regions[day] for day in ONS_day],label='region')


#for day in total_ages:
#    print(day,total_ages[day],total_regions[day],total_total[day])



