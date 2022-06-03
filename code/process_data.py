import pandas as pd
from datetime import datetime


#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

###### CASES ##########
df=pd.read_csv('../raw_data/England_cases.csv').fillna(0)
df=df[df['areaName']=='England'] 
# add a new column with header days since Mar1
time_in_days=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['date'].tolist()]

df['days_since_march1']=time_in_days
df=df.sort_values('days_since_march1',ascending=True)
df=df[df['days_since_march1']>=0]
cases=df['newCasesPCROnlyBySpecimenDate'].tolist()

LFD_only_cases=df['newCasesLFDOnlyBySpecimenDate'].tolist()
all_cases=df['newCasesBySpecimenDate'].tolist()
# calculate the proportion of positives that that are LFD only
LFD_proportion=[LFD_only_cases[i]/all_cases[i] for i in range(len(cases))]

dates=df['date'].tolist()

# Get variant proportions for England
# https://www.ons.gov.uk/peoplepopulationandcommunity/healthandsocialcare/conditionsanddiseases/datasets/coronaviruscovid19infectionsurveydata

# make it 0 before this date
start_date='1 November 2020'
start_day=(datetime.strptime(str(start_date), '%d %B %Y')-time_zero).days
end_day=len(cases)

#columns=['Week starting','N protein only','ORF1ab only','S protein only','ORF1ab + N protein','ORF1ab + S protein','N protein + S protein','ORF1ab + N protein + S protein','Mean','10th Percentile','25th Percentile','50th Percentile','75th Percentile','90th Percentile']

#old version
columns=['Week starting','N protein only','OR only','S only','OR+N','OR+S','N+S','OR+N+S','Mean','10th Percentile','25th Percentile','50th Percentile','75th Percentile','90th Percentile']  

#####
#xl = pd.ExcelFile('../raw_data/ONS_technical.xlsx')
#print(xl.sheet_names)


rows=[117,226]
df=pd.read_excel('../raw_data/ONS_technical.xlsx',sheet_name='1a',usecols='A:N',names=columns,skiprows=rows[0],nrows=rows[1]-rows[0]-1) 
# add a new column with header days since Mar1
time_in_days=[]
date=df['Week starting'].tolist()
# printthis to check its getting all the data
print('From',date[0],'to',date[-1])
time_in_days=[(datetime.strptime(str(d), '%Y-%m-%d %H:%M:%S')-time_zero).days for d in date] 

# add it to the dataframe
df['days_since_march1']=time_in_days      
df=df[df['days_since_march1']>start_day]

time_in_days=df['days_since_march1'].tolist()
#proportion=(df['OR+N']+df['N only']+df['OR only']).tolist()
#proportion=(df['OR+N']+df['OR only']).tolist()
# proportion=df['OR+N'].tolist()

# proportion is the ration of S+N+OR to OR+N
proportion=(100*df['OR+N']/(df['OR+N']+df['OR+N+S']+df['N+S']+df['OR+S']+df['S only'])).tolist()


# interpolation 
sgtf_proportion=[0 for i in range(start_day)]
t=start_day
p=0
prop=0
while time_in_days:
    last_t=t
    last_p=p
    
    t=time_in_days.pop(0)
    p=proportion.pop(0)/100
    delta_t=t-last_t
    delta_p=p-last_p
    #print(delta_p,delta_t)
    for i in range(delta_t):
        prop=prop+delta_p/delta_t
        
        #print(delta_p/delta_t)
        sgtf_proportion.append(prop)

# for the remaining time keep the final proportion
while len(sgtf_proportion)<end_day:
    sgtf_proportion.append(prop)

'''
#### PCR and LFD tests #######################
#### PCR TESTS #######
df=pd.read_csv('../raw_data/England_PCR.csv')
#print(df.head())

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
PCR_tests=[0 for i in range(min(time_in_days))]+df['uniquePeopleTestedBySpecimenDateRollingSum'].tolist()

print(len(PCR_tests))

#### LF TESTS #######
df=pd.read_csv('../raw_data/England_LFD.csv')
#print(df.head())

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
LFD_tests=[0 for i in range(min(time_in_days))]+df['newLFDTests'].tolist()
LFD_tests=LFD_tests[0:end_day]
print(len(LFD_tests))
'''
#end_day=min(len(cases),len(LFD_tests),len(PCR_tests))

end_day=len(cases)
print(end_day)

output_df=pd.DataFrame({'Date':dates[0:end_day],
                        'cases':cases[0:end_day],
                        'all_cases':all_cases[0:end_day],
                        #'PCR_tests':PCR_tests[0:end_day],
                        #'LFD_tests':LFD_tests[0:end_day],
                        'SGTF_proportion':sgtf_proportion[0:end_day],
                        'LFD_proportion':LFD_proportion[0:end_day]})

output_df.to_csv('../processed_data/England_daily_data.csv',index=False)
