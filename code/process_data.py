import pandas as pd
from datetime import datetime


#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

###### CASES ##########
df=pd.read_csv('../raw_data/England_cases.csv').fillna(0)
df=df[df['areaName']=='England'] 

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
cases=df['newCasesPCROnlyBySpecimenDate'].tolist()

LFD_only_cases=df['newCasesLFDOnlyBySpecimenDate'].tolist()
all_cases=df['newCasesBySpecimenDate'].tolist()
# calculate the proportion of positives that that are LFD only
LFD_proportion=[LFD_only_cases[i]/all_cases[i] for i in range(len(cases))]

case_time_in_days=df['days_since_march1']
dates=df['date'].tolist()
#print(len(dates))
#print(len(cases))
# Get thenew variant proportion for England
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


rows=[92,174]
df=pd.read_excel('../raw_data/ONS_technical.xlsx',sheet_name='1a',usecols='A:N',names=columns,skiprows=rows[0],nrows=rows[1]-rows[0]-1) 
# add a new column with header days since Mar1
time_in_days=[]
date=df['Week starting'].tolist()
# printthis to check its getting all the data
print('From',date[0],'to',date[-1])
while date:
    d=date.pop(0)
    print(d)
    #day_numerical=(datetime.strptime(str(d), '%d %B %Y')-time_zero).days
    day_numerical=(datetime.strptime(str(d), '%Y-%m-%d %H:%M:%S')-time_zero).days
    #day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
    
    # minus 4 to make it a mid-week estimate
    time_in_days.append(day_numerical)  
# add it to the dataframe
df['days_since_march1']=time_in_days      
df=df[df['days_since_march1']>start_day]

time_in_days=df['days_since_march1'].tolist()
#proportion=(df['OR+N']+df['N only']+df['OR only']).tolist()
#proportion=(df['OR+N']+df['OR only']).tolist()
# proportion=df['OR+N'].tolist()

# proportion is the ration of S+N+OR to OR+N
proportion=(100*df['OR+N']/(df['OR+N']+df['OR+N+S']+df['N+S']+df['OR+S']+df['S only'])).tolist()
# CT=df['Mean'].tolist()

# mean_CT=[]
# variant_proportion=[]
# next_change=time_in_days.pop(0)
# next_p=proportion.pop(0)
# next_ct=CT.pop(0)
# p=0
# ct=0
# for t in range(end_day):
#     mean_CT.append(ct)
#     variant_proportion.append(p/100)
#     if t==next_change:
#         p=next_p
#         ct=next_ct
#         if time_in_days:
#             next_change=time_in_days.pop(0)
#             next_p=proportion.pop(0)
#             next_ct=CT.pop(0)

variant_proportion=[0 for i in range(start_day)]
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
        variant_proportion.append(prop)
        
while len(variant_proportion)<end_day:
    variant_proportion.append(prop)

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

#end_day=min(len(cases),len(LFD_tests),len(PCR_tests))

end_day=len(cases)
print(end_day)

output_df=pd.DataFrame({'Date':dates[0:end_day],
                        #'Days_since_March1':case_time_in_days[0:end_day],
                        'cases':cases[0:end_day],
                        'all_cases':all_cases[0:end_day],
                        #'PCR_tests':PCR_tests[0:end_day],
                        #'LFD_tests':LFD_tests[0:end_day],
                        'NV_proportion':variant_proportion[0:end_day],
                        'LFD_proportion':LFD_proportion[0:end_day]})

output_df.to_csv('../processed_data/England_daily_data.csv',index=False)
