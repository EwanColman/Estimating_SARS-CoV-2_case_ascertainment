import pandas as pd
from datetime import datetime
from math import isnan
# for correlations use the last T points

folder='../raw_data/Regions'

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
population_of={'London':8961989,
               'SouthEast':9180135,
               'EastofEngland':6236072,
               'SouthWest':5624696,
               'EastMidlands':4835928,
               'WestMidlands':5934037,
               'YorkshireandTheHumber':5502967,
               'NorthEast':2669941,
               'NorthWest':7341196
               }

name_of={'London':'London',
        'SouthEast':'South East',
        'EastofEngland':'East of England',
        'SouthWest':'South West',
        'EastMidlands':'East Midlands',
        'WestMidlands':'West Midlands',
        'YorkshireandTheHumber':'Yorkshire and the Humber',
        'NorthEast':'North East',
        'NorthWest':'North West'
        }

#### Variant proportion df for correlations 


columns=['Week starting','N only','OR only','S only','OR+N','OR+S','N+S','OR+N+S','Mean','10th Percentile','25th Percentile','50th Percentile','75th Percentile','90th Percentile']  

regions=['NorthEast',
      'NorthWest',
      'YorkshireandTheHumber',
      'EastMidlands',
      'WestMidlands',
      'EastofEngland',
      'London',
      'SouthEast',
      'SouthWest']

variant_proportion={}

start_date='1 November 2020'
start_day=(datetime.strptime(str(start_date), '%d %B %Y')-time_zero).days
end_date='1 August 2021'
end_day=(datetime.strptime(str(end_date), '%d %B %Y')-time_zero).days
Sep20='20 September 2020'
Sep20=(datetime.strptime(str(Sep20), '%d %B %Y')-time_zero).days



rows={'NorthEast':[7,48],
      'NorthWest':[51,94],
      'YorkshireandTheHumber':[97,140],
      'EastMidlands':[143,186],
      'WestMidlands':[189,232],
      'EastofEngland':[235,278],
      'London':[281,324],
      'SouthEast':[327,370],
      'SouthWest':[373,412]}
#n=32
for region in regions:

    df=pd.read_excel('../raw_data/ONS_technical.xlsx',sheet_name='1b',usecols='A:N',names=columns,skiprows=rows[region][0],nrows=rows[region][1]-rows[region][0]-1)
    
    # # add a new column with header days since Mar1
    time_in_days=[]
    date=df['Week starting'].tolist()
    # printthis to check its getting all the data
    print(region,'From',date[0],'to',date[-1])
    while date:
        d=date.pop(0)
        #print(d)
        #day_numerical=(datetime.strptime(str(d), '%d %B %Y')-time_zero).days
        day_numerical=(datetime.strptime(str(d), '%Y-%m-%d %H:%M:%S')-time_zero).days
    
        # minus 4 to make it a mid-week estimate
        time_in_days.append(day_numerical)  
    # add it to the dataframe
    df['days_since_march1']=time_in_days     
    # get a df for the ct threshold
    sep_df=df[df['days_since_march1']>Sep20]
    # and one for VOC proportion 
    df=df[df['days_since_march1']>start_day]
    
    
    time_in_days=df['days_since_march1'].tolist()
    #proportion=df['OR+N'].tolist()
    proportion=(100*df['OR+N']/(df['OR+N']+df['OR+N+S']+df['N+S']+df['OR+S']+df['S only'])).tolist()
    
    for i in range(len(proportion)):
        if isnan(proportion[i]):
            proportion[i]=proportion[i-1]
    
    
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

        for i in range(delta_t):
            prop=prop+delta_p/delta_t           
            variant_proportion.append(prop)
 
            
    while len(variant_proportion)<end_day:
        variant_proportion.append(prop)                                                         

      


    time_in_days=sep_df['days_since_march1'].tolist()
    CT=sep_df['Mean'].tolist()
    OR=sep_df['OR only'].tolist()

    mean_CT=[0 for i in range(Sep20)]
    OR_only=[0 for i in range(Sep20)]
    
    t=Sep20
    c=0
    ct=0
    o=0
    ornf=0
    while time_in_days:
        last_t=t
        last_c=c
        last_o=o
        
        t=time_in_days.pop(0)
        c=CT.pop(0)
        o=OR.pop(0)
        
        delta_t=t-last_t
        delta_c=c-last_c
        delta_o=o-last_o
        
        for i in range(delta_t):
            ct=ct+delta_c/delta_t
            
            ornf=ornf+delta_o/delta_t
            
            mean_CT.append(ct)
            OR_only.append(ornf)
            
    while len(mean_CT)<end_day:                                                           
        mean_CT.append(ct)
        OR_only.append(ornf)
    
    
    ###### CASES ##########
    
    df=pd.read_csv(folder+'/'+region+'_cases.csv')
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
    cases=df['newCasesBySpecimenDate'].tolist()
    case_time_in_days=df['days_since_march1'].tolist()
    dates=df['date'].tolist()
    
    #### PCR TESTS #######
    df=pd.read_csv(folder+'/'+region+'_PCR.csv')
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
    
    ### LF TESTS #######
    df=pd.read_csv(folder+'/'+region+'_LF.csv')
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
     
    ###### Make sreadsheet ##########
    end_day=min(len(cases),len(PCR_tests),len(LFD_tests))
    #print(len(dates),len(cases),len(LFD_tests),len(PCR_tests),len(variant_proportion))

    output_df=pd.DataFrame({'Date':dates[0:end_day],
                        'Days_since_March1':case_time_in_days[0:end_day],
                        'cases':cases[0:end_day],
                        'PCR_tests':PCR_tests[0:end_day],
                        'LFD_tests':LFD_tests[0:end_day],
                        'NV_proportion':variant_proportion[0:end_day],
                        'mean_CT':mean_CT[0:end_day],
                        'OR_only':OR_only[0:end_day]})

    output_df.to_csv('../processed_data/'+region+'_daily_data.csv',index=False)



