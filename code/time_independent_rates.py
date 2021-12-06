
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
import pickle as pk
import scipy.optimize as optimize

x1_values=[100/i for i in range(1,100)]
x2_values=[100/i for i in range(1,100)]
#x3_values=[100/i for i in range(20,55)]


# change for sensitivity analysis
#delta=1
delta=2

stretch=1
max_day=380

fs=12

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

#measure days from this day
#time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Sep1='23 August 2020'
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days

Mar1='1 March 2021'
Mar1=(datetime.strptime(str(Mar1), '%d %B %Y')-time_zero).days



dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,3)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]
dates_words2=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
variant_proportion=df['NV_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()


testable_period=30#5+int(30*stretch)

# create the output dictionary
output={'Group':[],'Wild_Rate':[],'Wild_Lower':[],'Wild_Upper':[],'Alpha_Rate':[],'Alpha_Lower':[],'Alpha_Upper':[]}

### Lognormal incubation period distribution
# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean,dispersion=5.5,1.52
mu,sigma=np.log(mean),np.log(dispersion)

R=[]
for i in range(testable_period):
    x=i+1
    u=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    x=i
    l=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    p=u-l
    R.append(p)

#### Probability of positive test #################
# the value in position i is th probability of testing positive i days after infection
#from https://www.medrxiv.org/content/10.1101/2020.11.24.20229948v1.full.pdf
df=pd.read_csv('../raw_data/PCR_curve_summary.csv')
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]
#S=S+[0 for i in range(30,testable_period)]
# in for old variants it should just be 0 after the end

df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]


######### COMMENT THIS OUT #######################
# for sensitivity analysis multiply by 1/S_pcr
#rescale_factor=1/max(S_pcr)
#S_pcr=[s*rescale_factor for s in S_pcr]
#S_lfd=[s*rescale_factor for s in S_lfd]
##################################################


# make the probability function as a list
P_pcr=[]
for j in range(delta,testable_period):
    numerator=R[j-delta]*S_pcr[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i-delta]*S_pcr[i] for i in range(testable_period)])
    P_pcr.append(numerator/denominator)   

# make the probability function as a list
P_lfd=[]
for j in range(testable_period):
    numerator=R[j]*S_lfd[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i]*S_lfd[i] for i in range(testable_period)])
    P_lfd.append(numerator/denominator) 


######## AGES ####################################
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

population_of={'02_10':6254603,
               '11_15':3370248,
               '16_24':5950637,
               '25_34':7596145,
               '35_49':10853151,
               '50_69':13618246,
               '70+':7679719
               }

name_of={'02_10':'2 to 10',
         '11_15':'11 to 15',
         '16_24':'16 to 24',
         '25_34':'25 to 34',
         '35_49':'35 to 49',
         '50_69':'50 to 69',
         '70+':'70+'}

#ages=['00_10','11_15','16_24','25_34','35_49','50_69','70+']

i=0
for age in name_of:
    print(age)
    output['Group'].append(age)

    df=pd.read_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    
    date=df['Date'].tolist()
    
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
        # minus 7 to get midpoint of the 2-week estimate
        time_in_days.append(day_numerical)
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    df=df[(df['days_since_march1']>Sep1)&(df['days_since_march1']<Mar1)] 
    
    # convert data frame to lists for plotting
    rate=df.reset_index()
    ONS_day=df['days_since_march1'].tolist()
    #date=df['Date'].tolist()
    
    ###### CASES ##########
    df=pd.read_csv('../raw_data/Diagnostic/'+age+'_cases.csv')

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
    
    C=cases+[0 for i in range(100)]
        
    theta_times_I=[]
    for t in range(len(cases)):
        Z_pcr=sum([C[t+j]*P_pcr[j-delta] for j in range(delta,testable_period)])
        Z_lfd=sum([C[t+j]*P_lfd[j] for j in range(testable_period)])
        
        Z=(1-LFD_proportion[t])*Z_pcr+LFD_proportion[t]*Z_lfd
        theta_times_I.append(Z)
        
        # set a list for normal cases, and an alternative one for variant cases
    new_cases=[]
    nv_cases=[]
    proportion=[]
    # estimate the number of test-positives
    for t in ONS_day:
        # j is the time since exposure
        estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        new_cases.append(estimate)

        proportion.append(variant_proportion[t])
    
    
    
    ## loop over the three lines
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    
    #for interval in reporting_multiplier:
    
    for interval in ['Rate','Lower','Upper']:
        new_ONS=[population_of[age]*r/100 for r in rate[interval]]
        ############# END ############                           
    
        #### Calculate th multiplier ###########
        # x is the multiplier
        x_values=[100/i for i in range(1,101)]
        #error=[]
        e_min=10**10
        for x1 in x1_values:
            for x2 in x2_values:
                #for x3 in x3_values:
                #x2=x1 
                #### Calculate an error #######
                #absolute error
                #x=[x1 for i in range(delta_index)]#+[x3 for i in range(delta_index,len(new_ONS))]
                
                e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
                # squared error
                #e=(sum([(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i])**2 for i in range(len(new_ONS))]))**(1/2)
                #print(int(100/x1),int(100/x2),np.log(e))
                #error.append(e)
                if e<e_min:
                    e_min=e
                    best_x=(x1,x2)
                    
                    best_new_cases=new_cases.copy()
                    #best_nv_cases=nv_cases.copy()
        
        # take the value with the smallest error
        x1,x2=best_x
         
        adjusted_cases=[x1*best_new_cases[i]*(1-proportion[i])+x2*best_new_cases[i]*proportion[i] for i in range(len(new_ONS))]
        print(age,interval,x1,x2)
        print(adjusted_cases)
        print()
        
        output['Alpha_'+interval].append(round(100/x2))
        output['Wild_'+interval].append(round(100/x1))       
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
               'Wales':3152879,
               'Scotland':5463300,
               'NorthernIreland':1893667
               }

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


#region='London'
m=0
i=0
for region in population_of:
    print(region)
    output['Group'].append(region)

    # get surveillance data for the region    
    df=pd.read_csv('../raw_data/Surveillance/'+region+'_surveillance_1k.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    
    date=df['Date'].tolist()
    
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
        # minus 4 to make it a mid-week estimate
        if region in ['Wales','Scotland','NorthernIreland']:
            time_in_days.append(day_numerical)
        else:
            time_in_days.append(day_numerical)
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    df=df[(df['days_since_march1']>Sep1)&(df['days_since_march1']<Mar1)] 
    
    
    # convert data frame to lists for plotting
    rate=df.reset_index()
    ONS_day=df['days_since_march1'].tolist()
    date=df['Date'].tolist()
    
    ###### CASES ##########
    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    #print(df.head())
    #print(len(df))
    cases=df['cases'].tolist()
    #time_in_days=df['Days_since_March1'].tolist()
    variant_proportion=df['NV_proportion'].tolist()
    
    ############## COPIED IN ###################
    C=cases+[0 for i in range(100)]

    theta_times_I=[]
    for t in range(len(cases)):
        Z_pcr=sum([C[t+j]*P_pcr[j-delta] for j in range(delta,testable_period)])
        Z_lfd=sum([C[t+j]*P_lfd[j] for j in range(testable_period)])
        
        Z=(1-LFD_proportion[t])*Z_pcr+LFD_proportion[t]*Z_lfd
        theta_times_I.append(Z)

    new_cases=[]
    nv_cases=[]
    proportion=[]
    # estimate the number of test-positives
    for t in ONS_day:
        # j is the time since exposure
        estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        new_cases.append(estimate)

        # to use the variant proportion
        proportion.append(variant_proportion[t])
    
    for interval in ['Rate','Lower','Upper']:
        new_ONS=[population_of[region]*r/100 for r in rate[interval]]
        ############# END ############
                            
        #### Calculate th multiplier ###########
        # x is the multiplier
    
        #error=[]
        e_min=10**10
        for x1 in x1_values:
            for x2 in x2_values:
                #for x3 in x3_values:
                    #x2=x1 
                    #### Calculate an error #######
                    #absolute error
                    #x=[x1 for i in range(delta_index)]+[x3 for i in range(delta_index,len(new_ONS))]
            
                e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
                # squared error
                #e=(sum([(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i])**2 for i in range(len(new_ONS))]))**(1/2)
                #print(int(100/x1),int(100/x2),np.log(e))
                #error.append(e)
                if e<e_min:
                    e_min=e
                    best_x=(x1,x2)
                    
                    best_new_cases=new_cases.copy()
                #best_nv_cases=nv_cases.copy()
    
        # take the value with the smallest error
        x1,x2=best_x                
        #x=[x1 for i in range(delta_index)]+[x3 for i in range(delta_index,len(new_ONS))]                
        adjusted_cases=[x1*best_new_cases[i]*(1-proportion[i])+x2*best_new_cases[i]*proportion[i] for i in range(len(new_ONS))]
        
        print(region,interval,x1,x2)
        print(adjusted_cases)
        print()
        output['Alpha_'+interval].append(round(100/x2))
        output['Wild_'+interval].append(round(100/x1))  
        
        
#pd.DataFrame(output).to_csv('../output/time_independent_rates.csv',index=False)
#pd.DataFrame(output).to_csv('../output/time_independent_rates_(higher_sensitivity).csv',index=False)
pd.DataFrame(output).to_csv('../output/time_independent_rates_(longer_delay).csv',index=False)
