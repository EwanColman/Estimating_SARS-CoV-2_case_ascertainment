
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
import pickle as pk
import scipy.optimize as optimize

# finds the best fitting combination of multipliers
def hill_climb(new_cases,new_ONS):

    # initial starting point
    params={'Wild':10,'Alpha':10,'Delta':10,'BA.1':10,'BA.2':10}
    # to start the best oarams are the old params
    E_old=10**10
    
    E_best=E_old
    best_params=params.copy()
    # now look for something better
    maxima_found=False
    ###
    while not maxima_found:
        #loop over every free variable, find the best perturbation
        for variable in params:
            #print(variable,'of',len(var_list))
            # perturb one parameter
            perturbed_params=params.copy()
            # first try the down direction
            perturbed_params[variable]=max(0,params[variable]-1)
            # new residual error
            E_new=sum([abs(sum(new_cases[variant][i]*100/perturbed_params[variant] for variant in new_cases)-new_ONS[i]) for i in range(len(new_ONS))])
            # is it an improvement
            if E_new<E_best:
                E_best=E_new
                best_params=perturbed_params.copy()
        
            else:
                # if not then try the other direction
                perturbed_params=params.copy()
                # try the up direction
                perturbed_params[variable]=min(params[variable]+1,100)

                # new residual error
                E_new=sum([abs(sum(new_cases[variant][i]*100/perturbed_params[variant] for variant in new_cases)-new_ONS[i]) for i in range(len(new_ONS))])
               
                if E_new<E_best:
                    E_best=E_new
                    best_params=perturbed_params.copy()

        # after all that the best improvement an improvement?
        if E_old>E_best:
            params=best_params.copy()
            E_old=E_best
            
        # otherwise there is no way to improve
        else:
            maxima_found=True

    return params



# change for sensitivity analysis
#delta=2
delta=1

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

#measure days from this day
#time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
start_date='23 August 2020'
start_date=(datetime.strptime(str(start_date), '%d %B %Y')-time_zero).days

# create thestart and end points for each variant
variant_start={'Wild':0,
               'Alpha':(datetime.strptime('1 November 2020', '%d %B %Y')-time_zero).days,
               'Delta':(datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,
               'BA.1':(datetime.strptime('1 November 2021', '%d %B %Y')-time_zero).days,
               'BA.2':(datetime.strptime('9 January 2022', '%d %B %Y')-time_zero).days
               }

variant_end={'Wild':(datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,
               'Alpha':(datetime.strptime('1 November 2021', '%d %B %Y')-time_zero).days,
               'Delta':(datetime.strptime('9 January 2022', '%d %B %Y')-time_zero).days,
               'BA.1':(datetime.strptime('1 June 2022', '%d %B %Y')-time_zero).days,
               'BA.2':(datetime.strptime('1 June 2022', '%d %B %Y')-time_zero).days
               }


# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
SGTF_proportion=df['SGTF_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()


### get all the proportions for each variant across time
variant_proportion={}
for variant in ['Wild','Delta','BA.2']:
    variant_proportion[variant]=[0 for i in range(0,variant_start[variant])]\
                                +[1-SGTF_proportion[t] for t in range(variant_start[variant],variant_end[variant])]\
                                +[0 for i in range(variant_end[variant],len(SGTF_proportion))]
# for the S negatives
for variant in ['Alpha','BA.1']:
    variant_proportion[variant]=[0 for i in range(0,variant_start[variant])]\
                                +[SGTF_proportion[t] for t in range(variant_start[variant],variant_end[variant])]\
                                +[0 for i in range(variant_end[variant],len(SGTF_proportion))]


# create the output dictionary
output={'Group':[]}
for variant in variant_start:
    output[variant+'_Rate']=[]
    output[variant+'_Lower']=[]
    output[variant+'_Upper']=[]


testable_period=30#5+int(30*stretch)
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

population_of={'02_10':6254603,
               '11_15':3370248,
               '16_24':5950637,
               '25_34':7596145,
               '35_49':10853151,
               '50_69':13618246,
               '70+':7679719
               }

for age in population_of:
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
    df=df[(df['days_since_march1']>start_date)] 
    
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
    
    # calculate theta*I for cases of each variant (theta different for each variant)
    theta_times_I=[]
    for t in range(len(cases)):
        Z_pcr=sum([C[t+j]*P_pcr[j-delta] for j in range(delta,testable_period)])
        Z_lfd=sum([C[t+j]*P_lfd[j] for j in range(testable_period)])
        
        Z=(1-LFD_proportion[t])*Z_pcr+LFD_proportion[t]*Z_lfd
        theta_times_I.append(Z)
    
    
    new_cases={}
    # need to repeaat for each variant
    for variant in variant_proportion:
        temp_list=[]
        # estimate the number of test-positives of each variant type
        for t in ONS_day:
            # j is the time since exposure
            estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
            # now adjust for the variant
            temp_list.append(estimate*variant_proportion[variant][t])
            
        # store in the dictionary
        new_cases[variant]=temp_list


    
    ## loop over the three lines
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    
    #for interval in reporting_multiplier:
    
    for interval in ['Rate','Lower','Upper']:
        new_ONS=[population_of[age]*r/100 for r in rate[interval]]
        ############# END ############                           
    
        #### Calculate th multiplier ###########
        
        params=hill_climb(new_cases,new_ONS)
        
        for variant in params:
            output[variant+'_'+interval].append(params[variant])

        
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
    
    df=df[(df['days_since_march1']>start_date)]
    
    
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
    SGTF_proportion=df['SGTF_proportion'].tolist()
    
    
    ############## COPIED IN ###################
    C=cases+[0 for i in range(100)]

    theta_times_I=[]
    for t in range(len(cases)):
        Z_pcr=sum([C[t+j]*P_pcr[j-delta] for j in range(delta,testable_period)])
        Z_lfd=sum([C[t+j]*P_lfd[j] for j in range(testable_period)])
        
        Z=(1-LFD_proportion[t])*Z_pcr+LFD_proportion[t]*Z_lfd
        theta_times_I.append(Z)

    new_cases={}
    # need to repeaat for each variant
    for variant in variant_proportion:
        temp_list=[]
        # estimate the number of test-positives of each variant type
        for t in ONS_day:
            # j is the time since exposure
            estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
            # now adjust for the variant
            temp_list.append(estimate*variant_proportion[variant][t])
            
        # store in the dictionary
        new_cases[variant]=temp_list

    
    for interval in ['Rate','Lower','Upper']:
        new_ONS=[population_of[region]*r/100 for r in rate[interval]]
        ############# END ############
                            
        #### Calculate the multiplier ###########
        params=hill_climb(new_cases,new_ONS)
        
        for variant in params:
            output[variant+'_'+interval].append(params[variant])

        

pd.DataFrame(output).to_csv('../output/time_independent_rates.csv',index=False)
#pd.DataFrame(output).to_csv('../output/time_independent_rates_(higher_sensitivity).csv',index=False)
#pd.DataFrame(output).to_csv('../output/time_independent_rates_(longer_delay).csv',index=False)
