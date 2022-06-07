import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy import stats
import pickle as pk
# for correlations use the last T points

import scipy.optimize as optimize

def interpolate(X_series,Y_series):
    
    x_series=X_series.copy()
    y_series=Y_series.copy()
    
    #### Interpolation
    # et te first day that is not None
    first_x=x_series[0]
    # linear interpolation to get all days
    y_interpolated=[]
    
    last_x=first_x
    last_y=0
    
    y_int=y_series[0]
    
    while x_series:
        x=x_series.pop(0)
        y=y_series.pop(0)
        if y:
            delta_x=x-last_x
            delta_y=y-last_y
            #print(delta_p,delta_t)
            for i in range(delta_x):
                y_int=y_int+delta_y/delta_x        
                #print(delta_p/delta_t)
                y_interpolated.append(y_int)
                
            last_x=x
            last_y=y  
   
    return y_interpolated


def f(theta):
   
    time=[i*7 for i in range(len(theta))]
    daily_theta=[0 for i in range(time[0])]
    daily_theta=daily_theta+interpolate(time,list(theta))
    while len(daily_theta)<len(cases):
        daily_theta.append(daily_theta[-1])

    #print(len(daily_theta))
    Z=sum([abs(population_of[age]*rate[t]/100 - 100*sum([theta_times_I[ONS_day[t]-j]*S_pcr[j]/daily_theta[ONS_day[t]-j] for j in range(testable_period)])) for t in range(len(rate))])   

    return Z



# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
SGTF_proportion=df['SGTF_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()

folder='../raw_data/Surveillance'
max_day=825
delta=1
fs=15
#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Sep1='1 September 2020'
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days

# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
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
               

#### Probability of positive test #################
# the value in position i is th probability of testing positive i days after infection
### Lognormal incubation period distribution
# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean,dispersion=5.5,1.52
mu,sigma=np.log(mean),np.log(dispersion)
testable_period=30

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
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]

#S=S+[0 for i in range(30,testable_period)]
df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]


# make the probability function as a list
P_pcr=[]
for j in range(delta,testable_period):
    numerator=R[j-delta]#*S[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i-delta]*S_pcr[i] for i in range(testable_period)])
    P_pcr.append(numerator/denominator)   

# make the probability function as a list
P_lfd=[]
for j in range(testable_period):
    numerator=R[j]#*S[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i]*S_lfd[i] for i in range(testable_period)])
    P_lfd.append(numerator/denominator) 


output_data={}

for age in population_of:
    print(age)
    output_data['reporting_multiplier_'+age]=[] 
    output_data['incidence_'+age]=[]
    
    #### SURVEILLANCE ########
    ons_df=pd.read_csv(folder+'/'+age+'_surveillance.csv',sep=',')
   
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in ons_df['Date'].tolist()]
    
    # add it to the dataframe
    ons_df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    ons_df=ons_df.sort_values('days_since_march1',ascending=True)
    #df.drop(df.tail(1).index,inplace=True)
    ons_df=ons_df[ons_df['days_since_march1']>100]

    ONS_day=ons_df['days_since_march1'].tolist() 
    
    ###### CASES ##########
    df=pd.read_csv('../raw_data/Diagnostic/'+age+'_cases.csv')
    
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days for d in df['date'].tolist()]

    df['days_since_march1']=time_in_days
    df=df.sort_values('days_since_march1',ascending=True)
    df=df[df['days_since_march1']>=0]
    cases=df['cases'].tolist()
    case_time_in_days=df['days_since_march1'].tolist()
    dates=df['date'].tolist()
     
    
    C=cases+[0 for i in range(100)]
      
    theta_times_I=[]
    for t in range(len(cases)):
        Y_pcr=sum([C[t+j]*P_pcr[j-delta]*S_pcr[j] for j in range(delta,testable_period)])
        Y_lfd=sum([C[t+j]*P_lfd[j]*S_lfd[j] for j in range(testable_period)])
        Y=(1-LFD_proportion[t])*Y_pcr+LFD_proportion[t]*Y_lfd
        theta_times_I.append(Y)
    
    new_cases=[]

    # estimate the number of test-positives
    for t in ONS_day:#range(20,len(cases)):
       
        # do it again for test-seeking
        estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
        #cases_on_day[surveillance_day[t]]=estimate
        new_cases.append(estimate)
    
    
    # 200 resamplings
    for i in range(200):
        print(i)
        # the resampled rate
        rate=[]
        
        for j,row in ons_df.iterrows():
            # calculate mean and sd for that day
            #print(row)
            mu=row['Rate']
            sigma=(row['Upper']-row['Lower'])/(2*1.96)
            # sample from distribution for new rate
            rate.append(max(0,np.random.normal(mu, sigma)))

        #### calculate multiplier for every time point ##########
        reporting_multiplier=[]
          
        first_guess=[]
        for j in range(len(new_cases)):
        
            if rate[j]>0:
                first_guess.append(min(100,100*100*new_cases[j]/(population_of[age]*rate[j]))) 
            else:
                first_guess.append(100)
        
        weeks=int(len(cases)/7)
        #### calculate multiplier for every time point ##########
        reporting_multiplier=[]
        
    
        # initial guess has zeros at the start followed by the interpolated first estimate
        initial=[0 for i in range(ONS_day[0])]+interpolate(ONS_day,first_guess)
        # add the final value enough times to make it up-to-date with case data + the time shift
        initial=initial+[initial[-1] for i in range(len(initial),len(cases)+9)]
        # take weekly time points, shifted 9 days to map to date of exposure (approximately)
        initial=[initial[7*i+9] for i in range(weeks)]
        
        # optimization using scipy, define bounds 
        bnds = tuple([(0, None) for i in range(weeks)])
        # do the optimization (can try different methods)
        res = optimize.minimize(f,initial, method='COBYLA', bounds=bnds, options={'disp': False})
        # extract the result
        reporting_multiplier = [r for r in res['x']] 
    
        # and plot    
        days=[i*7 for i in range(weeks)]
        
        theta=reporting_multiplier
        time=[i*7 for i in range(len(theta))]
        daily_theta=[0 for i in range(time[0])]
        daily_theta=daily_theta+interpolate(time,list(theta))
        while len(daily_theta)<len(cases):
            daily_theta.append(daily_theta[-1])
    
        print()
        I=[0 for i in range(ONS_day[0])]
        I=I+[100*100*theta_times_I[t]/(daily_theta[t]*population_of[age]) for t in range(ONS_day[0],len(theta_times_I))]
    
    
        output_data['reporting_multiplier_'+age].append(reporting_multiplier)
        output_data['incidence_'+age].append(I)
  
    # find the daily incidence 




output_data['date']=days
pk.dump(output_data,open('../pickles/reporting_rates_(age)_200.p','wb'))


