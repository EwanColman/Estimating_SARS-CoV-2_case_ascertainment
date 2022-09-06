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
    # first convert to daily (with interplotation)
    # daily_theta=[]
    # for k in range(1,len(theta)):
    #     daily_theta=daily_theta+[theta[k-1]+i*(theta[k]-theta[k-1])/7 for i in range(7)]
    

    time=[i*7 for i in range(len(theta))]
    daily_theta=[0 for i in range(time[0])]
    daily_theta=daily_theta+interpolate(time,list(theta))
    while len(daily_theta)<len(cases):
        daily_theta.append(daily_theta[-1])

    #print(len(daily_theta))
    Z=sum([abs(population*rate[t]/100 - 100*sum([theta_times_I[ONS_day[t]-j]*S_pcr[j]/daily_theta[ONS_day[t]-j] for j in range(testable_period)])) for t in range(len(rate))])   

    return Z

population=10000000
max_day=600
delta=1

#### Probability of positive test #################
# add a bit of the future to prevent indexing errors

# choose th testable period (tail length)

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

df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]



#S=S+[0 for i in range(30,testable_period)]

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



#### SURVEILLANCE ########
ons_df=pd.read_csv('../synthetic_data/surveillance.csv',sep=',')
   
# convert data frame to lists for plotting
rate=ons_df.reset_index()

ONS_day=ons_df['Date'].tolist()
    
###### CASES ##########
df=pd.read_csv('../synthetic_data/cases.csv')


cases=df['cases'].tolist()
case_time_in_days=df['Date'].tolist()

   
C=cases+[0 for i in range(100)]
    
theta_times_I=[]
for t in range(len(cases)):
    Y_pcr=sum([C[t+j]*P_pcr[j-delta]*S_pcr[j] for j in range(delta,testable_period)])
    #Y_lfd=sum([C[t+j]*P_lfd[j]*S_lfd[j] for j in range(testable_period)])
    #Y=(1-LFD_proportion[t])*Y_pcr+LFD_proportion[t]*Y_lfd
    theta_times_I.append(Y_pcr)
   
############################################

new_cases=[]

# estimate the number of test-positives
for t in ONS_day:
    # j is the time since exposure
    estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
    
    new_cases.append(estimate)

output_data={}
output_data['reporting_multiplier']=[] 
output_data['incidence']=[]

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




    first_guess=[]
    for j in range(len(new_cases)):
    
        if rate[j]>0:
            first_guess.append(min(100,100*100*new_cases[j]/(population*rate[j]))) 
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
    
    # day estimate is for
    days=[i*7 for i in range(weeks)]
    
    theta=reporting_multiplier
    time=[i*7 for i in range(len(theta))]
    daily_theta=[0 for i in range(time[0])]
    daily_theta=daily_theta+interpolate(time,list(theta))
    while len(daily_theta)<len(cases):
        daily_theta.append(daily_theta[-1])
    
    I=[]

    theta=reporting_multiplier
    time=[i*7 for i in range(len(theta))]
    daily_theta=[0 for i in range(time[0])]
    daily_theta=daily_theta+interpolate(time,list(theta))
    while len(daily_theta)<len(cases):
        daily_theta.append(daily_theta[-1])

    
    I=[0 for i in range(ONS_day[0])]
    I=I+[100*100*theta_times_I[t]/(daily_theta[t]*population) for t in range(ONS_day[0],len(theta_times_I))]

    #print(I)
    #print()
    output_data['reporting_multiplier'].append(reporting_multiplier)
    #output_data['variant_proportion_'+region]=new_variant_proportion
    output_data['incidence'].append(I)

# dates for the region's data  
output_data['date']=days


pk.dump(output_data,open('../pickles/reporting_rates_(synthetic)_200.p','wb'))


