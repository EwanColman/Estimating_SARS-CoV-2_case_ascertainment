
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt

from scipy.special import erf
import pickle as pk
import scipy.optimize as optimize

def f(x):
    percentage_WT=x[0]
    percentage_VoC=x[1]
    s_max=x[2]
    stretch=x[3]
    
    theta_WT=100/percentage_WT
    theta_VoC=100/percentage_VoC
    # calculate te error as a function of the inputs
    #stretch=round(10/interval,2)
    testable_period=5+int(30*stretch)
    
    # add zeros to make the arrays the same length as the stretched versions
    So=S+[0 for i in range(30,testable_period+1)]
    Ro=R+[0 for i in range(30,testable_period+1)]
    
    # stretch the sensitivity function
    pi=0
    Sv=So[0:pi]
    for tau in range(pi,testable_period):
        f=(tau-pi)/stretch
        #print(tau,pi+int(f)+1,len(So))
        y=So[pi+int(f)]+(f-int(f))*(So[pi+int(f)+1]-So[pi+int(f)])    
        Sv.append(y)
    
    # stretch the incubation period function
    
    Rv=R[0:pi]
    for tau in range(pi,testable_period):
        f=(tau-pi)/stretch
        #print(tau,pi+int(f)+1,len(So))
        y=Ro[pi+int(f)]+(f-int(f))*(Ro[pi+int(f)+1]-Ro[pi+int(f)])    
        Rv.append(y)
           
    Po=[]
    for j in range(delta,testable_period):
            
        numerator=Ro[j-delta]*So[j]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([Ro[i-delta]*So[i] for i in range(testable_period)])
     
        Po.append(numerator/denominator)   
     
    Pv=[]                                                                                                                                          
    for j in range(delta,testable_period):
            
        numerator=Rv[j-delta]*Sv[j]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([Rv[i-delta]*Sv[i] for i in range(testable_period)])
           
        Pv.append(numerator/denominator)


                
    theta_times_I=[]
    theta_times_Iv=[]
    for t in range(len(cases)):
             
        Z=sum([C[t+j]*Po[j-delta] for j in range(delta,testable_period)])
        theta_times_I.append(Z)
        
        Zv=sum([C[t+j]*Pv[j-delta] for j in range(delta,testable_period)])
        theta_times_Iv.append(Zv)
  
    # set a list for normal cases, and an alternative one for variant cases
    new_cases=[]
    nv_cases=[]
    proportion=[]
    # estimate the number of test-positives
    for t in ONS_day:
        # j is the time since exposure
        estimate=sum([theta_times_I[t-j]*So[j] for j in range(testable_period)])
        new_cases.append(estimate)
        # now do the alternative
        #estimate=sum([nv_exposures[t-j]*nv_positive_probability[j] for j in range(testable_period)])
        estimate=sum([theta_times_Iv[t-j]*Sv[j] for j in range(testable_period)])
        nv_cases.append(estimate)
        
        # to use the variant proportion
        proportion.append(variant_proportion[t])
        
    new_ONS=[population_of_England*r/100 for r in rate]

    # calculate the total_error
    e=sum([abs(theta_WT*new_cases[i]*(1-proportion[i])+theta_VoC*nv_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
  
    #print('WT %:',percentage_WT,'VoC %:',percentage_VoC,'s max:',s_max,'interval:',interval,'error =',int(e))

    
    return e


delta=1
#max_day=425

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Mar1='1 March 2021'
Mar1=(datetime.strptime(str(Mar1), '%d %B %Y')-time_zero).days



population_of_England=56287000


### Get the ONS data 
df=pd.read_csv('../raw_data/ONS_incidence.csv')

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

# cut off before delta variant arrives
df=df[df['days_since_march1']<Mar1]

rate=df['Rate'].tolist()
upper=df['Upper'].tolist()
lower=df['Lower'].tolist()
ONS_day=df['days_since_march1'].tolist()
date=df['Date'].tolist()

###### CASES ##########
df=pd.read_csv('../processed_data/England_daily_data.csv')
cases=df['cases'].tolist()
time_in_days=df['Days_since_March1'].tolist()
variant_proportion=df['NV_proportion'].tolist()
#PCR_tests=df['PCR_tests'].tolist()
#LFD_tests=df['LFD_tests'].tolist()

# add a bit of the future to prevent indexing errors
C=cases+[0 for i in range(100)]

# choose th testable period (tail length)

### Lognormal incubation period distribution
# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean,dispersion=5.5,1.52
mu,sigma=np.log(mean),np.log(dispersion)

R=[]
for i in range(1000):
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
S=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]
# in for old variants it should just be 0 after the end


################


reporting_WT=[i for i in range(24,26)]
reporting_VoC=[i for i in range(1,101)]
s_max_values=[0.8]#[i/100 for i in range(1,100)]
stretch_values=[round(1+i/50,2) for i in range(100)]

n=0
error={}
#for delta in range(1,3):
for percentage_WT in reporting_WT:
    for percentage_VoC in reporting_VoC:
        for s_max in s_max_values:
            for stretch in stretch_values:
                
                x=(percentage_WT,percentage_VoC,s_max,stretch)
                error[x]=f(x)
                n=n+1
                #print(x,error[x])
                if n % 1000==0:
                    print(n,'of',len(reporting_WT)*len(reporting_VoC)*len(s_max_values)*len(stretch_values))
                          


                    
pk.dump(error,open('../pickles/error_(test_seeking).p','wb'))
