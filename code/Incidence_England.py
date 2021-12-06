import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy import stats
import pickle as pk
from datetime import timedelta
# 

cmap = plt.cm.get_cmap('winter')
max_day=506
min_day=100
delta=2
population_of_England=56287000
#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

### Lognormal incubation period distribution
# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean,dispersion=5.5,1.52
mu,sigma=np.log(mean),np.log(dispersion)
testable_period=60

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
S=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]

S=S+[0 for i in range(30,testable_period)]


# the value in position i is th probability of testing positive i days after infection




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
ONS_day=df['days_since_march1'].tolist()
date=df['Date'].tolist()


df=pd.read_csv('../processed_data/England_daily_data.csv')
print(len(df))
cases=df['cases'].tolist()
time_in_days=df['Days_since_March1'].tolist()
variant_proportion=df['NV_proportion'].tolist()
PCR_tests=df['PCR_tests'].tolist()
LFD_tests=df['LFD_tests'].tolist()
#CT_mean=df['mean_CT'].tolist()
#OR_only=df['OR_only'].tolist()

# make the fnction and add some 0s to the end
C=cases+[0 for i in range(100)]

theta_times_I=[]
psi_times_I=[]
for t in range(len(cases)):
    # Z is the total sum
    Z=0
    Y=0
    for j in range(delta,testable_period):      
        numerator=R[j-delta]*S[j]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([R[i-delta]*S[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        Z=Z+C[t+j]*(numerator/denominator)
        
        # repeat for test-seeking
        numerator=R[j-delta]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([R[i-delta]*S[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        Y=Y+C[t+j]*(numerator/denominator)

    theta_times_I.append(Z)
    psi_times_I.append(Y)
 

new_cases=[]
new_cases_ts=[]
new_variant_proportion=[]
# estimate the number of test-positives
for t in ONS_day:#range(20,len(cases)):
    # j is the time since exposure
    estimate=sum([theta_times_I[t-j]*S[j] for j in range(testable_period)])
    #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
    #cases_on_day[surveillance_day[t]]=estimate
    new_cases.append(estimate)
    
    # repeat for test-seeking
    estimate=sum([psi_times_I[t-j]*S[j] for j in range(testable_period)])
    #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
    #cases_on_day[surveillance_day[t]]=estimate
    new_cases_ts.append(estimate)
    
    new_variant_proportion.append(variant_proportion[t])
 
#### calculate multiplier for every time point ##########
reporting_multiplier=[]
testing_multiplier=[]

for j in range(len(new_cases)):
    if rate[j]>0:
        reporting_multiplier.append(100*100*new_cases[j]/(population_of_England*rate[j]))
        testing_multiplier.append(100*100*new_cases_ts[j]/(population_of_England*rate[j]))
    
    else:
        reporting_multiplier.append(None)
        testing_multiplier.append(None)

# interpolate
start_day=ONS_day[0]
# linear interpolation to get all days
daily_theta=[6 for i in range(start_day)]
last_t=start_day
last_v=0
theta=reporting_multiplier[0]
time_in_days=ONS_day.copy()
r_m=reporting_multiplier.copy()

while time_in_days:
    
    
    t=time_in_days.pop(0)

    v=r_m.pop(0)
    if v:
        delta_t=t-last_t
        delta_v=v-last_v
        #print(delta_p,delta_t)
        for i in range(delta_t):
            theta=theta+delta_v/delta_t
            
            #print(delta_p/delta_t)
            daily_theta.append(theta)
        last_t=t
        last_v=v

while len(daily_theta)<len(cases):
    daily_theta.append(theta)    

fig=plt.figure(figsize=(12,6))
plt.plot(daily_theta)

I=[]
for i in range(len(theta_times_I)):
    I.append(100*theta_times_I[i]/daily_theta[i])


# and plot   
fig=plt.figure(figsize=(12,6))
 
#plt.plot(daily_theta)
plt.plot(I,label='Estimated daily new infections')
plt.plot([c*4 for c in cases],label='Reported cases $\\times 4$')
plt.legend(loc=2)
#plt.scatter(ONS_day,reporting_multiplier,c=new_variant_proportion,s=10,cmap='winter',vmin=0,vmax=1)
 
dates=['01/0'+str(i)+'/2020' for i in range(6,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,7)]
day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]

daily_dates=[datetime.date((time_zero+timedelta(t))) for t in range(len(I))]
pd.DataFrame({'Date':daily_dates,'Infectins':I}).to_csv('Incidence_estimate_England.csv',index=False)

pk.dump(I,open('../pickles/England_incidence.p','wb'))

#plt.yticks([i/2 for i in range(5)])
#plt.ylim([0,3.1])
#plt.ylabel('Percentage of people')
#plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,dates_words,rotation=0)

#plt.legend(loc=2)
#plt.yscale('log')
plt.savefig('../admissions_figures/cases_and_estimated_new_infections.png',format='png',dpi=300,bbox_inches='tight')



df=pd.read_csv('../raw_data/admissions_England.csv')

time_in_days=[]
date=df['date'].tolist()
while date:
    d=date.pop(0)
    #print(d)
    day=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
    # minus 4 to make it a mid-week estimate
    time_in_days.append(day)
df['days_since_march1']=time_in_days

df=df.sort_values('days_since_march1',ascending=True)
df=df[df['days_since_march1']>=0]
time_in_days=df['days_since_march1'].tolist()
admissions=[0 for i in range(time_in_days[0])]+df['newAdmissions'].tolist()

#print(len(cases),len(admissions))
### CASES #######
fig=plt.figure(figsize=(8,3))
plt.xticks(day_numerical,dates_words,rotation=0)
 
for n in [7]:

    ratio=[np.mean([admissions[i-j]/cases[i-j-n] for j in range(7)]) for i in range(len(admissions))]

    plt.plot(ratio,label='n='+str(n))

plt.title('Number of hospital admissions divided by number of cases 7 days earlier \n (7-day moving average)')
plt.xlim([min_day,max_day])
plt.ylim([0,0.25])
#plt.legend(loc=2)
plt.savefig('../admissions_figures/case_admission_ratio.png',format='png',dpi=300,bbox_inches='tight')

# and plot   

#### Estimated infections ############
fig=plt.figure(figsize=(12,6))
plt.xticks(day_numerical,dates_words,rotation=0)
 
for n in [10,12,14,16,18]:

    ratio=[np.mean([admissions[i-j]/I[i-j-n] for j in range(7)]) for i in range(len(admissions))]

    plt.plot(ratio[:-6],label='n='+str(n))

plt.title('new admissions / new infections n days earlier (7-day average after calculating ratio)')
plt.xlim([min_day,max_day])
plt.ylim([0,0.1])
plt.legend(loc=2)
plt.savefig('../admissions_figures/infection_admission_ratio.png',format='png',dpi=300,bbox_inches='tight')

#### Admissions ############
fig=plt.figure(figsize=(12,6))
plt.xticks(day_numerical,dates_words,rotation=0)
 

plt.plot(admissions)

#plt.title('new admissions / new infections n days earlier (7-day average after calculating ratio)')
plt.xlim([min_day,max_day])
plt.yscale('log')
#plt.ylim([0,0.1])
#plt.legend(loc=2)
#plt.savefig('../admissions_figures/infection_admission_ratio.png',format='png',dpi=300,bbox_inches='tight')



print(sum(I),100*sum(I)/population_of_England)
