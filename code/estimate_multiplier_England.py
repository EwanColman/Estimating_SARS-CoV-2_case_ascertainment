
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# adjust the sensitivity for new variant with this parameter
#alpha=0.9
delta=2
max_day=492
testable_period=30

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Mar1='1 March 2021'
Mar1=(datetime.strptime(str(Mar1), '%d %B %Y')-time_zero).days

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

df=df[::2]

rate=df['Rate'].tolist()
upper=df['Upper'].tolist()
lower=df['Lower'].tolist()
ONS_day=df['days_since_march1'].tolist()
date=df['Date'].tolist()




print(date[-2:])
# mak the plot and ad the ONS data
plt.figure(figsize=(8.5,4))
#plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
#plt.scatter(ONS_day,rate,color='k',marker='^',label='ONS Surveillance: estimated number of test-positive people')



###### CASES ##########
df=pd.read_csv('../processed_data/England_daily_data.csv')
cases=df['all_cases'].tolist()
time_in_days=df['Days_since_March1'].tolist()
variant_proportion=df['NV_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()
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



e_min=10**10


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
    # now do the alternative
    #estimate=sum([nv_exposures[t-j]*nv_positive_probability[j] for j in range(testable_period)])
    #estimate=sum([theta_times_Iv[t-j]*Sv[j] for j in range(testable_period)])
    #nv_cases.append(estimate)

    # to use the variant proportion
    proportion.append(variant_proportion[t])


population_of_England=56287000
new_ONS=[population_of_England*r/100 for r in rate]

# find the first place where the delta replaces the non-alpha/delta variants
delta_index=0
while ONS_day[delta_index]<Mar1:
    delta_index=delta_index+1                   

plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
plt.scatter(ONS_day,rate,s=10,color='k',marker='^')

#### Calculate th multiplier ###########
# x is the multiplier
x_values=[100/i for i in range(1,101)]
x1_values=[100/i for i in range(5,35)]
x2_values=[100/i for i in range(10,45)]
x3_values=[100/i for i in range(10,65)]
#error=[]
e_min=10**10
for x1 in x1_values:
    for x2 in x2_values:
        for x3 in x3_values:
            #x2=x1 
            #### Calculate an error #######
            #absolute error
            x=[x1 for i in range(delta_index)]+[x3 for i in range(delta_index,len(new_ONS))]
            
            e=sum([abs(x[i]*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
            # squared error
            #e=(sum([(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i])**2 for i in range(len(new_ONS))]))**(1/2)
            #print(int(100/x1),int(100/x2),np.log(e))
            #error.append(e)
            if e<e_min:
                e_min=e
                best_x=(x1,x2,x3)
                
                best_new_cases=new_cases.copy()
            #best_nv_cases=nv_cases.copy()

# take the value with the smallest error
x1,x2,x3=best_x
#adjusted_cases=[x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i] for i in range(len(new_ONS))]
#adjusted_cases=[x1*best_new_cases[i]*(1-proportion[i])+x2*best_new_cases[i]*proportion[i] for i in range(len(new_ONS))]

x=[x1 for i in range(delta_index)]+[x3 for i in range(delta_index,len(new_ONS))]
            
adjusted_cases=[x[i]*best_new_cases[i]*(1-proportion[i])+x2*best_new_cases[i]*proportion[i] for i in range(len(new_ONS))]

plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
plt.scatter(ONS_day,rate,color='k',marker='^',label='ONS CIS Surveillance: estimated number of test-positive people')
plt.scatter(ONS_day,[100*c/population_of_England for c in adjusted_cases],s=10,facecolors='none', edgecolors='b',label='$\\theta_{o}='+str(round(100/x1))+'$%, $\\theta_{\\alpha}='+str(round(100/x2))+'$%, $\\theta_{\\delta}='+str(round(100/x3))+'$%')
# add the result to the plot
#plt.text(100,3,



dates=['01/0'+str(i)+'/2020' for i in range(6,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,6)]
day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]


plt.yticks([i/2 for i in range(5)])
plt.ylim([0,3.1])
plt.ylabel('Percentage of people')
#plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,dates_words,rotation=0)
#plt.xlim([60,260])

#plt.xlabel('Date')
plt.legend(loc=2)


plt.xlim([60,max_day])


#plt.savefig('../figures/England_different_multipliers.png',format='png',dpi=300,bbox_inches='tight')
