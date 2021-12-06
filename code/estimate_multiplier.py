
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# adjust the sensitivity for new variant with this parameter
#alpha=0.9
delta=2
max_day=492
testable_period=60

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

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


# adjust according to OR drop out 
#df=pd.read_csv('processed_data/England_daily_data.csv')
#OR_only=df['OR only'].tolist()
#lower=[lower[i]*(1-(OR_only[ONS_day[i]]/100)) for i in range(len(ONS_day))]
#upper=[upper[i]*(1-(OR_only[ONS_day[i]]/100)) for i in range(len(ONS_day))]
#rate=[rate[i]*(1-(OR_only[ONS_day[i]]/100)) for i in range(len(ONS_day))]
   


print(date[-2:])
# mak the plot and ad the ONS data
plt.figure(figsize=(8.5,4))
plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
plt.scatter(ONS_day,rate,color='k',marker='^',label='ONS Surveillance: estimated number of test-positive people')



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
# in for old variants it should just be 0 after the end


peak=3
s_max=0.61

e_min=10**10
#for s_max in [0.1*i for i in range(4,10)]:
#for peak in range(10):
#for testable_period in range(40,50):

S=S+[0 for i in range(30,testable_period+1)]
So=S+[0 for i in range(30,testable_period+1)]
Ro=R+[0 for i in range(30,testable_period+1)]

# for the variant, choose the peak and the tail value (same as)
#Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,120,10)]+[np.mean(df['median'].tolist()[120:130]) for i in range(120,300,10)]

#Sv=[i*s_max/peak for i in range(peak)]+[s_max*(testable_period-i)/(testable_period-peak) for i in range(peak,testable_period)]
#

#stretch=1.6

interval=5
print('stretch=',10/interval)
#testable_period=int(30*stretch)
#Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,40,10)]
#Sv=Sv+[np.mean(df['median'].tolist()[i:i+interval]) for i in range(40,300,interval)]
#Sv=Sv+[0 for i in range(len(Sv),testable_period)]
#Sv=[s*s_max/0.785 for s in Sv]


stretch=2.1
pi=4
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


#print()

# the time lag is the number of 0s in the first list
#R=[incubation_probability[i] for i in range(testable_period-1)]

# create the combined probability distribution (distribution of time from exposure to test)
#exp_to_positive_probability=[exposure_to_test_probability[i]*positive_probability[i] for i in range(testable_period)]


Po=[]
for j in range(delta,testable_period):
        
    numerator=Ro[j-delta]#*So[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([Ro[i-delta]*So[i] for i in range(testable_period)])
 
    Po.append(numerator/denominator)   
 
Pv=[]                                                                                                                                          
for j in range(delta,testable_period):
        
    numerator=Rv[j-delta]#*Sv[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([Rv[i-delta]*Sv[i] for i in range(testable_period)])
       
    Pv.append(numerator/denominator)
'''
theta_times_I=[]
theta_times_Iv=[]
for t in range(len(cases)):
         
    Z=sum([C[t+j]*Po[j-delta] for j in range(delta,testable_period)])
    theta_times_I.append(Z)
    
    Zv=sum([C[t+j]*Pv[j-delta] for j in range(delta,testable_period)])
    theta_times_Iv.append(Zv)

'''

theta_times_I=[]
theta_times_Iv=[]
for t in range(len(cases)):
    
    Z=0
    
    for j in range(delta,testable_period):
        
        #numerator=R[j-delta]#*S[j]
        # sum([C[t+k] for k in range(3,10)])*
        #denominator=sum([R[i-delta]*So[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        Z=Z+C[t+j]*Po[j-delta]#(numerator/denominator)
        
        #print(numerator/denominator,Po[j-delta])
    
    theta_times_I.append(Z)
    
    Zv=0
    
    for j in range(delta,testable_period):
        
        # numerator=R[j-delta]#*Sv[j]
        # sum([C[t+k] for k in range(3,10)])*
        # denominator=sum([R[i-delta]*Sv[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        
        
        Zv=Zv+C[t+j]*Pv[j-delta]#(numerator/denominator)
    
    theta_times_Iv.append(Zv)
    


# set a list for normal cases, and an alternative one for variant cases
new_cases=[]
nv_cases=[]
proportion=[]
# estimate the number of test-positives
for t in ONS_day:
    # j is the time since exposure
    estimate=sum([theta_times_I[t-j]*S[j] for j in range(testable_period)])
    new_cases.append(estimate)
    # now do the alternative
    #estimate=sum([nv_exposures[t-j]*nv_positive_probability[j] for j in range(testable_period)])
    estimate=sum([theta_times_Iv[t-j]*Sv[j] for j in range(testable_period)])
    nv_cases.append(estimate)
    
    # to use the variant proportion
    proportion.append(variant_proportion[t])
    #proportion.append(0)
    # to use the test proportion 
    #proportion.append(sum(LFD_tests[t-4:t+3])/(sum(LFD_tests[t-4:t+3])+PCR_tests[t+3]))
    



population_of_England=56287000
new_ONS=[population_of_England*r/100 for r in rate]


#### Calculate th multiplier ###########
# x is the multiplier
x_values=[100/i for i in range(1,100)]
#error=[]
#e_min=10**10
#for x1 in x_values:
x1=4
for x1 in x_values:
    for x2 in x_values:
        
        # calculate the total_error
        e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
        #error.append(e)
        #print(e)
        if e<e_min:
            e_min=e
            best_x=(x1,x2,peak,testable_period,s_max)
            
            best_new_cases=new_cases.copy()
            best_nv_cases=nv_cases.copy()   
# take the value with the smallest error
#x=x_values[error.index(min(error))]

print(best_x,e_min)
x1,x2,peak,testable_period,s_max=best_x
 
adjusted_cases=[x1*best_new_cases[i]*(1-proportion[i])+x2*best_nv_cases[i]*proportion[i] for i in range(len(new_ONS))]

plt.scatter(ONS_day,[100*c/population_of_England for c in adjusted_cases],facecolors='none', edgecolors='b',label='Estimated assuming '+str(round(100/x2))+'% of VOC infections reported, and '+str(round(100/x1))+'% for other variants')



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


day=(datetime.strptime('29/12/2020', '%d/%m/%Y')-time_zero).days

plt.xlim([60,max_day])


plt.savefig('../figures/England_different_multipliers_and_sensitivities.png',format='png',dpi=300,bbox_inches='tight')

# for i in range(len(lower)):
#     if lower[i]<[100*c/population_of_England for c in adjusted_cases][i] and upper[i]>[100*c/population_of_England for c in adjusted_cases][i]:
#         x='In'
#     else:
#         x='Out'
        
#     #print(x,lower[i],[100*c/population_of_England for c in adjusted_cases][i],upper[i])
#     print(abs(rate[i]-[100*c/population_of_England for c in adjusted_cases][i]))
