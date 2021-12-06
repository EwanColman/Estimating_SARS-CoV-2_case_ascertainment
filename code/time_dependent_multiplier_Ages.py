import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf
from scipy import stats
import pickle as pk
# for correlations use the last T points


def get_incidence(ONS_day,theta_times_I,reporting_multiplier):
    #### Interpolation
    
    # et te first day that is not None
    
    start_day=ONS_day[0]
    # linear interpolation to get all days
    daily_theta=[6 for i in range(start_day)]
    
    last_t=start_day
    last_v=0
    
    theta=reporting_multiplier[0]
    if theta==None:
        theta=0

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
    # add the last bit
    while len(daily_theta)<len(cases)+9:
        daily_theta.append(theta)    
    
    # compute incidence
    I=[]
    for i in range(len(theta_times_I)):
        I.append(100*theta_times_I[i]/daily_theta[i+9])
    return I



# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
variant_proportion=df['NV_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()

folder='../raw_data/Surveillance'
max_day=530
delta=1
fs=15
#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Sep1='1 September 2020'
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days

# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
population_of={'02_10':6264662,#7527576,
               '11_15':2664513,
               '16_24':5860347,
               '25_34':7567108,
               '35_49':10864046,
               '50_69':13688509,
               '70+':8114862
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
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]

#S=S+[0 for i in range(30,testable_period)]
df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]


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


# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,6)]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
   

# start making the latex table
table='\\begin{tabular}{l|cccc} \n \\toprule \n & PCR tests & LF Tests & VOC proportion & Cycle threshold \\\ \n \\midrule \n'



fig=plt.figure(figsize=(12,6))
plt.subplots_adjust(hspace=0,wspace=0)

#region='London'

output_data={}

m=0
i=0
print()
for age in population_of:

    
    i=i+1
#    n=i % 3
#    m=1+int(n/3)
    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,85,'Age '+name_of[age])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,100])
    if i%3==1:
        ax.set_yticks([0,20,40,60])
    else:
        ax.set_yticks([])
        
    if i>4:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i==4: 
        plt.ylabel('Percentage of infections reported')
    # if i==8:
    #     plt.xlabel('Date')

    # get te data for the region   
    
    #### SURVEILLANCE ########
    df=pd.read_csv(folder+'/'+age+'_surveillance.csv',sep=',')
    print(df.head())
    
    # add a new column with header days since Mar1
    time_in_days=[]
    date=df['Date'].tolist()
    while date:
        d=date.pop(0)
        #day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days

        
        time_in_days.append(day_numerical)
    
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1',ascending=True)
    #df.drop(df.tail(1).index,inplace=True)
    df=df[df['days_since_march1']>100] 
    # convert data frame to lists for plotting
    rate=df.reset_index()
    #rate=df['Rate'].tolist()
    #upper=df['Upper'].tolist()
    #lower=df['Lower'].tolist()
    ONS_day=df['days_since_march1'].tolist()
    
    # re-read it since it was popped empty
    date=df['Date'].tolist()
    #dates_simce_Nov1=dates_simce_Nov1[:-2]
    #ONS_day=[day for day in ONS_day if day>=Sep1]
    #rate=rate[-len(ONS_day):]

    
    ###### CASES ##########
    df=pd.read_csv('../raw_data/Diagnostic/'+age+'_cases.csv')
    print(df.head())
    print()
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
    
    
    
    #time_in_days=df['Days_since_March1'].tolist()
    
    C=cases+[0 for i in range(100)]
    
    theta_times_I=[]
    for t in range(len(cases)):
        Z_pcr=sum([C[t+j]*P_pcr[j-delta] for j in range(delta,testable_period)])
        Z_lfd=sum([C[t+j]*P_lfd[j] for j in range(testable_period)])
        
        Z=(1-LFD_proportion[t])*Z_pcr+LFD_proportion[t]*Z_lfd
        theta_times_I.append(Z)
    
    psi_times_I=[]
    for t in range(len(cases)):
        Y_pcr=sum([C[t+j]*P_pcr[j-delta]*S_pcr[j] for j in range(delta,testable_period)])
        Y_lfd=sum([C[t+j]*P_lfd[j]*S_lfd[j] for j in range(testable_period)])
        Y=(1-LFD_proportion[t])*Y_pcr+LFD_proportion[t]*Y_lfd
        psi_times_I.append(Y)
    
    
    
    new_cases=[]
    new_cases_ts=[]
    new_variant_proportion=[]
    # estimate the number of test-positives
    for t in ONS_day:#range(20,len(cases)):
        # j is the time since exposure
        estimate=sum([psi_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
        #cases_on_day[surveillance_day[t]]=estimate
        new_cases.append(estimate)
        
        # do it again for test-seeking
        estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
        #cases_on_day[surveillance_day[t]]=estimate
        new_cases_ts.append(estimate)
        new_variant_proportion.append(variant_proportion[t])
    
    #### calculate multiplier for every time point ##########
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    #testing_multiplier=[]
    
    for j in range(len(new_cases)):
        for interval in reporting_multiplier:
            if rate[interval][j]>0:
                reporting_multiplier[interval].append(min(100,100*100*new_cases[j]/(population_of[age]*rate[interval][j])))
    #        testing_multiplier.append(100*100*new_cases_ts[j]/(population_of[age]*rate[j]))
     
            else:
                reporting_multiplier[interval].append(100)
    #        testing_multiplier.append(None)
    
    
    plt.scatter(ONS_day,reporting_multiplier['Rate'],s=10,facecolors='#006666')
    plt.plot(ONS_day,reporting_multiplier['Lower'])
    plt.plot(ONS_day,reporting_multiplier['Upper'])
    
    I=get_incidence(ONS_day,psi_times_I,reporting_multiplier['Rate'])
        
    plt.fill_between(range(len(cases)),[0 for i in range(len(cases))],[100*100*i/population_of[age] for i in I],color='k',linewidth=0,alpha=0.2)

    
    #plt.scatter(ONS_day,testing_multiplier,s=10,facecolors='#006666')
    #output_data[age]=multiplier
    output_data['reporting_multiplier_'+age]=reporting_multiplier
    #output_data['testing_multiplier_'+age]=testing_multiplier
    output_data['incidence_'+age]=[100*i/population_of[age] for i in I]
    output_data['variant_proportion_'+age]=new_variant_proportion
    output_data['sgtf_proportion_'+age]=variant_proportion
    
    # find the daily incidence 
    

output_data['date']=ONS_day
pk.dump(output_data,open('../pickles/reporting_rates_(age).p','wb'))
   
plt.savefig('../figures/time_dependent_multiplier_ages.png',format='png', bbox_inches='tight',dpi=256)



