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
    Z=sum([abs(population_of[age]*rate[interval][t]/100 - 100*sum([theta_times_I[ONS_day[t]-j]*((1-vax_proportion[t])*S_pcr[j]+vax_proportion[t]*S_pcr_vax[j])/daily_theta[ONS_day[t]-j] for j in range(testable_period)])) for t in range(len(rate[interval]))])   

    return Z



# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
SGTF_proportion=df['SGTF_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()

folder='../raw_data/Surveillance'
max_day=825
delta=1
fs=15
trunc=10
efficacy=0.56


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

# truncate for vaccinated
S_pcr_vax=S_pcr[0:trunc]+[0 for i in range(trunc,testable_period)]
S_lfd_vax=S_lfd[0:trunc]+[0 for i in range(trunc,testable_period)]

# make the probability function as a list
P_pcr=[]
P_pcr_vax=[]

for j in range(delta,testable_period):
    numerator=R[j-delta]#*S[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i-delta]*S_pcr[i] for i in range(testable_period)])
    P_pcr.append(numerator/denominator)
    
    denominator=sum([R[i-delta]*S_pcr_vax[i] for i in range(testable_period)])
    P_pcr_vax.append(numerator/denominator)  
    

# make the probability function as a list
P_lfd=[]
P_lfd_vax=[]
for j in range(testable_period):
    numerator=R[j]#*S[j]
    # sum([C[t+k] for k in range(3,10)])*
    denominator=sum([R[i]*S_lfd[i] for i in range(testable_period)])
    P_lfd.append(numerator/denominator)

    denominator=sum([R[i]*S_lfd_vax[i] for i in range(testable_period)])
    P_lfd_vax.append(numerator/denominator) 

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
    print(age)
    
    i=i+1
#    n=i % 3
#    m=1+int(n/3)
    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,4,'Age '+name_of[age])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,5])
    #if i%3==1:
    #    ax.set_yticks([0,20,40,60])
    #else:
    #    ax.set_yticks([])
        
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
    
    #### vaccination ####   
    df=pd.read_csv('../raw_data/Vaccination/'+age+'_vaccinations.csv')
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
    
    # add 0s from march 1
    vaxed=[0 for i in range(min(df['days_since_march1']))]
    vaxed=vaxed+df['vaccinated'].tolist()
    
    # calculate as a proportion of pop
    vax_pop=[sum(vaxed[:i])/population_of[age] for i in range(len(vaxed))]
    # fill to the end
    vax_pop=vax_pop+[vax_pop[-1] for i in range(len(vax_pop),len(cases))] 
    ### convert to proportion of infections ###
    vax_proportion=[(1-efficacy)*v/(1-efficacy*v) for v in vax_pop]
    
    # add a bit extra to avoid index errors
    C=cases+[0 for i in range(100)]
     
    theta_times_I=[]
    for t in range(len(cases)):
      
        Y_pcr=sum([C[t+j]*P_pcr[j-delta]*S_pcr[j] for j in range(delta,testable_period)])
        Y_lfd=sum([C[t+j]*P_lfd[j]*S_lfd[j] for j in range(testable_period)])
        
        Y_pcr_vax=sum([C[t+j]*P_pcr_vax[j-delta]*S_pcr_vax[j] for j in range(delta,testable_period)])
        Y_lfd_vax=sum([C[t+j]*P_lfd_vax[j]*S_lfd_vax[j] for j in range(testable_period)])
        
        Y_unvax=(1-LFD_proportion[t])*Y_pcr+LFD_proportion[t]*Y_lfd
        Y_vax=(1-LFD_proportion[t])*Y_pcr_vax+LFD_proportion[t]*Y_lfd_vax
        Y=(1-vax_proportion[t])*Y_unvax+vax_proportion[t]*Y_vax
        theta_times_I.append(Y)
    
    
    new_cases=[]

    # estimate the number of test-positives
    for t in ONS_day:#range(20,len(cases)):
       
        # do it again for test-seeking
        estimate_unvax=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        estimate_vax=sum([theta_times_I[t-j]*S_pcr_vax[j] for j in range(testable_period)])
        estimate=(1-vax_proportion[t])*estimate_unvax+vax_proportion[t]*estimate_vax
        
        #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
        #cases_on_day[surveillance_day[t]]=estimate
        new_cases.append(estimate)
        
    
    #### calculate multiplier for every time point ##########
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    #testing_multiplier=[]
      
    first_guess={'Lower':[],'Upper':[],'Rate':[]}
    for j in range(len(new_cases)):
        for interval in first_guess:
            if rate[interval][j]>0:
                first_guess[interval].append(min(100,100*100*new_cases[j]/(population_of[age]*rate[interval][j]))) 
            else:
                first_guess[interval].append(100)
    
    weeks=int(len(cases)/7)
    #### calculate multiplier for every time point ##########
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    
    for interval in reporting_multiplier:
        # initial guess has zeros at the start followed by the interpolated first estimate
        initial=[0 for i in range(ONS_day[0])]+interpolate(ONS_day,first_guess[interval])
        # add the final value enough times to make it up-to-date with case data + the time shift
        initial=initial+[initial[-1] for i in range(len(initial),len(cases)+9)]
        # take weekly time points, shifted 9 days to map to date of exposure (approximately)
        initial=[initial[7*i+9] for i in range(weeks)]
        
        # optimization using scipy, define bounds 
        bnds = tuple([(0, None) for i in range(weeks)])
        # do the optimization (can try different methods)
        res = optimize.minimize(f,initial, method='COBYLA', bounds=bnds, options={'disp': False})
        # extract the result
        reporting_multiplier[interval] = [r for r in res['x']] 
    
    # and plot    
    days=[i*7 for i in range(weeks)]
    
    #print(first_guess['Rate'])
    #plt.plot(first_guess['Rate'])
    #plt.scatter(days,initial,s=10,vmin=0,vmax=1)
    #plt.scatter(days,reporting_multiplier['Rate'],s=10,vmin=0,vmax=1)
    
    
    theta=reporting_multiplier['Rate']
    time=[i*7 for i in range(len(theta))]
    daily_theta=[0 for i in range(time[0])]
    daily_theta=daily_theta+interpolate(time,list(theta))
    while len(daily_theta)<len(cases):
        daily_theta.append(daily_theta[-1])
    
    
    print([100*sum([theta_times_I[ONS_day[t]-j]*S_pcr[j]/daily_theta[ONS_day[t]-j] for j in range(testable_period)]) for t in range(len(rate['Rate']))][-5:])
    print([population_of[age]*rate['Rate'][t]/100 for t in range(len(rate['Rate']))][-5:])

    plt.plot(ONS_day,[100*100*sum([theta_times_I[t-j]*((1-vax_proportion[t])*S_pcr[j]+vax_proportion[t]*S_pcr_vax[j])/(daily_theta[t-j]*population_of[age]) for j in range(testable_period)]) for t in ONS_day])
    
    
    plt.scatter(ONS_day,rate['Rate'])

    print()
    I=[0 for i in range(ONS_day[0])]
    I=I+[100*100*theta_times_I[t]/(daily_theta[t]*population_of[age]) for t in range(ONS_day[0],len(theta_times_I))]
    
    '''
    
    plt.scatter(ONS_day,reporting_multiplier['Rate'],s=10,facecolors='#006666')
    plt.plot(ONS_day,reporting_multiplier['Lower'])
    plt.plot(ONS_day,reporting_multiplier['Upper'])
    
    I=get_incidence(ONS_day,psi_times_I,reporting_multiplier['Rate'])
        
    plt.fill_between(range(len(cases)),[0 for i in range(len(cases))],[100*100*i/population_of[age] for i in I],color='k',linewidth=0,alpha=0.2)

    '''
    
    #plt.scatter(ONS_day,testing_multiplier,s=10,facecolors='#006666')
    #output_data[age]=multiplier
    output_data['reporting_multiplier_'+age]=reporting_multiplier
    #output_data['testing_multiplier_'+age]=testing_multiplier
    output_data['incidence_'+age]=I
    output_data['sgtf_proportion_'+age]=SGTF_proportion


output_data['date']=days

pk.dump(output_data,open('../pickles/reporting_rates_(age)_VAX.p','wb'))
   
plt.savefig('../figures/time_dependent_multiplier_ages.png',format='png', bbox_inches='tight',dpi=256)



