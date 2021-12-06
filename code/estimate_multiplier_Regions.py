import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

folder='../raw_data/Regions'
max_day=450

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
population_of={'NorthWest':7341196,
               'YorkshireandTheHumber':5502967,
               'NorthEast':2669941,
               'WestMidlands':5934037,
               'EastMidlands':4835928,
               'EastofEngland':6236072,
               'SouthWest':5624696,
               'London':8961989,
               'SouthEast':9180135}

name_of={'London':'London',
        'SouthEast':'South East',
        'EastofEngland':'East of England',
        'SouthWest':'South West',
        'EastMidlands':'East Midlands',
        'WestMidlands':'West Midlands',
        'YorkshireandTheHumber':'Yorkshire and The Humber',
        'NorthEast':'North East',
        'NorthWest':'North West'
        }

#### Probability of positive test #################
# the value in position i is th probability of testing positive i days after infection

#from https://cmmid.github.io/topics/covid19/reports/pcr_profile/pcr_profile_preprint.pdf
#df=pd.read_csv('PCR_curve_summary.csv')
#positive_probability=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]
testable_period=60
delta=2
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

# the time lag is the number of 0s in the first list
#exposure_to_test_probability=[0]+[incubation_probability[i] for i in range(testable_period-1)]
# create the combined probability distribution (distribution of time from exposure to test)
#exp_to_positive_probability=[exposure_to_test_probability[i]*positive_probability[i] for i in range(testable_period)]

#total=sum(exp_to_positive_probability)
#exposure_to_positive_probability=[i/total for i in exp_to_positive_probability]
  

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,4)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]


fig=plt.figure(figsize=(8,4))
plt.subplots_adjust(hspace=0,wspace=0)

#region='London'

m=0
i=0

for region in population_of:
    
    i=i+1
#    n=i % 3
#    m=1+int(n/3)
    
    ax=fig.add_subplot(3,3,i) 
    ax.text(100,4,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,5.2])
    if i%3==1:
        ax.set_yticks([0,1,2,3,4])
    else:
        ax.set_yticks([])
        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=60)
    else:
        ax.set_xticks([])
   
    if i==4: 
        plt.ylabel('Percentage of people')
    #if i==8:
    #    plt.xlabel('Date')

    # get surveillance data for the region    
    df=pd.read_csv(folder+'/'+region+'_surveillance.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    
    date=df['Date'].tolist()
    
    while date:
        d=date.pop(0)
        #print(d)
        day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        # minus 4 to make it a mid-week estimate
        time_in_days.append(day_numerical)
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    
    # convert data frame to lists for plotting
    rate=df['Rate'].tolist()
    upper=df['Upper'].tolist()
    lower=df['Lower'].tolist()
    ONS_day=df['days_since_march1'].tolist()
    date=df['Date'].tolist()
    
    #cases_day=ONS_day#+[ONS_day[-1]+7,ONS_day[-1]+10]

    
    ###### CASES ##########
    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    #print(df.head())
    #print(len(df))
    cases=df['cases'].tolist()
    time_in_days=df['Days_since_March1'].tolist()
    variant_proportion=df['NV_proportion'].tolist()
    PCR_tests=df['PCR_tests'].tolist()
    LFD_tests=df['LFD_tests'].tolist()
    
    ############## COPIED IN ###################
    C=cases+[0 for i in range(100)]
    
    df=pd.read_csv('../raw_data/PCR_curve_summary.csv')
    S=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)] 
    # in for old variants it should just be 0 after the end    
    S=S+[0 for i in range(30,testable_period)]

    # for the variant, choose the peak and the tail value (same as)
    #Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,120,10)]+[np.mean(df['median'].tolist()[120:130]) for i in range(120,300,10)]
    
    #Sv=[i*s_max/peak for i in range(peak)]+[s_max*(testable_period-i)/(testable_period-peak) for i in range(peak,testable_period)]
    #
    
    #stretch=1.6
    
    interval=5
    print('stretch=',10/interval)
    #testable_period=int(30*stretch)
    Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,40,10)]
    Sv=Sv+[np.mean(df['median'].tolist()[i:i+interval]) for i in range(40,300,interval)]
    Sv=Sv+[0 for i in range(len(Sv),testable_period)]
    # s_max=0.39
    # Sv=[s*s_max/0.8 for s in Sv]
    
    #print()
    
    # the time lag is the number of 0s in the first list
    #R=[incubation_probability[i] for i in range(testable_period-1)]
    
    # create the combined probability distribution (distribution of time from exposure to test)
    #exp_to_positive_probability=[exposure_to_test_probability[i]*positive_probability[i] for i in range(testable_period)]
    
    
    theta_times_I=[]
    theta_times_Iv=[]
    for t in range(len(cases)):
        
        Z=0
        
        for j in range(delta,testable_period):
            
            numerator=R[j-delta]*S[j]
            # sum([C[t+k] for k in range(3,10)])*
            denominator=sum([R[i-delta]*S[i] for i in range(testable_period)])
            # sum([C[t+j-i+k] for k in range(3,10)])*
            Z=Z+C[t+j]*(numerator/denominator)
        
        theta_times_I.append(Z)
        
        Zv=0
        
        for j in range(delta,testable_period):
            
            numerator=R[j-delta]*Sv[j]
            # sum([C[t+k] for k in range(3,10)])*
            denominator=sum([R[i-delta]*Sv[i] for i in range(testable_period)])
            # sum([C[t+j-i+k] for k in range(3,10)])*
            Zv=Zv+C[t+j]*(numerator/denominator)
        
        theta_times_Iv.append(Z)
        
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
    
    
    new_ONS=[population_of[region]*r/100 for r in rate]
 
     
  
    ############# END ############
                            
    
    plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
    plt.scatter(ONS_day,rate,s=10,color='k',marker='^')

    #### Calculate th multiplier ###########
    # x is the multiplier
    x_values=[100/i for i in range(1,100)]
    #error=[]
    e_min=10**10
    for x1 in x_values:
        for x2 in x_values:
            #x2=x1 
            #### Calculate an error #######
            #e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
            
            e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
       
            #error.append(e)
            if e<e_min:
                e_min=e
                best_x=(x1,x2)
                
                best_new_cases=new_cases.copy()
                best_nv_cases=nv_cases.copy()
    
    # take the value with the smallest error
    x1,x2=best_x
    #adjusted_cases=[x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i] for i in range(len(new_ONS))]
    adjusted_cases=[x1*best_new_cases[i]*(1-proportion[i])+x2*best_nv_cases[i]*proportion[i] for i in range(len(new_ONS))]
      
     
    plt.scatter(ONS_day,[100*c/population_of[region] for c in adjusted_cases],s=10,facecolors='none', edgecolors='b')
    # add the result to the plot
    ax.text(100,3,'$\\theta_{OLD}='+str(round(1/x1,2))+'$, $\\theta_{VOC}='+str(round(1/x2,2))+'$')
    
plt.savefig('../figures/Regions_different_mutipliers.png',format='png', bbox_inches='tight',dpi=256)


