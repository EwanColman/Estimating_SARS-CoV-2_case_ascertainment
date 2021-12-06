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
        # use theta 9 days in the future as this is when 
        # those exposed on day i are most likely to show up
        # in the survey
        I.append(100*theta_times_I[i]/daily_theta[i+9])
    return I

folder='../raw_data/Surveillance'
cmap = plt.cm.get_cmap('winter')
max_day=530
delta=1

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Jan1='1 January 2021'
Jan1=(datetime.strptime(str(Jan1), '%d %B %Y')-time_zero).days
Nov1='1 November 2020'
Nov1=(datetime.strptime(str(Nov1), '%d %B %Y')-time_zero).days
May9='9 May 2020'
May9=(datetime.strptime(str(May9), '%d %B %Y')-time_zero).days
Sep20='20 September 2020'
Sep20=(datetime.strptime(str(Sep20), '%d %B %Y')-time_zero).days
Sep1='1 September 2020'
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days



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
               'SouthEast':9180135,
               'Wales':3152879,
               'Scotland':5463300,
               'NorthernIreland':1893667
               }

name_of={'London':'London',
        'SouthEast':'South East',
        'EastofEngland':'East of England',
        'SouthWest':'South West',
        'EastMidlands':'East Midlands',
        'WestMidlands':'West Midlands',
        'YorkshireandTheHumber':'Yorkshire and the Humber',
        'NorthEast':'North East',
        'NorthWest':'North West',
        'Wales':'Wales',
        'Scotland':'Scotland',
        'NorthernIreland':'Northern Ireland'
        }

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
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]

df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]



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


# the value in position i is th probability of testing positive i days after infection

# Get list of dates for axes 
dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,6)]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
   

# start making the latex table
table='\\begin{tabular}{l|cccc} \n \\toprule \n & PCR tests & LF Tests & VOC proportion & Cycle threshold \\\ \n \\midrule \n'

# get LFD prop for England
df=pd.read_csv('../processed_data/England_daily_data.csv')
LFD_proportion=df['LFD_proportion'].tolist()


fig=plt.figure(figsize=(12,8))
plt.subplots_adjust(hspace=0,wspace=0)

#region='London'

output_data={}

m=0
i=0
print()
for region in population_of:
    table=table+name_of[region]+' ' 
    
    i=i+1
#    n=i % 3
#    m=1+int(n/3)
    
    ax=fig.add_subplot(4,3,i) 
    ax.text(100,85,name_of[region])
    ax.set_xlim([90,max_day])
    ax.set_ylim([0,100])
    if i%3==1:
        ax.yaxis.set_label_position("left")
        ax.set_yticks([0,20,40,60])
    if i%3==0:
        ax.yaxis.set_label_position("right")
        ax.set_yticks([0,20,40,60])
        ax.set_yticklabels([0,0.2,0.4,0.6])
    else:
        ax.set_yticks([])
        
    if i>6:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words,ha='left',rotation=90)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
   
    if i==4: 
        plt.ylabel('Percentage of infections reported')
    

    # get te data for the region   
    
    #### SURVEILLANCE ########
    df=pd.read_csv(folder+'/'+region+'_surveillance_1k.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[]
    date=df['Date'].tolist()
    while date:
        d=date.pop(0)
        #day_numerical=(datetime.strptime(str(d), '%Y-%m-%d')-time_zero).days
        day_numerical=(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days
        
        if region in ['Wales','Scotland','NorthernIreland']:
            time_in_days.append(day_numerical)
        else:
            time_in_days.append(day_numerical)
    
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1',ascending=True)
    # remove days when the rate is 0
    df=df[df['Rate']>0]
    
    #df.drop(df.tail(1).index,inplace=True)
    
    # convert data frame to lists for plotting
    #rate=df['Rate'].tolist()
    rate=df.reset_index()
    #upper=df['Upper'].tolist()
    #lower=df['Lower'].tolist()
    ONS_day=df['days_since_march1'].tolist()
    date=df['Date'].tolist()
    
    #cut of the beginning 
    #ONS_day=[day for day in ONS_day if day>=Sep1]
    #rate=rate[-len(ONS_day):]
    #print(len(rate),len(ONS_day))
    
    # dates_simce_Nov1=[day for day in ONS_day if day>=Nov1 and day<Jan1]
    # dates_simce_May9=[day for day in ONS_day if day>=May9]
    # dates_simce_Sep20=[day for day in ONS_day if day>=Sep20]
    #dates_simce_Nov1=dates_simce_Nov1[:-2]
    
    ###### CASES ##########
    df=pd.read_csv('../processed_data/'+region+'_daily_data.csv')
    #  drop last one (remove this is not needed)
    
    
    
    cases=df['cases'].tolist()
    
    #time_in_days=df['Days_since_March1'].tolist()
    variant_proportion=df['NV_proportion'].tolist()
    #PCR_tests=df['PCR_tests'].tolist()
    #LFD_tests=df['LFD_tests'].tolist()
    #CT_mean=df['mean_CT'].tolist()
    #OR_only=df['OR_only'].tolist()
    
    # make the fnction and add some 0s to the end
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
    '''
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
   '''
    
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
        
        # repeat for test-seeking
        estimate=sum([theta_times_I[t-j]*S_pcr[j] for j in range(testable_period)])
        #estimate=sum([exposures[t-j]*positive_probability[j] for j in range(20)])     
        #cases_on_day[surveillance_day[t]]=estimate
        new_cases_ts.append(estimate)
        
        new_variant_proportion.append(variant_proportion[t])
  
    #### calculate multiplier for every time point ##########
    reporting_multiplier={'Lower':[],'Upper':[],'Rate':[]}
    #testing_multiplier=[]
     #testing_multiplier=[]
    for j in range(len(new_cases)):
        for interval in reporting_multiplier:
            if rate[interval][j]>0:
           
                reporting_multiplier[interval].append(min(100,100*100*new_cases[j]/(population_of[region]*rate[interval][j])))
     #       testing_multiplier.append(100*100*new_cases_ts[j]/(population_of[region]*rate[j]))
        
            else:
                reporting_multiplier[interval].append(100)
     #       testing_multiplier.append(None)
   # and plot    
    plt.scatter(ONS_day,reporting_multiplier['Rate'],c=new_variant_proportion,s=10,cmap='winter',vmin=0,vmax=1)
    plt.plot(ONS_day,reporting_multiplier['Lower'])
    plt.plot(ONS_day,reporting_multiplier['Upper'])
    
    I=get_incidence(ONS_day,psi_times_I,reporting_multiplier['Rate'])
    
    
    
    plt.fill_between(range(len(cases)),[0 for i in range(len(cases))],[100*100*i/population_of[region] for i in I],color='k',linewidth=0,alpha=0.2)

    
    
    #print(I)
    #print()
    
    
    
    
    
    
    
    
    
    
    
    #plt.scatter(ONS_day,testing_multiplier,c=new_variant_proportion,s=5,marker='x',cmap='winter',vmin=0,vmax=1,facecolors='none')
    #plt.scatter(ONS_day,testing_multiplier,c='w',s=5)
    #plt.bar(ONS_day,testing_multiplier, 1)
    #plt.bar(ONS_day,reporting_multiplier,1)
    # for l in range(2,len(ONS_day)):
    #     plt.plot(ONS_day[l-2:l],reporting_multiplier[l-2:l],color=cmap(new_variant_proportion[l]),linestyle=':',linewidth=1)
    #     plt.plot(ONS_day[l-2:l],testing_multiplier[l-2:l],color=cmap(new_variant_proportion[l]),linewidth=0.5)
 
    # plt.fill_between(ONS_day,reporting_multiplier,testing_multiplier,color='k',linewidth=0,alpha=0.1)
    # plt.fill_between(ONS_day,[0 for j in ONS_day],reporting_multiplier,color='k',linewidth=0,alpha=0.15)
    
    # print(region)
    # print(reporting_multiplier)
    # print()
    # print(new_variant_proportion)
    # print()
    # print(rate)
    # print()
    # print(multiplier)
    # print()
    ### Correlations ######
    #multiplier=[multiplier[i] for i in range(len(multiplier)) if ONS_day[i] in dates_simce_Nov1]
    #### case estimate #######
    # new_cases=[new_cases[i] for i in range(len(new_cases)) if ONS_day[i] in dates_simce_Nov1]
  
    # pearson, p=stats.pearsonr(new_cases,multiplier)
    # print()
    # print('Est, '+name_of[region]+', Pearson r=',pearson,'p=',p)
    # table=table+' & $'+str("%.2f" % pearson)+' $'
    # if p<0.05:
    #     table=table+'*'
    # if p<0.01:
    #     table=table+'*'
    # check since the start of november
    
    #### PCR TESTS #######
    
    # weekly sum with date in midle of week
    # weekly_tests=[sum(PCR_tests[t-4:t+3]) for t in dates_simce_May9]
    # new_multiplier=[multiplier[i] for i in range(len(multiplier)) if ONS_day[i] in dates_simce_May9]
    # # print(len(weely_tests))
    # # print(new_multiplier)
    
    # # check since the start of november
    
    # pearson, p=stats.pearsonr(weekly_tests,new_multiplier)
    # print()
    # print('PCR, '+name_of[region]+', Pearson r=',pearson,'p=',p)
    # table=table+' & $'+str("%.2f" % pearson)+' $'
    # if p<0.05:
    #     table=table+'*'
    # if p<0.01:
    #     table=table+'*'
    # # check since the start of november
    
    # #### LF TESTS #######
    # # weekly sum with date in midle of week
    # weekly_tests=[sum(LFD_tests[t-4:t+3]) for t in dates_simce_Nov1]
    # new_multiplier=[multiplier[i] for i in range(len(multiplier)) if ONS_day[i] in dates_simce_Nov1]
  
    
    # #pearson, p=stats.pearsonr(new_tests,multiplier)
    # #print(name_of[region]+', Pearson r=',pearson,'p=',p)
    # # check since the start of november
    # pearson, p=stats.pearsonr(weekly_tests,new_multiplier)
    # print('LFD, '+name_of[region]+', Pearson r=',pearson,'p=',p)
    # table=table+' & $'+str("%.2f" % pearson)+' $'
    # if p<0.05:
    #     table=table+'*'
    # if p<0.01:
    #     table=table+'*'         

    # ### New Variant proportion ###############
    # v_prop=[variant_proportion[t] for t in dates_simce_Nov1]
    # new_multiplier=[multiplier[i] for i in range(len(multiplier)) if ONS_day[i] in dates_simce_Nov1]
  
    # # check since the start of november
    # pearson, p=stats.pearsonr(v_prop,new_multiplier)
    # print('NVP, '+name_of[region]+', Pearson r=',pearson,'p=',p)
    # table=table+' & $'+str("%.2f" % pearson)+' $' 
    # if p<0.05:
    #     table=table+'*'
    # if p<0.01:
    #     table=table+'*'       
        
    # ### Cycle thresholds ###############
    # ct=[CT_mean[t] for t in dates_simce_Sep20]
    # new_multiplier=[multiplier[i] for i in range(len(multiplier)) if ONS_day[i] in dates_simce_Sep20]
  
    # # check since the start of november
    # pearson, p=stats.pearsonr(ct,new_multiplier)
    # print('CT, '+name_of[region]+', Pearson r=',pearson,'p=',p)
    # table=table+' & $'+str("%.2f" % pearson)+' $' 
    # if p<0.05:
    #     table=table+'*'
    # if p<0.01:
    #     table=table+'*'    
 
   # ### OR only PCR results ###############
   # o=[OR_only[t] for t in dates_simce_Sep20]
   # # check since the start of november
   # pearson, p=stats.pearsonr(o,multiplier)
   # print('OR, '+name_of[region]+', Pearson r=',pearson,'p=',p)
   # table=table+' & $'+str("%.2f" % pearson)+' $' 
   # if p<0.05:
   #     table=table+'*'
   # if p<0.01:
   #     table=table+'*'      
 
    # gend this row of table
   # table=table+'\\\ \n'
    
    output_data['reporting_multiplier_'+region]=reporting_multiplier
    #output_data['testing_multiplier_'+region]=testing_multiplier
    output_data['variant_proportion_'+region]=new_variant_proportion
    output_data['incidence_'+region]=[100*i/population_of[region] for i in I]
    output_data['date_'+region]=ONS_day
    output_data['sgtf_proportion_'+region]=variant_proportion

pk.dump(output_data,open('../pickles/reporting_rates_(region).p','wb'))
 
    

#table=table+'\\bottomrule \n \end{tabular} \n'
#print()
#print(table)
    
# add variant proportion as colour

cax = plt.axes([0.92, 0.3, 0.01, 0.4])
cbar=plt.colorbar(cax=cax,ticks=[0,1])
cbar.set_label('SGTF proportion')
cbar.outline.set_visible(False)
plt.savefig('../figures/time_dependent_multiplier_regions.png',format='png', bbox_inches='tight',dpi=256)



