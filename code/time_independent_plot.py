
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
import pickle as pk
import scipy.optimize as optimize


### read in the multipliers
rate_df=pd.read_csv('../output/time_independent_rates.csv',index_col='Group')
print(rate_df.head())


delta=1

stretch=1
max_day=825

fs=12

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

#measure days from this day
#time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')
Sep1='23 August 2020'
Sep1=(datetime.strptime(str(Sep1), '%d %B %Y')-time_zero).days

Mar1='1 March 2021'
Mar1=(datetime.strptime(str(Mar1), '%d %B %Y')-time_zero).days


dates=['01/0'+str(i)+'/2020' for i in range(9,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,10)]+['01/'+str(i)+'/2021' for i in range(10,12)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]
dates_words2=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

# get variant prop for ages
df=pd.read_csv('../processed_data/England_daily_data.csv')
SGTF_proportion=df['SGTF_proportion'].tolist()
LFD_proportion=df['LFD_proportion'].tolist()



fig = plt.figure(figsize=(12,16))
gs = fig.add_gridspec(8, 3, height_ratios=[3,2,2,1,2,2,2,2])
plt.subplots_adjust(hspace=0,wspace=0)

testable_period=30#5+int(30*stretch)

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

df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,300,10)]

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

######## AGES ####################################
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]

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

variants=['Wild','Alpha','Delta','BA.1','BA.2']

i=0
for age in name_of:
    print(age)
    #table=table+name_of[age]+' '
#    n=i % 3
#    m=1+int(n/3)
    
    #ax=fig.add_subplot(3,3,i) 
    #ax = fig.add_subplot(gs[2*int(i/3):2+2*int(i/3),(4*(i%3)):(4*(i%3)+4)])
    ax = fig.add_subplot(gs[int(i/3),i%3])
    
    ax.set_xlim([Sep1,max_day])
    
    if i<3:      
        ax.text(Sep1+10,7,'Age '+name_of[age],size=fs)
        ax.set_ylim([0,8])
    else:
        ax.set_ylim([0,5])
        ax.text(Sep1+10,4,'Age '+name_of[age],size=fs)
    if i==0:
        ax.set_yticks([0,1,2,3,4,5,6,7])
        ax.set_yticklabels([0,1,2,3,4,5,6,7],size=fs)
    elif i%3==0:
        ax.set_yticks([0,1,2,3,4])
        ax.set_yticklabels([0,1,2,3,4],size=fs)
    else:
        ax.set_yticks([])   
    if i>3:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words2,ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]
    else:
        ax.set_xticks([])
    if i==3: 
        plt.ylabel('Percentage of people')
        plt.text(Sep1-45,12,'A',size=20)

    # increment for next region
    i=i+1

    #plt.scatter(date,reporting_multiplier,s=10,facecolors='#006666')
    
    #if i==0:
    leg2='Observed test-positive prevalence (ONS CIS)'
    leg1='95% confidence interval'
    leg3='Modelled test-positive prevalence'
    #else:
    #    leg1=None
    #    leg2=None
    #    leg3=None
    df=pd.read_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    df=df[(df['days_since_march1']>Sep1)]#&(df['days_since_march1']<Mar1)] 
    
    # convert data frame to lists for plotting
    rate=df.reset_index()
    ONS_day=df['days_since_march1'].tolist()
    #date=df['Date'].tolist()

 
    #for interval in reporting_multiplier:
    
    interval='Rate'
    new_ONS=[population_of[age]*r/100 for r in rate[interval]]
    ############# END ############
                            
    ax.scatter(ONS_day,rate[interval],s=10,color='k',marker='^',label=leg2)
    ax.fill_between(ONS_day,rate['Lower'],rate['Upper'],color='k',linewidth=0,alpha=0.1,label=leg1) 
   

    #ax.scatter(ONS_day,[100*c/population_of[age] for c in adjusted_cases],s=10,facecolors='none', edgecolors='b',label=leg3)
    # add the result to the plot
    #ax.text(100,3,'$\\theta_{o}='+str(round(100/x1))+'$%, $\\theta_{A}='+str(round(100/x2))+'$%, $\\theta_{\\Delta}='+str(round(100/x3))+'$%',size=fs)
    #table=table+' & $'+str(round(100/x1))+'$ & $'+str(round(100/x2))+'$ \\\ \n'

ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))



######## REGIONS   ##############################



# use different axis labels


# ONS opulation estimates
# pop_df=pd.read_excel('ukmidyearestimates20192020ladcodes',sheet_name='MYE2 - Persons',skiprows=4)
population_of={'SouthWest':5624696,
               'London':8961989,
               'SouthEast':9180135,
               'WestMidlands':5934037,
               'EastMidlands':4835928,
               'EastofEngland':6236072,
               'NorthWest':7341196,
               'YorkshireandTheHumber':5502967,
               'NorthEast':2669941,
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


#region='London'
m=0
i=0
for region in population_of:
    print(region)
    #table=table+name_of[region]+' '
    #ax = fig.add_subplot(gs[11-2*int(i/3):13-2*int(i/3),(4*(i%3)):(4*(i%3)+4)])
    ax = fig.add_subplot(gs[7-int(i/3),i%3])
    #ax=fig.add_subplot(3,3,i) 
    ax.text(Sep1+10,4,name_of[region],size=fs)
    ax.set_xlim([Sep1,max_day])
    ax.set_ylim([0,5])
    if i%3==0:
        ax.set_yticks([0,1,2,3,4])
        ax.set_yticklabels([0,1,2,3,4],size=fs)
    else:
        ax.set_yticks([])
        
    if i<3:
        ax.set_xticks(dates_numerical)
        ax.set_xticklabels(dates_words2,ha='left',size=fs,rotation=90)
    else:
        ax.set_xticks([])
   
    if i==6: 
        plt.ylabel('Percentage of people')
        plt.text(Sep1-45,10,'B',size=20)

    # increment for next region
    i=i+1
    
   #if i==8:
    #    plt.xlabel('Date')

    # get surveillance data for the region    
    df=pd.read_csv('../raw_data/Surveillance/'+region+'_surveillance_1k.csv',sep=',')
    
    # add a new column with header days since Mar1
    time_in_days=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in df['Date'].tolist()]
    # add it to the dataframe
    df['days_since_march1']=time_in_days
    # needs to be in order earliest to latest
    df=df.sort_values('days_since_march1')
    #df.drop(df.tail(1).index,inplace=True)
    df=df[(df['days_since_march1']>Sep1)]#&(df['days_since_march1']<Mar1)] 
    
    
    # convert data frame to lists for plotting
    rate=df['Rate'].tolist()
    upper=df['Upper'].tolist()
    lower=df['Lower'].tolist()
    ONS_day=df['days_since_march1'].tolist()
    date=df['Date'].tolist()
    
    
    
    new_ONS=[population_of[region]*r/100 for r in rate]
    ############# END ############
                            
    
    plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
    plt.scatter(ONS_day,rate,s=10,color='k',marker='^')

    #plt.scatter(ONS_day,[100*c/population_of[region] for c in adjusted_cases],s=10,facecolors='none', edgecolors='b')
    # add the result to the plot
    #ax.text(100,3,'$\\theta_{o}='+str(round(100/x1))+'$%, $\\theta_{A}='+str(round(100/x2))+'$%, $\\theta_{\\Delta}='+str(round(100/x3))+'$%',size=fs)
    # add to the latex output table 
    #table=table+' & $'+str(round(100/x1))+'$ & $'+str(round(100/x2))+'$ \\\ \n'

    # if region=='NorthEast':
    #     print(best_x)
    #     print(adjusted_cases)
        
plt.savefig('../figures/supplementary_figure_time_independant.pdf',format='pdf',dpi=300,bbox_inches='tight')
