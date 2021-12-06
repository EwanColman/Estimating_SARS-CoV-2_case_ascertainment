
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.special import erf
import pickle as pk
import scipy.optimize as optimize

delta=1


fs=12

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

Mar1='1 March 2021'
Mar1=(datetime.strptime(str(Mar1), '%d %B %Y')-time_zero).days
max_day=Mar1+10

dates=['01/0'+str(i)+'/2020' for i in range(7,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,4)]
dates_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]
dates_words2=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b') for d in dates]


fig = plt.figure(figsize=(20,5))
# gs = fig.add_gridspec(18, 12)
gs = gridspec.GridSpec(ncols=8, nrows=1, figure=fig)
plt.subplots_adjust(hspace=0,wspace=0)

########### HEATMAP ######################################
reporting_WT=[i for i in range(1,101)]
reporting_VoC=[i for i in range(1,101)]
s_max_values=[i/100 for i in range(1,100)]
stretch_values=[round(1+i/50,2) for i in range(100)]

error=pk.load(open('../pickles/error_(test_seeking).p','rb'))
#error=pk.load(open('../pickles/error_1.p','rb'))
# find optimum
best_value=10**10
best_x=None
for key, value in error.items(): 
    if value<best_value:
        best_value=value
        best_x=key
print(best_x,np.log(best_value))
percentage_WT=25
interval=5

percentage_WT,percentage_VoC,s_max,stretch=best_x

M=[]

n=1
#for delta in range(1,3):
for p_VoC in reporting_VoC:
    m=[]
    for s in stretch_values:
        x=(percentage_WT,p_VoC,s_max,s)
        m.append(np.log(error[x]))
        #m.append(n)
        #n=n+1
    M.append(m)


Values=np.array(M)



ax = fig.add_subplot(gs[0,5:7])
plt.text(-16,92,'B',size=20)
plt.ylabel('Alpha variant reporting rate, $\\theta_{v}$',size=fs)

plt.xlabel('Alpha variant strech',size=fs)
plt.xticks([20,40,60,80],[1.4,1.8,2.2,2.6],size=fs)
plt.xlim([0,99])
plt.yticks([20,40,60,80],size=fs)


print(max(max(M)),min(min(M)))

plt.contourf(range(100), range(100), Values, cmap='Greys',levels=[13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20])


#plt.imshow(Values,cmap='Greys',origin='lower') # displays in color
cbar=plt.colorbar(fraction=0.046, pad=0.04,ticks=[14,20],orientation='vertical')
cbar.outline.set_visible(False)
cbar.set_label('log(ablosute error)',size=fs,labelpad=10)
cbar.ax.tick_params(labelsize=fs, length=5)
cbar.ax.set_yticklabels([14, 20])

plt.scatter([int(50*(stretch-1))],[percentage_VoC],c='r')

 ######## ENGLAND (ALL) #######################################

 #s_max=0.85
 #interval=5
#x1=100/percentage_WT
#x2=100/percentage_VoC
#gradient=100/(x2*s_max)
 #print('grad=',gradient)
 #plt.plot([0,len(M)],[0,gradient],':g',linewidth=5,label='$\\theta_{v}=0.66\\times s_{max}$')
 #plt.legend(loc=2,prop={'size':fs})

ax = fig.add_subplot(gs[0,0:4])
plt.text(25,2.9,'A',size=20)


#stretch=2
testable_period=5+int(30*stretch)

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



print(date[-2:])
# mak the plot and ad the ONS data





plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
plt.scatter(ONS_day,rate,color='k',marker='^',label='ONS Surveillance: estimated number of test-positive people in England')



###### CASES ##########
df=pd.read_csv('../processed_data/England_daily_data.csv')
cases=df['cases'].tolist()


variant_proportion=df['NV_proportion'].tolist()
# add a bit of the future to prevent indexing errors
C=cases+[0 for i in range(100)]
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
#S=S+[0 for i in range(30,testable_period)]
# in for old variants it should just be 0 after the end


So=S+[0 for i in range(30,testable_period+1)]
Ro=R+[0 for i in range(30,testable_period+1)]
    
pi=0
Sv=So[0:pi]
for tau in range(pi,testable_period):
    f=(tau-pi)/stretch
    #print(tau,pi+int(f)+1,len(So))
    y=So[pi+int(f)]+(f-int(f))*(So[pi+int(f)+1]-So[pi+int(f)])    
    Sv.append(y)

pi=0
Rv=R[0:pi]
for tau in range(pi,testable_period):
    f=(tau-pi)/stretch
    #print(tau,pi+int(f)+1,len(So))
    y=Ro[pi+int(f)]+(f-int(f))*(Ro[pi+int(f)+1]-Ro[pi+int(f)])    
    Rv.append(y)

# Sv=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,40,10)]
# Sv=Sv+[np.mean(df['median'].tolist()[i:i+interval]) for i in range(40,300,interval)]
# Sv=Sv+[0 for i in range(len(Sv),testable_period)]
# Sv=[s*s_max/0.785 for s in Sv]

#Sv=So.copy()
########################################

theta_times_I=[]
theta_times_Iv=[]
for t in range(len(cases)):
    Z=0
    for j in range(delta,testable_period):    
        numerator=Ro[j-delta]*So[j]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([Ro[i-delta]*So[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        Z=Z+C[t+j]*(numerator/denominator)
    
    theta_times_I.append(Z)
    
    Zv=0
    for j in range(delta,testable_period):    
        numerator=Rv[j-delta]*Sv[j]
        # sum([C[t+k] for k in range(3,10)])*
        denominator=sum([Rv[i-delta]*Sv[i] for i in range(testable_period)])
        # sum([C[t+j-i+k] for k in range(3,10)])*
        Zv=Zv+C[t+j]*(numerator/denominator)
    
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
 
population_of_England=56287000
new_ONS=[population_of_England*r/100 for r in rate]



#x_values=[100/i for i in range(1,101)]
#e_min=10**10
#for x1 in x_values:
#    for x2 in x_values:
#        #x2=x1 
#        #### Calculate an error #######
#        #e=sum([abs(x1*new_cases[i]*(1-proportion[i])+x2*new_cases[i]*proportion[i]-new_ONS[i]) for i in range(len(new_ONS))])
#        
#        e=(sum([(x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i]-new_ONS[i])**2 for i in range(len(new_ONS))]))**(1/2)
#   
#        #error.append(e)
#        if e<e_min:
#            e_min=e
#            best_x=(x1,x2)
#            
#            best_new_cases=new_cases.copy()
#            best_nv_cases=nv_cases.copy()
#
## take the value with the smallest error
#x1,x2=best_x

x1,x2=100/percentage_WT,100/percentage_VoC

adjusted_cases=[x1*new_cases[i]*(1-proportion[i])+x2*nv_cases[i]*proportion[i] for i in range(len(new_ONS))]
# plot the manuppulated case numbers
plt.scatter(ONS_day,[100*c/population_of_England for c in adjusted_cases],facecolors='none', edgecolors='b',label='Estimated with $\\theta_{o}='+str(round(100/x1))+'$% and $\\theta_{\\alpha}='+str(round(100/x2))+'$%')

plt.yticks([i/2 for i in range(5)])
plt.ylim([0,3.3])
plt.ylabel('Percentage of people',size=fs)
#plt.title('Estmated coronavirus test-positive people in England')
ax.set_xticks(dates_numerical)
ax.set_xticklabels(dates_words,size=fs,rotation=60)
#plt.xlim([60,260])

#plt.xlabel('Date')
plt.legend(loc=2,prop={'size':fs})

plt.savefig('../figures/figure4.pdf',format='pdf',dpi=300,bbox_inches='tight')
