
import pandas as pd
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# adjust the sensitivity for new variant with this parameter
#alpha=0.9

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

max_day=420
# #################

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

# make the plot of the ONS data
plt.figure(figsize=(8.5,4))
plt.fill_between(ONS_day,lower,upper,color='k',linewidth=0,alpha=0.1) 
plt.scatter(ONS_day,rate,color='k',marker='^',label='ONS Surveillance: estimated number of test-positive people')
plt.yticks([i/2 for i in range(5)])
plt.ylim([0,2.2])
plt.ylabel('Percentage of people')
dates=['01/0'+str(i)+'/2020' for i in range(3,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,5)]
day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]
dates_words=[datetime.strptime(str(d), '%d/%m/%Y').strftime('%b')+' 1' for d in dates]



#plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,dates_words,rotation=0)
plt.xlim([60,max_day])
plt.legend(loc=2)
#plt.xlabel('Date')
plt.savefig('../raw_data_figures/ONS_number_test-positive.pdf',format='pdf', bbox_inches='tight',dpi=256)


###### Read file ##########
df=pd.read_csv('../processed_data/England_daily_data.csv')
cases=df['cases'].tolist()
time_in_days=df['Days_since_March1'].tolist()
variant_proportion=df['NV_proportion'].tolist()
PCR_tests=df['PCR_tests'].tolist()
LFD_tests=df['LFD_tests'].tolist()
All_tests=(df['PCR_tests']+df['LFD_tests']).tolist()
# LFD adjustment (remove a fraction of cases)

# make the plot of the Case data
plt.figure(figsize=(8.5,4))

plt.bar(time_in_days,cases,width=1,color='k',label='Daily Pillar 1 & 2 cases')
plt.yticks([i*20000 for i in range(4)],[i*20 for i in range(4)])
plt.ylim([0,75000])
plt.ylabel('Cases')
plt.text(60,76000,'$\\times 10^{3}$')
# dates=['01/0'+str(i)+'/2020' for i in range(3,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,3)]
# day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]

# #plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,[d[0:5] for d in dates],rotation=0)
plt.xlim([60,max_day])
plt.legend(loc=2)
plt.xlabel('Date')
plt.savefig('../raw_data_figures/Pillar_cases.pdf',format='pdf', bbox_inches='tight',dpi=256)


# make the plot of the PCR_test data
plt.figure(figsize=(8.5,4))

plt.bar(time_in_days,PCR_tests,width=1,color='k',label='Pillar 1 & 2 PCR tests (7 day rolling sum)')
plt.yticks([i*500000 for i in range(6)],[i*0.5 for i in range(6)])
plt.ylim([0,2600000])
plt.ylabel('Tests')
plt.text(60,2650000,'$\\times 10^{6}$')
# dates=['01/0'+str(i)+'/2020' for i in range(3,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,5)]
# day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]

# #plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,[d[0:5] for d in dates],rotation=0)
plt.xlim([60,max_day])
plt.legend(loc=2)
plt.xlabel('Date')
plt.savefig('../raw_data_figures/PCR_tests.pdf',format='pdf', bbox_inches='tight',dpi=256)

# make the plot of the LFD_test data
plt.figure(figsize=(8.5,4))

plt.bar(time_in_days,LFD_tests,width=1,color='k',label='Daily lateral flow device tests')
plt.yticks([i*500000 for i in range(4)],[i*0.5 for i in range(4)])
plt.ylim([0,2000000])
plt.ylabel('Tests')
plt.text(60,2050000,'$\\times 10^{6}$')
# dates=['01/0'+str(i)+'/2020' for i in range(3,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,3)]
# day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]

# #plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,[d[0:5] for d in dates],rotation=0)
plt.xlim([60,max_day])
plt.legend(loc=2)
plt.xlabel('Date')
plt.savefig('../raw_data_figures/LFD_tests.pdf',format='pdf', bbox_inches='tight',dpi=256)


# make the plot of te new variant proportion
plt.figure(figsize=(8.5,4))


x_values=[time_in_days[i] for i in range(0,len(variant_proportion),7)]
proportion=[100*variant_proportion[i] for i in range(0,len(variant_proportion),7)]
          
plt.plot(x_values,proportion,color='k',label='Estimated percentage of cases that are the new variant')
plt.yticks([i*20 for i in range(5)])
plt.ylim([0,100])
plt.ylabel('Percentage')
# dates=['01/0'+str(i)+'/2020' for i in range(3,10)]+['01/'+str(i)+'/2020' for i in range(10,13)]+['01/0'+str(i)+'/2021' for i in range(1,3)]
# day_numerical=[(datetime.strptime(str(d), '%d/%m/%Y')-time_zero).days for d in dates]

# #plt.title('Estmated coronavirus test-positive people in England')
plt.xticks(day_numerical,[d[0:5] for d in dates],rotation=0)
plt.xlim([60,max_day])
plt.legend(loc=2)
plt.xlabel('Date')
plt.savefig('../raw_data_figures/new_variant_proportion.pdf',format='pdf', bbox_inches='tight',dpi=256)

