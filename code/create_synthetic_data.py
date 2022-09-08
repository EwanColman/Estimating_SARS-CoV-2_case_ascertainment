
import matplotlib.pyplot as plt
from numpy.random import poisson
from numpy.random import binomial
from numpy.random import normal
from numpy.random import random
import numpy as np
from scipy.special import erf

import pandas as pd


R0=2
infectious_period=5
latent_period=2

S=[10000000] # initial conditions
E=[100]
I=[100] # seed infections
R=[0] 
new_infections=[]

N=S[-1]+I[-1]+R[-1]

R0_changes={(0,250):1.2,
            (250,400):1,
            (400,450):2,
            (450,500):1,
            (500,550):2.5,
            (550,600):2.5}

R0=[]
for interval in R0_changes:
    R0=R0+[R0_changes[interval] for i in range(interval[0],interval[1])]



iterations=600
for t in range(iterations):
    
    prev_S=S[-1]
    prev_E=E[-1]
    prev_I=I[-1]
    prev_R=R[-1]
    
    R_rand=normal(R0[t],0)
    
    print((1/infectious_period)*R_rand*prev_S*prev_I/N)
    
    # prev_S*prev_I trials with probability (1/infectious_period)*R0*/N
    #infections=poisson(max(0,(1/infectious_period)*R_rand*prev_S*prev_I/N))
    infections=round((1/infectious_period)*R_rand*prev_S*prev_I/N)
    
    
    #transitions=poisson(max(0,(1/latent_period)*prev_E))
    transitions=round((1/latent_period)*prev_E)
    
    #recoveries=poisson(max(0,(1/infectious_period)*prev_I))
    recoveries=round((1/infectious_period)*prev_I)
    
    next_S=prev_S-infections
    next_E=prev_E+infections-transitions
    next_I=prev_I+transitions-recoveries
    next_R=prev_R+recoveries

    S.append(next_S)
    E.append(next_E)
    I.append(next_I)
    R.append(next_R)

    new_infections.append(infections)

# plt.plot(S,label='S')
# plt.plot(I,label='E')
# plt.plot(I,label='I')
# plt.plot(R,label='R')


testable_period=30
#### Probability of positive test #################
# the value in position i is th probability of testing positive i days after infection
#from https://www.medrxiv.org/content/10.1101/2020.11.24.20229948v1.full.pdf
df=pd.read_csv('../raw_data/PCR_curve_summary.csv')
S_pcr=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]

#S=S+[0 for i in range(30,testable_period)]
df=pd.read_csv('../raw_data/LFD_curve_summary.csv')
S_lfd=[np.mean(df['median'].tolist()[i:i+10]) for i in range(0,10*testable_period,10)]

# from https://www.acpjournals.org/doi/full/10.7326/M20-0504
mean=5.5
dispersion=1.52

sigma=np.log(dispersion)
mu=np.log(mean)    

incubation_probability=[]
for i in range(testable_period):
    x=i+1
    u=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    x=i
    l=(1/2)*(1+erf((np.log(x)-mu)/(sigma*(2**(1/2)))))
    p=u-l
    incubation_probability.append(p) 

# for the synthetic ONS, take a random sample of the population
# take 10000 people - same as take anyone with probability

# find the number of days since their infection
# there are new_infections[t] people who were infected t days ago 
# poisson with mean new_infections[t]*p 
#sample_size=10000
#p=sample_size/N
p=0.0027
sample_size=p*N
rate=[]
upper=[]
lower=[]
for day in range(0,iterations,7):
    test_positives=sum(binomial(new_infections[day-t],p*S_pcr[t]) for t in range(min(30,day)))

    sd=test_positives**(1/2)

    rate.append(100*test_positives/sample_size)
    upper.append(100*(test_positives+1.96*sd)/sample_size)
    lower.append(100*(test_positives-1.96*sd)/sample_size)


# for the case numbers project the infections forward to the time of test and 
# account for under ascertainment

test_seeking_rate=[0.4,0.4,0.4,0.4,0.4,0.4,0.4]
delta=1

cases=[0 for i in range(iterations)]
for x in range(iterations):
    
    for infection in range(new_infections[x]):
        # do they get a test?
        if random()<test_seeking_rate[x%7]:
            # for each one find a time until test
            r=random()
            # to draw a number from the incubation period distribution function
            # slide up through incubation_probability distribution 
            slide=0
            incubation_period=0
            while r>slide and incubation_period+delta+1<testable_period:
                
                slide=slide+incubation_probability[incubation_period]
                incubation_period=incubation_period+1
            
            # total time includes a delay
            time_to_test=incubation_period+delta

            # do they get a positive
            if random()<S_pcr[time_to_test] and x+time_to_test<len(cases):
                # then add them to cases
                cases[x+time_to_test]=cases[x+time_to_test]+1



incidence=[100*i/N for i in new_infections]

pd.DataFrame({'Date':range(iterations),'cases':cases}).to_csv('../synthetic_data/cases.csv')
pd.DataFrame({'Date':range(iterations),'Infections':incidence}).to_csv('../synthetic_data/incidence.csv')
pd.DataFrame({'Date':range(0,iterations,7),'Rate':rate,'Upper':upper,'Lower':lower}).to_csv('../synthetic_data/surveillance.csv')

### Plotting

fig=plt.figure(figsize=(11,2))
plt.subplots_adjust(wspace=0.4)

ax = fig.add_subplot(1,3,1)
ax.plot([i/1000 for i in new_infections])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Infections (thousands)')
ax.set_title('Synthetic underlying incidence')
ax.text(-0.2,1.1,'A',transform=ax.transAxes,size=15)


ax = fig.add_subplot(1,3,2)
ax.plot([c/1000 for c in cases])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Cases (thousands)')
ax.set_title('Synthetic diagnostic data')
ax.text(-0.2,1.1,'B',transform=ax.transAxes,size=15)


ax = fig.add_subplot(1,3,3)
ax.plot(range(0,iterations,7),rate,linewidth=0.5)
ax.fill_between(range(0,iterations,7),lower,upper,color='k',alpha=0.3)
ax.set_xlabel('Time (days)')
ax.set_ylabel('% testing positive')
ax.set_title('Synthetic surveillance data')
ax.text(-0.2,1.1,'C',transform=ax.transAxes,size=15)


plt.savefig('../figures/synthetic_data.pdf',format='pdf',dpi=256,bbox_inches='tight')




