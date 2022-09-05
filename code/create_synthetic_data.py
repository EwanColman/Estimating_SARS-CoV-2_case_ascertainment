
import matplotlib.pyplot as plt
from numpy.random import poisson
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
            (400,450):1.8,
            (450,500):1,
            (500,550):2,
            (550,600):0.7}

R0=[]
for interval in R0_changes:
    R0=R0+[R0_changes[interval] for i in range(interval[0],interval[1])]



iterations=600
for t in range(iterations):
    
    prev_S=S[-1]
    prev_E=E[-1]
    prev_I=I[-1]
    prev_R=R[-1]
    
    R_rand=normal(R0[t],0.2)
    
    print((1/infectious_period)*R_rand*prev_S*prev_I/N)
    
    # prev_S*prev_I trials with probability (1/infectious_period)*R0*/N
    infections=poisson(max(0,(1/infectious_period)*R_rand*prev_S*prev_I/N))
    
    transitions=poisson(max(0,(1/latent_period)*prev_E))
    
    recoveries=poisson(max(0,(1/infectious_period)*prev_I))
    
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

plt.figure()
plt.plot(new_infections)
#plt.legend()
# project

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
for i in range(len(S)):
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
sample_size=10000
p=sample_size/N
rate=[]
upper=[]
lower=[]
for day in range(0,iterations,7):
    test_positives=sum(poisson(new_infections[day-t]*p)*S_pcr[t] for t in range(min(30,day)))

    sd=(new_infections[day-t]*p)**(1/2)

    rate.append(100*test_positives/sample_size)
    upper.append(100*(test_positives+1.96*sd)/sample_size)
    lower.append(100*(test_positives-1.96*sd)/sample_size)

plt.figure()
plt.plot(range(0,iterations,7),rate)
plt.fill_between(range(0,iterations,7),lower,upper,alpha=0.3)


# for the case numbers project the infections forward to the time of test and 
# account for under ascertainment

test_seeking_rate=[0.5,0.5,0.5,0.5,0.5,0.3,0.3]
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
            while r>slide:
                slide=slide+incubation_probability[incubation_period]
                incubation_period=incubation_period+1
            
            # total time includes a delay
            time_to_test=incubation_period+delta
            # do they get a positive
            if random()<S[time_to_test] and x+time_to_test<len(cases):
                # then add them to cases
                cases[x+time_to_test]=cases[x+time_to_test]+1

plt.figure()
plt.plot(cases)

incidence=[100*i/N for i in new_infections]

pd.DataFrame({'Date':range(iterations),'cases':cases}).to_csv('../synthetic_data/cases.csv')
pd.DataFrame({'Date':range(iterations),'Infections':incidence}).to_csv('../synthetic_data/incidence.csv')
pd.DataFrame({'Date':range(0,iterations,7),'Rate':rate,'Upper':upper,'Lower':lower}).to_csv('../synthetic_data/surveillance.csv')

