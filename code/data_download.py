# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 10:49:03 2021

@author: ewanc
"""

import pandas as pd




code={'Wales':'W92000004',
        'Scotland':'S92000003',
        'NorthernIreland':'N92000002'}

for nation in code:
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode='+code[nation]+'&metric=newCasesBySpecimenDate&format=csv'

    df=pd.read_csv(url)
    df.to_csv('../raw_data/Diagnostic/'+nation+'_cases.csv',index=False)
    
    # get admissions
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode='+code[nation]+'&metric=newAdmissions&format=csv'
    df=pd.read_csv(url)
    df.to_csv('../raw_data/Healthcare/'+nation+'_admissions.csv',index=False)
   
    # vaccinations
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode='+code[nation]+'&metric=newPeopleVaccinatedFirstDoseByPublishDate&format=csv'
            
    #url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode='+code[nation]+'&metric=newPeopleVaccinatedFirstDoseByVaccinationDate&format=csv'
    #print(nation)
    df=pd.read_csv(url)
    #print(df.head())
    df.to_csv('../raw_data/Vaccination/'+nation+'_vaccinations.csv',index=False)


code={'SouthWest':'E12000009',
        'London':'E12000007',
        'SouthEast':'E12000008',
        'WestMidlands':'E12000005',
        'EastMidlands':'E12000004',
        'EastofEngland':'E12000006',
        'NorthWest':'E12000002',
        'YorkshireandTheHumber':'E12000005',
        'NorthEast':'E12000001'}

for region in code: 
    
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=region&areaCode='+code[region]+'&metric=newCasesBySpecimenDate&format=csv'
    df=pd.read_csv(url)
    df.to_csv('../raw_data/Diagnostic/'+region+'_cases.csv',index=False)

    # vacination
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=region&areaCode='+code[region]+'&metric=newPeopleVaccinatedFirstDoseByVaccinationDate&format=csv'
    df=pd.read_csv(url)
    df.to_csv('../raw_data/Vaccination/'+region+'_vaccinations.csv',index=False)

    

# cases in England PCR only, LFD only, and all
url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&metric=newCasesLFDConfirmedPCRBySpecimenDate&metric=newCasesLFDOnlyBySpecimenDate&metric=newCasesPCROnlyBySpecimenDate&format=csv'

df=pd.read_csv(url)
df.to_csv('../raw_data/England_cases.csv',index=False)



#### age stuff ####
# categories are ot the same so we make adjustments
url="https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv"

df=pd.read_csv(url)

age_categories=['02_10','11_15','16_24','25_34','35_49','50_69']#,'70+']

groups=['00_04','05_09','10_14','15_19','20_24',
'25_29',
'30_34',
'35_39',
'40_44',
'45_49',
'50_54',
'55_59',
'60_64',
'65_69']

last_groups=['70_74',
'75_79',
'80_84',
'85_89','90+']

dates=list(set(df['date'].tolist()))

output={}
for age in age_categories:
    output[age]=[0 for i in range(len(dates))]
    
output['70+']=[0 for i in range(len(dates))]


# compare every group and every age category
age_groups=[]
for group in groups:
    group_df=df[df['age']==group]
    cases=group_df['cases'].tolist()
    dates=group_df['date'].tolist()

    for age in age_categories:
       
        # get the intervals of the age and the group
        x1,x2=int(group[0:2]),int(group[3:5])+1

        y1,y2=int(age[0:2]),int(age[3:5])+1
        
        # find if they overlap (and by how much)
        overlap=min(x2,y2)-max(x1,y1)
        #print(age,'and',group,', overlap=',overlap)
        # if they do, 
        if overlap>0:
            # add a portion of the cases in the overlapping agegroup to the caase for the age category
            output[age]=[output[age][i]+round((overlap/5)*cases[i]) for i in range(len(dates))]
    
for group in last_groups:
    cases=group_df['cases'].tolist()
    output['70+']=[output['70+'][i]+cases[i] for i in range(len(dates))]

for age in output:
    
    df=pd.DataFrame()
    df['date']=dates
    df['cases']=output[age]
    df.to_csv("../raw_data/Diagnostic/"+age+'_cases.csv',index=False)


### Deaths ###

url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newDeaths28DaysByDeathDateAgeDemographics&format=csv'

df=pd.read_csv(url)
dates=list(set(df['date'].tolist()))

output={}
for age in age_categories:
    output[age]=[0 for i in range(len(dates))]
    
output['70+']=[0 for i in range(len(dates))]

# compare every group and every age category
age_groups=[]
for group in groups:
    group_df=df[df['age']==group]
    deaths=group_df['deaths'].tolist()
    dates=group_df['date'].tolist()

    for age in age_categories:
       
        # get the intervals of the age and the group
        x1,x2=int(group[0:2]),int(group[3:5])+1

        y1,y2=int(age[0:2]),int(age[3:5])+1
        
        # find if they overlap (and by how much)
        overlap=min(x2,y2)-max(x1,y1)
        #print(age,'and',group,', overlap=',overlap)
        # if they do, 
        if overlap>0:
            # add a portion of the cases in the overlapping agegroup to the caase for the age category
            output[age]=[output[age][i]+round((overlap/5)*deaths[i]) for i in range(len(dates))]
    
for group in last_groups:
    deaths=group_df['deaths'].tolist()
    output['70+']=[output['70+'][i]+deaths[i] for i in range(len(dates))]


for age in output:
    
    df=pd.DataFrame()
    df['date']=dates
    df['deaths']=output[age]
    df.to_csv("../raw_data/Death/"+age+'_deaths.csv',index=False)


###### Vaccinations #############

url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=vaccinationsAgeDemographics&format=csv'

df=pd.read_csv(url)

print(df.columns)

df=df.sort_values('date',ascending=True)
#df=df[df['days_since_march1']>=0]

dates=list(set(df['date'].tolist()))

vaccinated={}
for age in age_categories:
    vaccinated[age]=[0 for i in range(len(dates))]
vaccinated['70+']=[0 for i in range(len(dates))]


ages_of={'02_10':[],
         '11_15':[],
         '16_24':['16_17','18_24'],
         '25_34':['25_29','30_34'],
         '35_49':['35_39','40_44','45_49'],
         '50_69':['50_54','55_59','60_64','65_69'],
         '70+':['70_74','75_79','80_84','85_89','90+']}

for group in vaccinated:
    for age in ages_of[group]:
        age_df=df[df['age']==age]
        dates=age_df['date'].tolist()
        new_vaccinations=age_df['newPeopleVaccinatedFirstDoseByVaccinationDate'].tolist()
        
        vaccinated[group]=[vaccinated[group][i]+new_vaccinations[i] for i in range(len(dates))]
        
for age in vaccinated:
    df=pd.DataFrame()
    df['date']=dates
    df['vaccinated']=vaccinated[age]
    df.to_csv("../raw_data/Vaccination/"+age+'_vaccinations.csv',index=False)

