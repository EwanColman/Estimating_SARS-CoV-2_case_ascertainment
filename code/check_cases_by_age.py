# -*- coding: utf-8 -*-
"""
Created on Sat Jul 24 10:49:03 2021

@author: ewanc
"""

import pandas as pd



'''
code={'Wales':'W92000004',
        'Scotland':'S92000003',
        'NorthernIreland':'N92000002'}

for nation in code:
    url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode='+code[nation]+'&metric=newCasesBySpecimenDate&format=csv'
    df=pd.read_csv(url)

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

'''
# cases in England PCR only, LFD only, and all
url='https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDate&metric=newCasesLFDConfirmedPCRBySpecimenDate&metric=newCasesLFDOnlyBySpecimenDate&metric=newCasesPCROnlyBySpecimenDate&format=csv'


df=pd.read_csv(url)
df.to_csv('../raw_data/England_cases.csv',index=False)

print(df.columns)

#### age stuff ####
# categories are ot the same so we make adjustments
url="https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=newCasesBySpecimenDateAgeDemographics&format=csv"

df=pd.read_csv(url)

print(df.columns)

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


total_cases_ages=[0 for i in range(1000)]

# compare every group and every age category
age_groups=[]
for group in groups:
    group_df=df[df['age']==group]
    cases=group_df['cases'].tolist()

    total_cases_ages=[total_cases_ages[i]+cases[i] for i in range(len(cases))]
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
    group_df=df[df['age']==group]
    cases=group_df['cases'].tolist()
    total_cases_ages=[total_cases_ages[i]+cases[i] for i in range(len(cases))]
    output['70+']=[output['70+'][i]+cases[i] for i in range(len(dates))]


total_cases2=[0 for i in range(1000)]
for age in age_categories+['70+']:
    total_cases2=[total_cases2[i]+output[age][i] for i in range(len(cases))]


for c in range(len(total_cases_ages)):
    print(total_cases_ages[c],total_cases2[c])
