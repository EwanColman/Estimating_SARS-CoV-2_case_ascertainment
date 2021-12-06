# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 11:13:38 2021

@author: ewanc
"""

import pandas as pd

groups=['00_01','02_10','11_15','16_24','25_34','35_49','50_69']


df=pd.read_excel('../raw_data/Population/ukpopestimatesmid2020on2021geography.xls',sheet_name='MYE2 - Persons',skiprows=7)

row=df[df['Name']=='ENGLAND']
total=0
for group in groups:
    
    p=sum([int(row[str(i)]) for i in range(int(group[0:2]),1+int(group[3:5]))])
   
        
    print(group,p)
    total=total+p

group='70+'
p=sum([int(row[str(i)]) for i in range(70,90)])+int(row['90+'])
print(group,p)
total=total+p

print(total)
#,'70+'

#print(row[])