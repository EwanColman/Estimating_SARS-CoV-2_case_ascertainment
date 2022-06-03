import pandas as pd
from datetime import datetime
from datetime import timedelta
from math import isnan
# for correlations use the last T points

folder='../raw_data/Regions'

#measure days from this day
time_zero=datetime.strptime('2020-03-01', '%Y-%m-%d')

### Regions ###
columns={'NorthEast':'B:D',
      'NorthWest':'I:K',
      'YorkshireandTheHumber':'P:R',
      'EastMidlands':'W:Y',
      'WestMidlands':'AD:AF',
      'EastofEngland':'AK:AM',
      'London':'AR:AT',
      'SouthEast':'AY:BA',
      'SouthWest':'BF:BH'}

rows=[6,60]
column_names=['Time']
df=pd.read_excel('../raw_data/raw_surveillance/England.xlsx',sheet_name='1k',usecols='A',names=column_names,skiprows=rows[0],nrows=rows[1]-rows[0]-1).replace('-',0)

time=[]
for t in df['Time'].tolist():
    start=(datetime.strptime(t[:t.find('to ')-1].strip(),'%d %B %Y')-time_zero).days
    end=(datetime.strptime(t[t.find('to ')+3:].strip(),'%d %B %Y')-time_zero).days
    # get middl of interval
    mid_point=start+int((end-start)/2)
    
    new_date=datetime.date((time_zero+timedelta(mid_point))).strftime('%d/%m/%Y')
    time.append(new_date)

column_names=['Rate','Lower','Upper']
for region in columns:
    print(region)

    df=pd.read_excel('../raw_data/raw_surveillance/England.xlsx',sheet_name='1k',usecols=columns[region],names=column_names,skiprows=rows[0],nrows=rows[1]-rows[0]-1).replace('-',0)
    df['Date']=time
    df.to_csv('../raw_data/Surveillance/'+region+'_surveillance_1k.csv',index=False)
    

### ages ###
columns={'02_10':'B:D',
         '11_15':'G:I',
         '16_24':'L:N',
         '25_34':'Q:S',
         '35_49':'V:X',
         '50_69':'AA:AC',
         '70+':'AF:AH'}


rows=[6,60]
column_names=['Time']
df=pd.read_excel('../raw_data/raw_surveillance/England.xlsx',sheet_name='1j',usecols='A',names=column_names,skiprows=rows[0],nrows=rows[1]-rows[0]-1).replace('-',0)
   
time=[]
for t in df['Time'].tolist():
    start=(datetime.strptime(t[:t.find('to ')-1].strip(),'%d %B %Y')-time_zero).days
    end=(datetime.strptime(t[t.find('to ')+3:].strip(),'%d %B %Y')-time_zero).days
    # get middl of interval
    mid_point=start+int((end-start)/2)
    
    new_date=datetime.date((time_zero+timedelta(mid_point))).strftime('%d/%m/%Y')
    time.append(new_date)

column_names=['Rate','Lower','Upper']
for age in columns:
    print(age)
        
    df=pd.read_excel('../raw_data/raw_surveillance/England.xlsx',sheet_name='1j',usecols=columns[age],names=column_names,skiprows=rows[0],nrows=rows[1]-rows[0]-1).replace('-',0)
    df['Date']=time
    df.to_csv('../raw_data/Surveillance/'+age+'_surveillance.csv',index=False)



######## Nations #########


rows={'NorthernIreland':[12,95],
      'Scotland':[11,92],
      'Wales':[7,102]
      }

for nation in rows:
    print(nation)
    df=pd.read_excel('../raw_data/raw_surveillance/'+nation+'.xlsx',sheet_name='1a',usecols='A',names=['Time'],skiprows=rows[nation][0],nrows=rows[nation][1]-rows[nation][0]-1).replace('-',0)
    time=[]
    for t in df['Time'].tolist():
        
        start=(datetime.strptime(t[:t.find('to ')-1].strip(),'%d %B %Y')-time_zero).days
        end=(datetime.strptime(t[t.find('to ')+3:].strip(),'%d %B %Y')-time_zero).days
        # get middl of interval
        mid_point=start+int((end-start)/2)
    
        new_date=datetime.date((time_zero+timedelta(mid_point))).strftime('%d/%m/%Y')
        time.append(new_date)    

    column_names=['Rate','Lower','Upper']          
    df=pd.read_excel('../raw_data/raw_surveillance/'+nation+'.xlsx',sheet_name='1a',usecols='B:D',names=column_names,skiprows=rows[nation][0],nrows=rows[nation][1]-rows[nation][0]-1).replace('-',0)
    df['Date']=time
    df.to_csv('../raw_data/Surveillance/'+nation+'_surveillance_1k.csv',index=False)
    

