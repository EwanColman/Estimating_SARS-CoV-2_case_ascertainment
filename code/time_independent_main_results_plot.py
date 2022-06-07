
import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv('../output/time_independent_rates.csv')

print(df.head())

print(df.iloc[0:7])

print(df.iloc[7:19])

ages=['Age 2 to 10',
      '11 to 15',
      '16 to 24',
      '25 to 34',
      '35 to 49',
      '50 to 69',
      '70+']


regions=['South West',
         'London',
         'South East',
         'West Midlands',
         'East Midlands',
         'East of England',
         'North West',
         'Yorkshire and The Humber',
         'North East',
         'Wales',
         'Scotland',
         'Northern Ireland']

fs=12
offset=0.2

fig=plt.figure(figsize=(12,10))



x_labels=df['Group'].tolist()
y_loc=[]

gs = fig.add_gridspec(2, 1,height_ratios=[7.6,12.6])
plt.subplots_adjust(hspace=0.05,wspace=0)



for i,row in df.iloc[0:7].iterrows():
    
    if i==0:
        lab1,lab2,lab3,lab4,lab5='Wild type','Alpha','Delta','Omicron BA.1','Omicron BA.2'
    else:
        lab1,lab2,lab3,lab4,lab5=None,None,None,None,None
        
    ax = fig.add_subplot(gs[0,0])
    ax.tick_params(length=0)
    plt.plot([row['Wild_Lower'],row['Wild_Upper']],[7-i+offset,7-i+offset],c='k',linewidth=0.5)
    plt.scatter([row['Wild_Rate']],[7-i+offset],color='k',facecolors='k',marker='o', label=lab1)
    y_loc.append(7-i)
    
    plt.plot([row['Alpha_Lower'],row['Alpha_Upper']],[7-i,7-i],c='b',linewidth=0.5,zorder=0)
    plt.scatter([row['Alpha_Rate']],[7-i],color='b',facecolors='w',marker='s',label=lab2,zorder=1)
    
    plt.plot([row['Delta_Lower'],row['Delta_Upper']],[7-i-offset,7-i-offset],c='r',linewidth=0.5,zorder=0)
    plt.scatter([row['Delta_Rate']],[7-i-offset],color='r',facecolors='w',marker='^',label=lab3,zorder=2)

    plt.plot([row['BA.1_Lower'],row['BA.1_Upper']],[7-i-offset,7-i-offset],c='g',linewidth=0.5,zorder=0)
    plt.scatter([row['BA.1_Rate']],[7-i-offset],color='g',facecolors='w',marker='v',label=lab4,zorder=2)
    
    plt.plot([row['BA.2_Lower'],row['BA.2_Upper']],[7-i-offset,7-i-offset],c='m',linewidth=0.5,zorder=0)
    plt.scatter([row['BA.2_Rate']],[7-i-offset],color='m',facecolors='w',marker='<',label=lab5,zorder=2)
    
    
    
    
    if not i==6:
        plt.axhline(7-i-0.5,c='k',linestyle=':',linewidth=0.5)

plt.ylim([0.2,7.8])
plt.yticks(y_loc,ages)
plt.xticks([])
plt.legend(prop={'size':fs})
plt.xlim([0,100])


y_loc=[]
for i,row in df.iloc[7:19].iterrows():
    
        
    ax = fig.add_subplot(gs[1,0]) 
    ax.tick_params(length=0)
    plt.plot([row['Wild_Lower'],row['Wild_Upper']],[19-i+offset,19-i+offset],c='k',linewidth=0.5)
    plt.scatter([row['Wild_Rate']],[19-i+offset],color='k',facecolors='k',marker='o', label=lab1)
    y_loc.append(19-i)
    
    plt.plot([row['Alpha_Lower'],row['Alpha_Upper']],[19-i,19-i],c='b',linewidth=0.5,zorder=0)
    plt.scatter([row['Alpha_Rate']],[19-i],color='b',facecolors='w',marker='s',label=lab2,zorder=1)
    
    plt.plot([row['Delta_Lower'],row['Delta_Upper']],[19-i-offset,19-i-offset],c='r',linewidth=0.5,zorder=0)
    plt.scatter([row['Delta_Rate']],[19-i-offset],color='r',facecolors='w',marker='^',label=lab3,zorder=2)

    plt.plot([row['BA.1_Lower'],row['BA.1_Upper']],[19-i-offset,19-i-offset],c='g',linewidth=0.5,zorder=0)
    plt.scatter([row['BA.1_Rate']],[19-i-offset],color='g',facecolors='w',marker='v',label=lab4,zorder=2)
    
    plt.plot([row['BA.2_Lower'],row['BA.2_Upper']],[19-i-offset,19-i-offset],c='m',linewidth=0.5,zorder=0)
    plt.scatter([row['BA.2_Rate']],[19-i-offset],color='m',facecolors='w',marker='<',label=lab5,zorder=2)
    
    
    if not i==18:
        plt.axhline(19-i-0.5,c='k',linestyle=':',linewidth=0.5)

plt.ylim([0.2,12.9])
plt.yticks(y_loc,regions)


plt.xlim([0,100])

plt.xlabel('Percentage of cases ascertained',size=fs)    

plt.savefig('../figures/figure3.pdf',format='pdf',dpi=300,bbox_inches='tight')
