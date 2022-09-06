import matplotlib.pyplot as plt
import pickle as pk


output_data=pk.load(open('../pickles/reporting_rates_(synthetic)_200.p','rb'))

days=output_data['date']

fs=15

fig=plt.figure()
#ax=fig.add_subplot(3,3,i) 
ax = fig.add_subplot(111)

ax.set_xlim([-10,610])
#ax.axvline((datetime.strptime('1 March 2021', '%d %B %Y')-time_zero).days,linewidth=1)
ax.set_ylim([0,100])

ax.set_yticks([0,20,40,60,80])
ax.set_yticklabels([0,20,40,60,80],size=fs)

ax.set_xticks([100,200,300,400,500])
ax.set_xticklabels([100,200,300,400,500],ha='left',rotation=90,size=fs)#[d[0:5] for d in dates]

plt.ylabel('Percentage of infections',size=fs)
      

# medians and confidence intervals
reporting_multiplier_list=output_data['reporting_multiplier']
reporting_multiplier={'Lower':[],'Rate':[],'Upper':[]}
for t in range(len(reporting_multiplier_list[0])):
    rate_list=sorted([reporting_multiplier_list[j][t] for j in range(200)])
    
    reporting_multiplier['Lower'].append(rate_list[4])
    reporting_multiplier['Rate'].append(rate_list[99])
    reporting_multiplier['Upper'].append(rate_list[194])

incidence_list=output_data['incidence']
incidence={'Lower':[],'Rate':[],'Upper':[]}
for t in range(len(incidence_list[0])):  
    inci_list=sorted([incidence_list[j][t] for j in range(200)])
    incidence['Lower'].append(inci_list[4])
    incidence['Rate'].append(inci_list[99])
    incidence['Upper'].append(inci_list[194])
##########
    
#if i==1:
leg1='Ascertainment rate'
leg2='Ascertainment rate 95% confidence interval'


ax.plot(days,reporting_multiplier['Rate'],linestyle=':',marker='o',markersize=2,color='k',linewidth=1,zorder=1,label=leg1)
ax.fill_between(days,reporting_multiplier['Lower'],reporting_multiplier['Upper'],color='m',alpha=0.1,linewidth=1,label=leg2)


ax2 = ax.twinx()
ax2.set_xlim([-10,610])
ax2.set_ylim([0,2])

ax2.set_yticks([0,0.5,1,1.5])
ax2.set_yticklabels([0,0.5,1,1.5],size=fs)

ax2.set_ylabel('Percentage of population',size=fs)

#I=data['incidence_'+age]
#plt.fill_between(range(Sep1,len(I)),[0 for i in range(Sep1,len(I))],I[Sep1:],color='k',linewidth=0,alpha=0.1)

I=incidence['Rate']
ax2.plot(range(len(I)),I,color='g',linewidth=0.75,label='Daily new infections')


# fake plot
#ax.fill_between([0,0],[0,0],[0,0],color='m',linewidth=0,edgecolor='k',label=leg2)
ax.legend(loc=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.5))
ax2.legend(loc=2,ncol=2,prop={'size':fs},frameon=False,bbox_to_anchor=(1.2, 0.2))


