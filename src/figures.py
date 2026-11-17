"""Generate Figures 1-7, 9 and S1-S3 from shared numerical results."""
from pathlib import Path
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter
from settings import INPUT, OUTPUT
O=OUTPUT; F=O/'figures'; T=O/'tables'; W=O; D=O/'data'; P=T
stage=F
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
def save(fig,name):
 fig.savefig(F/(name+'.pdf'),bbox_inches='tight')
 fig.savefig(F/(name+'.png'),dpi=300,bbox_inches='tight');plt.close(fig)
s=pd.read_csv(T/'TableS_verified_summary_1990_2021.csv')
sd=pd.read_csv(O/'data/SDI_support_1990_2021.csv')
assoc=pd.read_csv(O/'data/SDI_ASR_country_region_1990_2021.csv')
cor=pd.read_csv(T/'Table3_SDI_Spearman_cross_sections.csv')
front=pd.read_csv(T/'Table4_frontier_gaps_2021.csv')
curve=pd.read_csv(T/'TableS_frontier_reference_curves.csv')
aw=pd.read_csv(T/'TableS0_standard_weights_55plus.csv');ages=aw.age_id.tolist();colors=plt.cm.viridis(np.linspace(.05,.95,9))
global_age=[];sdi_age=[]
for z in pd.read_csv(O/'data/verified_age_specific_1990_2021.csv',chunksize=60000):
 global_age.append(z[z.location_id.eq(1)&z.sex_id.eq(3)])
 sdi_age.append(z[z.level.eq('SDI')&z.sex_id.eq(3)])
ga=pd.concat(global_age);sa=pd.concat(sdi_age)

groups=[44635,44634,44639,44636,44637]
fig,ax=plt.subplots(2,5,figsize=(16,7.7))
for row,m in enumerate([5,2]):
 for col,l in enumerate(groups):
  a=ax[row,col];d=sa[(sa.measure_id==m)&sa.location_id.eq(l)]
  for age,c in zip(ages,colors):
   z=d[d.age_id.eq(age)].sort_values('year');a.plot(z.year,z.val_3,color=c,lw=1.1)
  a.set_title(d.location_name.iloc[0]);a.set_xticks([1990,2005,2021]);a.tick_params(labelsize=8);a.set_xlabel('Year')
  if col==0:a.set_ylabel(('A  Prevalence' if m==5 else 'B  DALY')+' age-specific rate\nper 100,000')
fig.legend([Line2D([0],[0],color=c) for c in colors],aw.age_name,title='Age group (years)',loc='lower center',ncol=9,frameon=False,fontsize=9)
fig.suptitle('Age-specific trends within five SDI groups, both sexes, 1990–2021',fontsize=14)
fig.subplots_adjust(top=.89,bottom=.17,wspace=.32,hspace=.4);save(fig,'Fig2_SDI_age_specific_trends')

fig,ax=plt.subplots(1,2,figsize=(14,9.2),layout='constrained')
for a,m,title in zip(ax,[5,2],['A  ASPR','B  ASDR']):
 z=s[s.level.isin(['Region','Global'])&s.year.eq(2021)&s.measure_id.eq(m)]
 order=z[z.sex_id.eq(3)].sort_values('ASR').location_id.tolist()
 for sex,c,off in [(1,'#0072B2',-.19),(2,'#D55E00',.19)]:
  q=z[z.sex_id.eq(sex)].set_index('location_id').loc[order]
  a.barh(np.arange(len(order))+off,q.ASR,height=.35,color=c,alpha=.8,label='Male' if sex==1 else 'Female')
  a.errorbar(q.ASR,np.arange(len(order))+off,xerr=np.vstack([q.ASR-q.ASR_lower_approx,q.ASR_upper_approx-q.ASR]),fmt='none',ecolor='#444444',elinewidth=.7,capsize=1.5)
 a.set_yticks(range(len(order)),z[z.sex_id.eq(3)].set_index('location_id').loc[order].location_name,fontsize=8)
 a.set(title=title,xlabel='Truncated age-standardized rate per 100,000');a.legend(frameon=False,fontsize=9)
fig.suptitle('Regional and global sex differences, ages 55+, 2021',fontsize=14);save(fig,'Fig3_regional_sex_standardized_rates')

chunks=[]
for d in pd.read_csv(O/'data/verified_age_specific_1990_2021.csv',chunksize=60000):
 chunks.append(d[d.location_id.eq(1)&d.year.eq(2021)&d.sex_id.isin([1,2])])
age=pd.concat(chunks);aw=pd.read_csv(O/'tables/TableS0_standard_weights_55plus.csv');ids=aw.age_id.tolist()
age.to_csv(O/'tables/Figure1_age_sex_source.csv',index=False)
fig,axes=plt.subplots(2,1,figsize=(10.5,9.5),layout='constrained')
scale_records=[]
for a,m,label in zip(axes,[5,2],['A  Prevalence','B  DALYs']):
 z=age[age.measure_id.eq(m)];k=z.upper_1.max()/z.upper_3.max()
 right=a.secondary_yaxis('right',functions=(lambda x,k=k:x/k,lambda x,k=k:x*k))
 right.set_ylabel('Age-specific '+('prevalence' if m==5 else 'DALY')+' rate per 100,000')
 right.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x:,.0f}'))
 for sex,c,offset in [(1,'#0072B2',-.19),(2,'#D55E00',.19)]:
  q=z[z.sex_id.eq(sex)].set_index('age_id').loc[ids]
  x=np.arange(9)+offset
  a.fill_between(x,q.lower_3*k,q.upper_3*k,color=c,alpha=.13,zorder=1)
  a.bar(x,q.val_1,width=.35,color=c,alpha=.65,edgecolor=c,linewidth=.6,zorder=2)
  a.errorbar(x,q.val_1,yerr=np.vstack([q.val_1-q.lower_1,q.upper_1-q.val_1]),fmt='none',ecolor='#333333',elinewidth=.8,capsize=2,zorder=3)
  a.plot(x,q.val_3*k,color=c,marker='o',ms=4,lw=1.8,zorder=4)
 a.set_xticks(range(9),aw.age_name);a.set(xlabel='Age group (years)',ylabel='Prevalent cases' if m==5 else 'DALYs',title=label,xlim=(-.65,8.65),ylim=(0,z.upper_1.max()*1.12))
 a.yaxis.set_major_formatter(FuncFormatter(lambda x,pos:f'{x/1000:g}k' if x else '0'))
 a.grid(axis='y',alpha=.15,lw=.5);a.spines['right'].set_visible(True)
 handles=[Patch(facecolor='#0072B2',alpha=.65,label='Male: number'),Patch(facecolor='#D55E00',alpha=.65,label='Female: number'),Line2D([0],[0],color='#0072B2',marker='o',label='Male: rate'),Line2D([0],[0],color='#D55E00',marker='o',label='Female: rate')]
 a.legend(handles=handles,ncol=4,frameon=False,fontsize=9,loc='upper left')
 scale_records.append(dict(measure_id=m,display_multiplier=k,left_axis='Number',right_axis='Age-specific rate per 100000',right_to_left='left_display=right_rate*multiplier',mean_roundtrip_max_error=float(np.max(np.abs(z.val_3*k/k-z.val_3)))))
fig.suptitle('Global glaucoma-related vision loss by age and sex, 2021',fontsize=14)
fig.supxlabel('Bars and whiskers: numbers and 95% UIs. Lines and shading: age-specific rates and 95% UIs.',fontsize=10)
save(fig,'Fig1_age_sex_dual_axis')
assoc=pd.read_csv(D/'SDI_ASR_country_region_1990_2021.csv')
sd=pd.read_csv(D/'SDI_support_1990_2021.csv')
sd=sd[['location_id','year','sdi']].drop_duplicates()
z=s[s.level.isin(['Region','Global'])&s.sex_id.eq(3)].merge(sd,on=['location_id','year'],validate='many_to_one')
out=[]
for m in [5,2]:
 for y in [1990,2021]:
  q=z[(z.measure_id==m)&(z.year==y)&z.level.eq('Region')];assert len(q)==21
  rho,p=stats.spearmanr(q.sdi,q.ASR);out.append(dict(measure_id=m,year=y,n=21,rho=rho,p_value=p))
pd.DataFrame(out).to_csv(P/'Figure4_1990_2021_regional_correlations.csv',index=False)
z.to_csv(P/'Figure4_annual_source.csv',index=False)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,ax=plt.subplots(1,2,figsize=(14,9.4))
locations=z[['location_id','location_name']].drop_duplicates().sort_values('location_name')
palette=plt.cm.turbo(np.linspace(.02,.98,len(locations)))
for a,m,title in zip(ax,[5,2],['A  ASPR','B  ASDR']):
 for (l,name),c in zip(locations.itertuples(index=False,name=None),palette):
  q=z[(z.measure_id==m)&(z.location_id==l)].sort_values('year')
  a.plot(q.sdi,q.ASR,color=c,lw=1.1,ls='--' if l==1 else '-',label=name)
  first,last=q.iloc[0],q.iloc[-1]
  a.scatter(first.sdi,first.ASR,marker='o',facecolors='white',edgecolors=[c],s=36,zorder=4)
  a.scatter(last.sdi,last.ASR,marker='^',color=c,s=30,zorder=4)
 a.set(xlabel='SDI',ylabel='Truncated age-standardized rate per 100,000',title=title)
 for i,row in enumerate([r for r in out if r['measure_id']==m]):
  a.text(.98,.98-i*.055,f"{row['year']}: rho={row['rho']:.3f}, P={row['p_value']:.3g} (n=21)",ha='right',va='top',transform=a.transAxes,fontsize=9)
handles,labels=ax[0].get_legend_handles_labels()
handles += [Line2D([],[],marker='o',mfc='white',mec='black',ls='none'),Line2D([],[],marker='^',color='black',ls='none')]
labels += ['1990 (open circle)','2021 (filled triangle)']
fig.legend(handles,labels,loc='lower center',ncol=4,fontsize=8,frameon=False)
fig.suptitle('Regional SDI trajectories, ages 55+, both sexes, 1990–2021')
fig.subplots_adjust(top=.90,bottom=.29,wspace=.25)
for ext in ['pdf','png']:fig.savefig(W/'figures'/f'Fig4_regional_SDI_trajectories.{ext}',dpi=220,bbox_inches='tight')
plt.close(fig)
fig,ax=plt.subplots(1,2,figsize=(12,5.3),layout='constrained')
for a,m,title in zip(ax,[5,2],['A  ASPR','B  ASDR']):
 z=assoc[assoc.level.eq('Country')&assoc.measure_id.eq(m)&assoc.year.eq(2021)]
 a.scatter(z.sdi,z.ASR,s=20,color='#0072B2',alpha=.75)
 q=cor[cor.level.eq('Country')&cor.measure_id.eq(m)&cor.year.eq(2021)].iloc[0]
 a.set(xlabel='SDI',ylabel=title[3:]+' per 100,000',title=f'{title}: n=204, rho={q.rho:.3f}, P={q.p_value:.3g}')
 for idx in [z.ASR.idxmax(),z.ASR.idxmin()]:
  p=z.loc[idx];a.annotate(p.location_name,(p.sdi,p.ASR),xytext=(5,8),textcoords='offset points',fontsize=8)
 a.margins(.10)
save(fig,'Fig6_country_SDI_associations')

fig,ax=plt.subplots(1,2,figsize=(13,6),layout='constrained')
z=assoc[assoc.level.eq('Country')&assoc.measure_id.eq(2)];q=curve[curve.measure_id.eq(2)]
sc=ax[0].scatter(z.sdi,z.ASR,c=z.year,cmap='Blues',s=7,alpha=.6,rasterized=True)
fig.colorbar(sc,ax=ax[0],label='Year',shrink=.7)
ax[0].set_title('A  Country-year observations, 1990–2021')
z1=front[front.measure_id.eq(2)].merge(s[s.level.eq('Country')&s.measure_id.eq(2)&s.sex_id.eq(3)&s.year.eq(1990)][['location_id','ASR']].rename(columns={'ASR':'ASR_1990'}),on='location_id',validate='one_to_one')
for trend,col,mask in [('Decrease or no change','#0072B2',z1.ASR<=z1.ASR_1990),('Increase','#D55E00',z1.ASR>z1.ASR_1990)]:
 h=z1[mask];ax[1].scatter(h.sdi,h.ASR,color=col,s=22,label=trend,zorder=3)
positive=z1[z1.gap>0]
ax[1].vlines(positive.sdi,positive.frontier_rate,positive.ASR,color='gray',alpha=.2,lw=.6)
for p in z1.nlargest(5,'gap').itertuples():ax[1].annotate(p.location_name,(p.sdi,p.ASR),xytext=(4,5),textcoords='offset points',fontsize=8)
ax[1].set_title('B  2021 rates and gaps; trend from 1990');ax[1].legend(frameon=False,fontsize=8,loc='lower left')
for a in ax:
 a.plot(q.sdi,q.reference_rate,color='black',lw=1.8);a.set(xlim=(0,1),xlabel='SDI',ylabel='ASDR per 100,000',ylim=(0,z.ASR.max()*1.06));a.invert_yaxis()
fig.supxlabel('Reversed vertical axes: zero at the top. Gaps are descriptive rate differences.',fontsize=10)
save(fig,'Fig7_DALY_reference_frontier_and_gaps')

ap=pd.read_csv(T/'TableS_BAPC_age_specific_1990_2050.csv')
fig,ax=plt.subplots(3,6,figsize=(17,10))
for outcome,m in enumerate([5,2]):
 for i,age in enumerate(ages):
  a=ax[i//3,outcome*3+i%3]
  obs=ga[ga.measure_id.eq(m)&ga.age_id.eq(age)].sort_values('year')
  future=ap[ap.measure_id.eq(m)&ap.age_id.eq(age)&ap.sex_id.eq(3)&ap.year.ge(2021)].sort_values('year')
  a.plot(obs.year,obs.val_3,color='#0072B2',lw=1.2);a.plot(future.year,future.rate,color='#D55E00',lw=1.3);a.fill_between(future.year,future.rate_lower,future.rate_upper,color='#D55E00',alpha=.14)
  a.axvline(2021,color='gray',ls=':',lw=.8);a.set_title(aw.iloc[i].age_name,fontsize=10);a.set_xticks([1990,2021,2050]);a.tick_params(labelsize=8)
  if i%3==0:a.set_ylabel(('Prevalence' if m==5 else 'DALY')+' rate per 100,000',fontsize=9)
  if i>=6:a.set_xlabel('Year')
fig.text(.27,.96,'A  Age-specific prevalence rates',ha='center',fontsize=14)
fig.text(.76,.96,'B  Age-specific DALY rates',ha='center',fontsize=14)
fig.subplots_adjust(top=.90,bottom=.10,wspace=.45,hspace=.42)
fig.text(.5,.025,'Global, both sexes. Blue: GBD 2021 observations (1990–2021); orange: conditional projections (2022–2050), with 95% marginal CrIs.',ha='center',fontsize=10)
save(fig,'Fig9_age_specific_BAPC_projections')


summary=s
groups=[1]+[44635,44634,44639,44636,44637]
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for row,(m,title) in enumerate([(5,'Prevalence'),(2,'DALYs')]):
 for col,(v,label) in enumerate([('number','Number'),('ASR','Truncated ASR per 100,000')]):
  a=ax[row,col]; z=s[(summary.measure_id==m)&summary.sex_id.eq(3)]
  for i,y in enumerate([1990,2021]):
   q=z[z.year.eq(y)].set_index('location_id').loc[groups]
   a.bar(np.arange(6)+(i-.5)*.36,q[v],width=.36,label=str(y),color=['#56B4E9','#D55E00'][i])
  a.set_xticks(range(6),['Global','High','High-middle','Middle','Low-middle','Low'],rotation=25);a.set_ylabel(label);a.set_title(title);a.legend(frameon=False)
fig.suptitle('Glaucoma-related vision loss at ages 55+: 1990 and 2021 (GBD 2021)');save(fig,'FigS1_1990_2021_summary')
pr=pd.read_csv(T/'Table7_BAPC_global_annual_1990_2050.csv')
va=pd.read_csv(T/'TableS_temporal_validation_yearly_2012_2021.csv')
fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained')
for i,m in enumerate([5,2]):
 for j,(v,lo,hi,ob,label) in enumerate([('ASR','ASR_lower','ASR_upper','observed_ASR','Truncated ASR per 100,000'),('number','number_lower_approx','number_upper_approx','observed_number','Prevalent cases' if m==5 else 'DALYs')]):
  a=ax[i,j]
  for sex,col in [(3,'#222222'),(1,'#0072B2'),(2,'#D55E00')]:
   z=pr[(pr.measure_id==m)&(pr.sex_id==sex)].sort_values('year');old=z[z.year<=2021];new=z[z.year>=2021]
   a.plot(old.year,old[ob],color=col,label={3:'Both',1:'Male',2:'Female'}[sex]);a.plot(new.year,new[v],color=col,ls='--');a.fill_between(new.year,new[lo],new[hi],color=col,alpha=.12)
  a.axvline(2021,color='gray',ls=':',lw=.8);a.set(xlabel='Year',ylabel=label,title=('Prevalence' if m==5 else 'DALYs'))
ax[0,0].legend(frameon=False);save(fig,'FigS2_global_projected_rates_and_counts')
fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
for a,m in zip(ax,[5,2]):
 z=va[(va.measure_id==m)&(va.sex_id==3)].sort_values('year');a.plot(z.year,z.ASR,color='#0072B2',ls='--',label='BAPC prediction');a.fill_between(z.year,z.ASR_lower,z.ASR_upper,color='#0072B2',alpha=.2);a.scatter(z.year,z.observed_ASR,c='black',label='Held-out GBD estimates');a.plot(z.year,z.naive_ASR,ls=':',color='#D55E00',label='Fixed 2011 age-specific rates')
 a.set(xlabel='Year',ylabel='Truncated ASR per 100,000',title=('Prevalence' if m==5 else 'DALYs')+': fit 1990–2011, predict 2012–2021')
ax[0].legend(frameon=False,fontsize=8);save(fig,'FigS3_temporal_validation')
subprocess.run(['Rscript',str(Path(__file__).with_name('maps.R')),str(F)],check=True)
subprocess.run(['pdfunite',str(F/'map_overview.pdf')]+[str(F/f'map_panel_{x}.pdf') for x in 'ABCDEF']+[str(F/'Fig5_country_maps.pdf')],check=True)
for p in [F/'map_overview.pdf',*[F/f'map_panel_{x}.pdf' for x in 'ABCDEF']]:p.unlink()
