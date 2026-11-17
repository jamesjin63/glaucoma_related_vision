"""Calculate historical burden, trends, SDI associations and bootstrap frontiers."""
from pathlib import Path
import hashlib, json, itertools, platform
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from settings import INPUT, OUTPUT
O=OUTPUT; T=O/'tables'; D=O/'data'; F=O/'figures'
for p in [T,D,F]: p.mkdir(exist_ok=True,parents=True)
AGES=[16,17,18,19,20,30,31,32,235]
REGIONS=[5,9,21,32,42,56,65,70,73,96,100,104,120,124,134,138,159,167,174,192,199]
SDI=[44635,44634,44639,44636,44637]
COUNTRIES=None
KEY=['measure_id','location_id','sex_id','age_id','cause_id','metric_id','year']
files=sorted((INPUT/'cause').glob('*.csv')); parts=[]; inventory=[]
for p in files:
 x=pd.read_csv(p); x['source_file']=p.name
 inventory.append(dict(file=p.name,sha256=hashlib.sha256(p.read_bytes()).hexdigest(),rows=len(x),locations=x.location_id.nunique(),year_min=int(x.year.min()),year_max=int(x.year.max()),ages=','.join(map(str,sorted(x.age_id.unique())))))
 parts.append(x)
country_export=next(x for x in parts if set(x.age_id)=={228} and x.location_id.nunique()==210 and 'rei_id' not in x)
COUNTRIES=sorted(set(country_export.location_id)-set([1]+SDI))
assert len(COUNTRIES)==204

cause=pd.concat([x for x in parts if 'rei_id' not in x],ignore_index=True)
risk=pd.concat([x for x in parts if 'rei_id' in x],ignore_index=True)
# Check value equality before removing overlapping exports. Location names are not keys.
groups=cause.groupby(KEY,dropna=False)
spread=groups[['val','lower','upper']].max()-groups[['val','lower','upper']].min()
assert (spread.to_numpy()==0).all(), 'Conflicting overlapping cause estimates'
duplicates=len(cause)-len(spread)
cause=cause.drop_duplicates(KEY).copy()
names=cause.dropna(subset=['location_name']).drop_duplicates('location_id').set_index('location_id').location_name
cause['location_name']=cause.location_id.map(names)
assert not cause[['val','lower','upper']].isna().any().any()
assert (cause.lower<=cause.val).all() and (cause.val<=cause.upper).all()
aw=pd.read_csv(D/'standard_weights_source.csv').set_index('age_id').loc[AGES].reset_index()
aw['age_name']=['55–59','60–64','65–69','70–74','75–79','80–84','85–89','90–94','95+']
aw['weight']=aw.Population_weight/aw.Population_weight.sum()
assert np.isclose(aw.weight.sum(),1)
aw.to_csv(T/'TableS0_standard_weights_55plus.csv',index=False)
selected=[1]+SDI+REGIONS+COUNTRIES
age=cause[cause.age_id.isin(AGES)&cause.location_id.isin(selected)].copy()
assert len(age)==len(selected)*2*3*9*32*2
idx=['measure_id','location_id','sex_id','age_id','year']
wide=age.pivot(index=idx,columns='metric_id',values=['val','lower','upper'])
wide.columns=[f'{v}_{m}' for v,m in wide.columns]; wide=wide.reset_index()
wide['population']=wide.val_1/wide.val_3*1e5
assert np.isfinite(wide.population).all() and (wide.population>0).all()
pp=wide.pivot(index=['location_id','sex_id','age_id','year'],columns='measure_id',values='population')
pop_discrep=(pp[2]/pp[5]-1).abs().max()
assert pop_discrep<1e-8
wide=wide.merge(aw[['age_id','age_name','weight']],on='age_id',validate='many_to_one')
wide['location_name']=wide.location_id.map(names)
wide['level']=np.select([wide.location_id.eq(1),wide.location_id.isin(SDI),wide.location_id.isin(REGIONS)],['Global','SDI','Region'],default='Country')
wide['outcome']=wide.measure_id.map({2:'DALYs',5:'Prevalence'})
wide['sex_name']=wide.sex_id.map({1:'Male',2:'Female',3:'Both'})
wide['cause_id']=670
wide.to_csv(D/'verified_age_specific_1990_2021.csv',index=False)
for v in ['val','lower','upper']: wide['weighted_'+v]=wide[v+'_3']*wide.weight
keys=['measure_id','outcome','location_id','location_name','level','sex_id','sex_name','year']
summary=wide.groupby(keys).agg(number=('val_1','sum'),number_lower_sum_marginals=('lower_1','sum'),number_upper_sum_marginals=('upper_1','sum'),population=('population','sum'),ASR=('weighted_val','sum'),ASR_lower_approx=('weighted_lower','sum'),ASR_upper_approx=('weighted_upper','sum'),n_ages=('age_id','nunique')).reset_index()
assert summary.n_ages.eq(9).all()
summary['crude_rate']=summary.number/summary.population*1e5
old=cause[cause.age_id.eq(228)&cause.location_id.isin(selected)]
oldwide=old.pivot(index=['measure_id','location_id','sex_id','year'],columns='metric_id',values=['val','lower','upper'])
oldwide.columns=[f'GBD55_{v}_{m}' for v,m in oldwide.columns]; oldwide=oldwide.reset_index()
summary=summary.merge(oldwide,on=['measure_id','location_id','sex_id','year'],validate='one_to_one')
# Official age-228 count bounds are used rather than summed marginal age bounds.
summary['number_lower']=summary.GBD55_lower_1; summary['number_upper']=summary.GBD55_upper_1
summary['count_vs_age228_relerr']=summary.number/summary.GBD55_val_1-1
summary['crude_vs_age228_relerr']=summary.crude_rate/summary.GBD55_val_3-1
assert summary.count_vs_age228_relerr.abs().max()<1e-8
assert summary.crude_vs_age228_relerr.abs().max()<1e-8
summary['ASR_vs_crude_pct']=100*(summary.ASR/summary.crude_rate-1)
summary['cause_id']=670; summary['age_range']='55+ (nine age groups)'
summary.to_csv(T/'TableS_verified_summary_1990_2021.csv',index=False)
summary[summary.year.isin([1990,2021])].to_csv(T/'Diagnostic_crude_vs_ASR_1990_2021.csv',index=False)

def trend(g,col):
 g=g.sort_values('year'); X=sm.add_constant(g.year-g.year.min()); fit=sm.OLS(np.log(g[col]),X).fit()
 b=fit.params.iloc[1]; ci=fit.conf_int().iloc[1]
 hac=fit.get_robustcov_results(cov_type='HAC',maxlags=3,use_t=True).conf_int()[1]
 return {'EAPC':100*np.expm1(b),'EAPC_CI_lower':100*np.expm1(ci.iloc[0]),'EAPC_CI_upper':100*np.expm1(ci.iloc[1]),'EAPC_HAC_CI_lower':100*np.expm1(hac[0]),'EAPC_HAC_CI_upper':100*np.expm1(hac[1]),'n_years':len(g)}
tk=['measure_id','location_id','sex_id']
tr=pd.DataFrame([{**dict(zip(tk,k)),**trend(g,'ASR')} for k,g in summary.groupby(tk)])
trc=pd.DataFrame([{**dict(zip(tk,k)),**{f'crude_{a}':b for a,b in trend(g,'crude_rate').items()}} for k,g in summary.groupby(tk)])
tr=tr.merge(trc,on=tk)
tr.to_csv(T/'TableS_EAPC_ASR_vs_crude_and_HAC.csv',index=False)
a=summary[summary.year.eq(1990)].drop(columns='year'); b=summary[summary.year.eq(2021)].drop(columns='year')
tab=a.merge(b,on=tk+['outcome','location_name','level','sex_name','cause_id','age_range'],suffixes=('_1990','_2021')).merge(tr,on=tk)
tab['number_change_pct']=100*(tab.number_2021/tab.number_1990-1)
tab['ASR_change_pct']=100*(tab.ASR_2021/tab.ASR_1990-1)
for measure,name in [(5,'Prevalence'),(2,'DALYs')]:
 z=tab[tab.measure_id.eq(measure)]
 z[(z.sex_id.eq(3)&z.level.ne('Country'))|(z.level.eq('Global'))].to_csv(T/f'Table{1 if measure==5 else 2}_{name}_Global_SDI_21regions.csv',index=False)
 z[z.level.eq('Country')].to_csv(T/f'TableS_country_{name}_204.csv',index=False)

# Exact additivity checks; all comparisons use point counts, not quantiles.
checks=[]
for m in [2,5]:
 for s in [1,2,3]:
  g=summary[(summary.measure_id==m)&(summary.sex_id==s)]
  glob=g[g.location_id.eq(1)].set_index('year').number
  for lev in ['Country','Region','SDI']:
   sums=g[g.level.eq(lev)].groupby('year').number.sum()
   checks.append(dict(check='sum_'+lev+'_vs_Global',measure_id=m,sex_id=s,max_relative_error=float((sums/glob-1).abs().max())))
sx=summary.pivot(index=['measure_id','location_id','year'],columns='sex_id',values='number')
checks.append(dict(check='male_plus_female_vs_both',max_relative_error=float(((sx[1]+sx[2])/sx[3]-1).abs().max())))
pd.DataFrame(checks).to_csv(T/'Diagnostic_additivity.csv',index=False)

# Country and regional SDI analyses: one cross-section at a time, never pooled-year P values.
sdi=pd.read_csv(D/'SDI_support_1990_2021.csv')
assert not sdi.duplicated(['location_id','year']).any()
assoc=summary[summary.level.isin(['Country','Region'])&summary.sex_id.eq(3)].merge(sdi[['location_id','year','sdi']],on=['location_id','year'],validate='many_to_one')
assoc.to_csv(D/'SDI_ASR_country_region_1990_2021.csv',index=False)
cor=[]; reg=[]
for (m,lev,year),g in assoc[assoc.year.isin([1990,2021])].groupby(['measure_id','level','year']):
 n=len(g); rho,p=stats.spearmanr(g.sdi,g.ASR)
 se=np.sqrt((1+rho*rho/2)/(n-3)); ci=np.tanh(np.arctanh(rho)+np.array([-1,1])*1.96*se)
 cor.append(dict(measure_id=m,level=lev,year=year,n=n,rho=rho,CI_lower=ci[0],CI_upper=ci[1],p_value=p))
 fit=sm.OLS(np.log(g.ASR),sm.add_constant(-10*g.sdi)).fit(); b=fit.params.iloc[1]; lim=fit.conf_int().iloc[1]
 reg.append(dict(measure_id=m,level=lev,year=year,n=n,percent_difference_per_0_1_lower_SDI=100*np.expm1(b),CI_lower=100*np.expm1(lim.iloc[0]),CI_upper=100*np.expm1(lim.iloc[1]),R_squared=fit.rsquared))
pd.DataFrame(cor).to_csv(T/'Table3_SDI_Spearman_cross_sections.csv',index=False)
pd.DataFrame(reg).to_csv(T/'TableS_SDI_loglinear_effect.csv',index=False)

# Country-cluster bootstrap empirical non-increasing FDH frontier. Seed fixed.
rng=np.random.default_rng(20261113); frontier=[]; curves=[]
for measure,g in assoc[assoc.level.eq('Country')].groupby('measure_id'):
 assert len(g)==204*32
 grid=np.linspace(0,1,401); ids=sorted(g.location_id.unique()); bs=[]
 for it in range(1000):
  chosen=rng.choice(ids,len(ids),replace=True)
  h=g[g.location_id.isin(chosen)].sort_values(['sdi','ASR'],ascending=[True,True])
  xx=h.sdi.to_numpy(); yy=np.minimum.accumulate(h.ASR.to_numpy())
  pos=np.searchsorted(xx,grid,side='right')-1
  bs.append(yy[np.clip(pos,0,len(yy)-1)])
 arr=np.array(bs); ref=np.mean(arr,axis=0)
 # No LOESS needed: the empirical step frontier remains monotone by construction.
 curves.append(pd.DataFrame({'measure_id':measure,'sdi':grid,'reference_rate':ref,'bootstrap_lower':np.quantile(arr,.025,axis=0),'bootstrap_upper':np.quantile(arr,.975,axis=0)}))
 h=g[g.year.eq(2021)].copy(); h['frontier_rate']=np.interp(h.sdi,grid,ref); h['gap']=np.maximum(h.ASR-h.frontier_rate,0)
 frontier.append(h)
front=pd.concat(frontier); curve=pd.concat(curves)
front.to_csv(T/'Table4_frontier_gaps_2021.csv',index=False); curve.to_csv(T/'TableS_frontier_reference_curves.csv',index=False)

# Aggregate risk exports supply original age-228 counts, bounds and crude rates.
assert set(risk.age_id)=={228}
ri=['location_id','location_name','sex_id','sex_name','year']
rr=risk.pivot(index=ri,columns='metric_id',values=['val','lower','upper']);rr.columns=[f'HFPG_{a}_{b}' for a,b in rr.columns];rr=rr.reset_index()
rr=rr.merge(summary[summary.measure_id.eq(2)][['location_id','sex_id','year','number','crude_rate']],on=['location_id','sex_id','year'],validate='one_to_one')
rr['PAF_percent']=100*rr.HFPG_val_1/rr.number
rr['risk_rate_type']='crude rate, age 55+'
rr['cause_id']=670;rr['rei_id']=105;rr['age_id']=228;rr['measure_id']=2
assert np.allclose(rr.HFPG_val_1/(rr.number/rr.crude_rate),rr.HFPG_val_3)
rr.to_csv(T/'Table5_HFPG_crude_rates_counts_PAF.csv',index=False)
rt=pd.DataFrame([{**dict(zip(['location_id','sex_id'],k)),**trend(g,'HFPG_val_3')} for k,g in rr.groupby(['location_id','sex_id'])])
rt['rate_type']='crude';rt.to_csv(T/'TableS_HFPG_crude_EAPC.csv',index=False)
