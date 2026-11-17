"""Calculate HFPG attributable rates, fractions and Figure 8."""
from pathlib import Path
import hashlib,json,itertools,shutil
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from settings import INPUT, OUTPUT
W=OUTPUT; T=W/'tables'; D=W/'data'
SRC=INPUT/'hfpg'
AGES=[16,17,18,19,20,30,31,32,235]; LOCS=[1,44635,44634,44639,44636,44637]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
 manifest=json.loads((SRC/'metadata/hfpg_download_manifest.json').read_text())
 assert manifest['study']=='GBD 2021' and manifest['version']==8016 and manifest['status']=='complete'
 assert len(manifest['files'])==32 and {f['year'] for f in manifest['files']}==set(range(1990,2022))
 parts=[]; inputs=[]
 for entry in sorted(manifest['files'],key=lambda f:f['year']):
  p=SRC/entry['path'];assert sha(p)==entry['sha256'],p
  a=pd.read_csv(p);assert len(a)==12474 and a.year.eq(entry['year']).all()
  parts.append(a);inputs.append({'path':entry['path'],'rows':len(a),'sha256':sha(p)})
 raw=pd.concat(parts,ignore_index=True).rename(columns={k:k+'_id' for k in ['location','sex','age','cause','rei','measure','metric']})
 s=pd.read_csv(T/'TableS_verified_summary_1990_2021.csv')
 locations=sorted(s.location_id.unique());assert len(locations)==231
 keys=['location_id','sex_id','age_id','metric_id','year']
 expected=set(itertools.product(locations,[1,2,3],AGES,[1,3],range(1990,2022)))
 actual=list(raw[keys].itertuples(index=False,name=None))
 assert len(actual)==len(set(actual))==len(expected) and set(actual)==expected
 assert raw.cause_id.eq(670).all() and raw.rei_id.eq(105).all() and raw.measure_id.eq(2).all()
 assert np.isfinite(raw[['val','lower','upper']]).all().all()
 assert ((raw.lower<=raw.val)&(raw.val<=raw.upper)).all()
 names=s[['location_id','location_name','level']].drop_duplicates()
 raw=raw.merge(names,on='location_id',validate='many_to_one')
 # Match aggregate risk counts and disease population denominators.
 old=pd.read_csv(T/'Table5_HFPG_crude_rates_counts_PAF.csv')
 count=raw[raw.metric_id.eq(1)].groupby(['location_id','sex_id','year']).val.sum().rename('new_number').reset_index()
 qc=count.merge(old,on=['location_id','sex_id','year'],validate='one_to_one')
 count_rel=float((qc.new_number/qc.HFPG_val_1-1).abs().max());assert count_rel<1e-9
 age=pd.read_csv(D/'verified_age_specific_1990_2021.csv');age=age[age.measure_id.eq(2)]
 idx=['location_id','age_id','sex_id','year']
 rate=raw[raw.metric_id.eq(3)].merge(age[idx+['population','val_1']],on=idx,validate='one_to_one').rename(columns={'val_1':'total_glaucoma_DALYs'})
 rate=rate.merge(raw[raw.metric_id.eq(1)][idx+['val']].rename(columns={'val':'attributable_DALYs'}),on=idx,validate='one_to_one')
 rate_err=float((rate.val-rate.attributable_DALYs/rate.population*1e5).abs().max());assert rate_err<1e-7
 rate['PAF_percent']=100*rate.attributable_DALYs/rate.total_glaucoma_DALYs
 assert rate.PAF_percent.between(0,100).all()
 weights=pd.read_csv(T/'TableS0_standard_weights_55plus.csv').set_index('age_id').weight
 assert np.isclose(weights.sum(),1)
 for v in ['val','lower','upper']:rate['weighted_'+v]=rate[v]*rate.age_id.map(weights)
 gkeys=['location_id','location_name','level','sex_id','year']
 annual=rate.groupby(gkeys).agg(ASDR=('weighted_val','sum'),ASDR_lower_approx=('weighted_lower','sum'),ASDR_upper_approx=('weighted_upper','sum'),attributable_DALYs=('attributable_DALYs','sum'),total_glaucoma_DALYs=('total_glaucoma_DALYs','sum'),population=('population','sum')).reset_index()
 annual['crude_rate']=annual.attributable_DALYs/annual.population*1e5
 annual['PAF_percent']=100*annual.attributable_DALYs/annual.total_glaucoma_DALYs
 assert ((annual.ASDR_lower_approx<=annual.ASDR)&(annual.ASDR<=annual.ASDR_upper_approx)).all()
 trends=[]
 for (loc,sex),z in annual.groupby(['location_id','sex_id']):
  z=z.sort_values('year');x=sm.add_constant(z.year-1990);fit=sm.OLS(np.log(z.ASDR),x).fit()
  ci=fit.conf_int().iloc[1].to_numpy();hac=fit.get_robustcov_results(cov_type='HAC',maxlags=3,use_t=True).conf_int()[1]
  independent=stats.linregress(z.year,np.log(z.ASDR))
  assert np.isclose(fit.params.iloc[1],independent.slope,atol=1e-12)
  trends.append({'location_id':loc,'sex_id':sex,'EAPC':100*np.expm1(fit.params.iloc[1]),'EAPC_CI_lower':100*np.expm1(ci[0]),'EAPC_CI_upper':100*np.expm1(ci[1]),'EAPC_HAC_CI_lower':100*np.expm1(hac[0]),'EAPC_HAC_CI_upper':100*np.expm1(hac[1]),'n_years':len(z)})
 tr=pd.DataFrame(trends)
 for l,u in [('EAPC_CI_lower','EAPC_CI_upper'),('EAPC_HAC_CI_lower','EAPC_HAC_CI_upper')]:assert ((tr[l]<=tr.EAPC)&(tr.EAPC<=tr[u])).all()
 summary=annual[annual.year.eq(1990)].drop(columns='year').merge(annual[annual.year.eq(2021)].drop(columns='year'),on=['location_id','location_name','level','sex_id'],suffixes=('_1990','_2021'),validate='one_to_one').merge(tr,on=['location_id','sex_id'],validate='one_to_one')
 summary['ASDR_change_pct']=100*(summary.ASDR_2021/summary.ASDR_1990-1)
 for name,frame in [('TableS_HFPG_annual_ASDR_crude_PAF.csv',annual),('TableS_HFPG_ASDR_EAPC_1990_2021.csv',summary),('TableS_HFPG_age_specific_rates_PAF.csv',rate)]:frame.to_csv(T/name,index=False)
 source=rate[rate.location_id.isin(LOCS)&rate.sex_id.eq(3)&rate.year.eq(2021)].copy()
 source['age_name']=source.age_id.map(dict(zip(AGES,['55–59','60–64','65–69','70–74','75–79','80–84','85–89','90–94','95+'])))
 source.to_csv(T/'TableS_Fig8B_54_age_specific_rates_UI_PAF.csv',index=False)
 six=summary[summary.location_id.isin(LOCS)&summary.sex_id.eq(3)].set_index('location_id').loc[LOCS].reset_index()
 six.to_csv(T/'TableS_HFPG_Global_5SDI_summary.csv',index=False)
 old=old.drop(columns=['ASDR','ASDR_lower_approx','ASDR_upper_approx'],errors='ignore').merge(annual[['location_id','sex_id','year','ASDR','ASDR_lower_approx','ASDR_upper_approx']],on=['location_id','sex_id','year'],validate='one_to_one')
 old['risk_rate_type']='Original age228 crude rates plus independently derived truncated ASDR; ASDR bounds are weighted marginals'
 old.to_csv(T/'Table5_HFPG_crude_rates_counts_PAF.csv',index=False)
 # Single complete figure, matching the original A-trend/B-age purposes.
 plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
 fig,axes=plt.subplots(1,2,figsize=(14,5.8));colors=['#222222','#0072B2','#56B4E9','#009E73','#E69F00','#D55E00']
 for loc,color in zip(LOCS,colors):
  z=annual[annual.location_id.eq(loc)&annual.sex_id.eq(3)].sort_values('year');label=z.location_name.iloc[0]
  axes[0].plot(z.year,z.ASDR,label=label,color=color,lw=1.8)
  q=source[source.location_id.eq(loc)].set_index('age_id').loc[AGES]
  axes[1].plot(range(9),q.val,label=label,color=color,lw=1.6,marker='o',ms=3)
 axes[0].set(title='A  Truncated age-standardized attributable DALY rates',xlabel='Year',ylabel='HFPG-attributable ASDR per 100,000',xticks=[1990,2000,2010,2021])
 axes[1].set(title='B  Age-specific attributable DALY rates, 2021',xlabel='Age group (years)',ylabel='HFPG-attributable DALY rate per 100,000',xticks=range(9),xticklabels=['55–59','60–64','65–69','70–74','75–79','80–84','85–89','90–94','95+'])
 axes[1].tick_params(axis='x',rotation=35)
 for ax in axes:ax.set_ylim(bottom=0)
 fig.legend(*axes[0].get_legend_handles_labels(),loc='lower center',ncol=6,frameon=False,fontsize=9,bbox_to_anchor=(.5,.055))
 fig.suptitle('HFPG-attributable glaucoma burden, ages 55+, both sexes, GBD 2021',fontsize=14)
 fig.text(.5,.015,'Point estimates. Original age-specific UIs and approximate standardized bounds, including negative lower limits, are retained in electronic tables.',ha='center',fontsize=9)
 fig.subplots_adjust(top=.84,bottom=.25,wspace=.25)
 for ext in ['pdf','png']:fig.savefig(W/'figures'/f'Fig8_HFPG_standardized_trends_and_age_specific.{ext}',dpi=240,bbox_inches='tight')
 plt.close(fig)
 print('HFPG estimates and Figure 8 written to', W)

if __name__=='__main__':main()
