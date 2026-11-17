"""Generate country and HFPG supplementary tables from original exports."""
from pathlib import Path
import json, hashlib
import numpy as np
import pandas as pd
from scipy import stats
from settings import INPUT, OUTPUT
OUT=OUTPUT
AGES=[16,17,18,19,20,30,31,32,235]
LOCS=[1,44635,44634,44639,44636,44637]
SEX={3:'Both sexes',1:'Male',2:'Female'}
def read(p):return pd.read_csv(p)
def check(name,value,limit=1e-8):
    if float(value)>limit:raise ValueError(f'{name}: {value} exceeds {limit}')
def compare(name,a,b):
    check(name,np.max(np.abs(np.asarray(a)-np.asarray(b))/np.maximum(np.abs(np.asarray(b)),1)))
def trends(df, group, rate):
    records=[]
    for key,z in df.groupby(group):
        z=z.sort_values('year'); assert z.year.tolist()==list(range(1990,2022))
        assert (z[rate]>0).all()
        fit=stats.linregress(z.year,np.log(z[rate]))
        lo,hi=fit.slope+np.array([-1,1])*stats.t.ppf(.975,30)*fit.stderr
        records.append(dict(zip(group,key if isinstance(key,tuple) else [key]),EAPC=100*np.expm1(fit.slope),EAPC_CI_lower=100*np.expm1(lo),EAPC_CI_upper=100*np.expm1(hi),n_years=32))
    return pd.DataFrame(records)

parts=[read(p) for p in sorted((INPUT/'cause').glob('*.csv'))]
raw=pd.concat([x for x in parts if 'rei_id' not in x],ignore_index=True)
key=['location_id','cause_id','measure_id','metric_id','age_id','sex_id','year']
groups=raw.groupby(key)[['val','lower','upper']]
spread=groups.max()-groups.min()
check('conflicting raw disease duplicates',spread.to_numpy().max(),0)
raw=raw.drop_duplicates(key)
countries=read(INPUT/'support/countries.csv')
assert countries.location_id.nunique()==len(countries)==204
raw=raw[raw.location_id.isin(countries.location_id.tolist()+LOCS)&raw.age_id.isin(AGES+[228])].copy()
assert raw.cause_id.eq(670).all()
assert set(raw.measure_id)=={2,5} and set(raw.metric_id)=={1,3}
assert len(raw)==210*2*2*10*3*32
assert ((raw.lower<=raw.val)&(raw.val<=raw.upper)).all()
w=read(OUT/'tables/TableS0_standard_weights_55plus.csv').set_index('age_id').loc[AGES]
check('normalized standard weight sum',abs(w.weight.sum()-1),1e-14)
compare('weight normalization',w.weight,w.Population_weight/w.Population_weight.sum())
w.reset_index().to_csv(OUT/'tables/Standard_weights_55plus.csv',index=False)
idx=['location_id','measure_id','sex_id','age_id','year']
age=raw[raw.age_id.isin(AGES)].pivot(index=idx,columns='metric_id',values=['val','lower','upper'])
age.columns=[f'{a}_{b}' for a,b in age.columns];age=age.reset_index()
age['population']=age.val_1/age.val_3*1e5
assert np.isfinite(age.population).all() and (age.population>0).all()
pop=age.pivot(index=['location_id','sex_id','age_id','year'],columns='measure_id',values='population')
compare('population inferred independently from prevalence and DALYs',pop[2],pop[5])
for v in ['val','lower','upper']:age[f'weighted_{v}']=age[f'{v}_3']*age.age_id.map(w.weight)
g=['location_id','measure_id','sex_id','year']
annual=age.groupby(g).agg(number=('val_1','sum'),population=('population','sum'),ASR=('weighted_val','sum'),ASR_lower_weighted_marginal=('weighted_lower','sum'),ASR_upper_weighted_marginal=('weighted_upper','sum'),n_ages=('age_id','nunique')).reset_index()
assert annual.n_ages.eq(9).all()
agg=raw[raw.age_id.eq(228)].pivot(index=g,columns='metric_id',values=['val','lower','upper'])
agg.columns=[f'{a}_{b}' for a,b in agg.columns];agg=agg.reset_index()
annual=annual.merge(agg,on=g,validate='one_to_one')
compare('nine age counts sum versus original age228 counts',annual.number,annual.val_1)
annual['crude_rate']=annual.number/annual.population*1e5
compare('crude rates versus original age228 rates',annual.crude_rate,annual.val_3)
tr=trends(annual,['location_id','measure_id','sex_id'],'ASR')
country_tables={}
for m,label in [(5,'Prevalence'),(2,'DALYs')]:
    z=annual[annual.measure_id.eq(m)&annual.location_id.isin(countries.location_id)].copy()
    columns=['val_1','lower_1','upper_1','population','crude_rate','ASR','ASR_lower_weighted_marginal','ASR_upper_weighted_marginal']
    endpoints=[]
    for year in [1990,2021]:
        q=z[z.year.eq(year)][['location_id','sex_id']+columns].rename(columns={'val_1':'number','lower_1':'number_lower_95UI','upper_1':'number_upper_95UI'})
        q=q.rename(columns={c:f'{c}_{year}' for c in q if c not in ['location_id','sex_id']});endpoints.append(q)
    result=endpoints[0].merge(endpoints[1],on=['location_id','sex_id'],validate='one_to_one').merge(tr[tr.measure_id.eq(m)],on=['location_id','sex_id'],validate='one_to_one').merge(countries,on='location_id',validate='many_to_one')
    result['sex']=result.sex_id.map(SEX);result['cause_id']=670
    result['sex_order']=result.sex_id.map({3:0,1:1,2:2})
    result=result.sort_values(['location_name','sex_order']).drop(columns='sex_order')
    result=result[['location_name','location_id','sex','sex_id','cause_id','measure_id']+[c for c in result if c not in ['location_name','location_id','sex','sex_id','cause_id','measure_id']]]
    assert len(result)==612 and not result.duplicated(['location_id','sex_id']).any()
    country_tables[f'Country_{label}']=result
    result.to_csv(OUT/f'tables/eTable{1 if m==5 else 2}_Country_{label}.csv',index=False)
hf=[]
meta=json.loads((INPUT/'hfpg/metadata/hfpg_download_manifest.json').read_text())
assert meta['study']=='GBD 2021' and meta['version']==8016 and len(meta['files'])==32
for entry in meta['files']:
    p=INPUT/'hfpg'/entry['path']
    assert hashlib.sha256(p.read_bytes()).hexdigest()==entry['sha256']
    x=read(p).rename(columns={c:c+'_id' for c in ['location','sex','age','cause','rei','measure','metric']})
    assert len(x)==12474 and x.year.eq(entry['year']).all()
    hf.append(x)
hf=pd.concat(hf,ignore_index=True)
assert len(hf)==399168 and not hf.duplicated(['location_id','sex_id','age_id','metric_id','year']).any()
assert hf.cause_id.eq(670).all() and hf.rei_id.eq(105).all() and hf.measure_id.eq(2).all()
assert ((hf.lower<=hf.val)&(hf.val<=hf.upper)).all()
hf=hf[hf.location_id.isin(LOCS)&hf.sex_id.eq(3)]
h=hf.pivot(index=['location_id','sex_id','age_id','year'],columns='metric_id',values=['val','lower','upper'])
h.columns=[f'{a}_{b}' for a,b in h.columns];h=h.reset_index()
den=age[age.measure_id.eq(2)][['location_id','sex_id','age_id','year','population','val_1']].rename(columns={'val_1':'total_glaucoma_DALYs'})
h=h.merge(den,on=['location_id','sex_id','age_id','year'],validate='one_to_one')
assert len(h)==6*9*32
compare('HFPG original rates versus counts/population',h.val_3,h.val_1/h.population*1e5)
for v in ['val','lower','upper']:h['weighted_'+v]=h[f'{v}_3']*h.age_id.map(w.weight)
names=raw[['location_id','location_name']].dropna().drop_duplicates('location_id').set_index('location_id').location_name
h['location_name']=h.location_id.map(names);h['age_group']=h.age_id.map(w.age_name)
a=h.groupby(['location_id','location_name','sex_id','year']).agg(attributable_DALYs=('val_1','sum'),population=('population','sum'),total_glaucoma_DALYs=('total_glaucoma_DALYs','sum'),ASDR=('weighted_val','sum'),ASDR_lower_weighted_marginal=('weighted_lower','sum'),ASDR_upper_weighted_marginal=('weighted_upper','sum'),n_ages=('age_id','nunique')).reset_index()
a['PAF_percent']=100*a.attributable_DALYs/a.total_glaucoma_DALYs
a['crude_rate']=a.attributable_DALYs/a.population*1e5
ht=trends(a,['location_id','sex_id'],'ASDR');a=a.merge(ht,on=['location_id','sex_id'],validate='many_to_one')
a['cause_id']=670;a['rei_id']=105;a['measure_id']=2
a=a.sort_values(['location_id','year'])
age2021=h[h.year.eq(2021)].copy()
age2021=age2021.rename(columns={'val_1':'attributable_DALYs','lower_1':'attributable_DALYs_lower_original_95UI','upper_1':'attributable_DALYs_upper_original_95UI','val_3':'age_specific_rate','lower_3':'rate_lower_original_95UI','upper_3':'rate_upper_original_95UI'})
age2021['cause_id']=670;age2021['rei_id']=105;age2021['measure_id']=2
ac=['location_name','location_id','sex_id','year','age_group','age_id','cause_id','rei_id','measure_id','attributable_DALYs','attributable_DALYs_lower_original_95UI','attributable_DALYs_upper_original_95UI','population','age_specific_rate','rate_lower_original_95UI','rate_upper_original_95UI','total_glaucoma_DALYs']
age2021=age2021[ac].sort_values(['location_id','age_id'])
assert len(age2021)==54 and len(a)==192
age2021.to_csv(OUT/'tables/TableS4_HFPG_Age_specific_2021.csv',index=False)
a.to_csv(OUT/'tables/TableS5_HFPG_Annual_1990_2021.csv',index=False)
wide=age2021.assign(value=age2021.apply(lambda r:f'{r.age_specific_rate:.2f} ({r.rate_lower_original_95UI:.2f}, {r.rate_upper_original_95UI:.2f})',axis=1)).pivot(index='age_id',columns='location_name',values='value').loc[AGES,names.loc[LOCS].tolist()]
wide.index=w.age_name;wide.to_csv(OUT/'tables/TableS4_HFPG_display_9x6.csv')
payload={**country_tables,'HFPG_Age_specific_2021':age2021,'HFPG_Annual_1990_2021':a}

map_values=[]
for m,label in [(5,'Prevalence'),(2,'DALYs')]:
    z=country_tables[f'Country_{label}'].query('sex_id==3')
    for value in ['number_2021','ASR_2021','EAPC']:
        map_values.append(z[['location_id','location_name',value]].rename(columns={value:'value'}).assign(measure_id=m,metric=value))
pd.concat(map_values).to_csv(OUT/'tables/Figure5_six_panel_source.csv',index=False)
six=a[a.year.eq(1990)].merge(a[a.year.eq(2021)],on=['location_id','location_name','sex_id'],suffixes=('_1990','_2021'))
six['ASDR_change_percent']=100*(six.ASDR_2021/six.ASDR_1990-1)
six.to_csv(OUT/'tables/HFPG_six_group_endpoint_and_EAPC_summary.csv',index=False)
for label in ['Prevalence','DALYs']:
    z=country_tables[f'Country_{label}'].query('sex_id==3').sort_values('ASR_2021',ascending=False)
    z['rank_descending']=np.arange(1,len(z)+1)
    z.to_csv(OUT/f'tables/Country_{label}_ranking_2021.csv',index=False)
