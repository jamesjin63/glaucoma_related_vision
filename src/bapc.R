suppressPackageStartupMessages({library(data.table); library(BAPC);library(INLA)})
inla.setOption(num.threads="2:1")
O <- Sys.getenv("GLAUCOMA_OUTPUT"); input <- Sys.getenv("GLAUCOMA_INPUT"); dir.create(file.path(O,"data/bapc_models"),showWarnings=FALSE)
w <- fread(file.path(O,"data/verified_age_specific_1990_2021.csv"))
aw <- fread(file.path(O,"tables/TableS0_standard_weights_55plus.csv")); ages <- aw$age_id
an <- c("55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+")
ss <- fread(file.path(O,"tables/TableS_verified_summary_1990_2021.csv"))
p <- fread(file.path(input,"support/population_support_global_2011_2050.csv"))
p[,sex_id:=match(sex_name,c("Male","Female","Both"))]
obs <- unique(w[location_id==1 & measure_id==5,.(sex_id,age_id,year,population)])
sc <- merge(obs[year==2021],p[year==2021],by=c("sex_id","age_id","year"))[,.(sex_id,age_id,scale=population/val)]
fwrite(sc,file.path(O,"tables/Diagnostic_future_population_baseline_scaling.csv"))
fut <- merge(p[year>2021],sc,by=c("sex_id","age_id"))[,.(sex_id,age_id,year,population=val*scale)]
pop <- rbind(obs,fut);fwrite(pop,file.path(O,"tables/TableS_population_global_1990_2050.csv"))
mat <- function(d,v) {
 z <- dcast(d,year~age_id,value.var=v);setorder(z,year)
 m <- as.data.frame(z)[,as.character(ages),drop=FALSE];rownames(m)<-z$year;colnames(m)<-an;m
}
run <- function(cs,pm,np,period="rw2",b=.00005){
 BAPC(APCList(cs,pm,gf=5),predict=list(npredict=np,retro=TRUE),stdweight=aw$weight,verbose=FALSE,
  model=list(age=list(model="rw2",prior="loggamma",param=c(1,b)),
   period=list(include=TRUE,model=period,prior="loggamma",param=c(1,b)),
   cohort=list(include=TRUE,model="rw2",prior="loggamma",param=c(1,b)),
   overdis=list(include=TRUE,model="iid",prior="loggamma",param=c(1,.005))))
}
extract <- function(res,yrs,pm){
 res<-qapc(res,percentiles=c(.025,.975))
 a<-as.data.table(agestd.rate(res));setnames(a,c("ASR","ASR_sd","ASR_lower","ASR_upper"));a[,year:=yrs]
 a[,c("ASR","ASR_sd","ASR_lower","ASR_upper"):=lapply(.SD,`*`,1e5),.SDcols=c("ASR","ASR_sd","ASR_lower","ASR_upper")]
 ar<-rbindlist(lapply(seq_along(an),function(i){z<-as.data.table(agespec.rate(res)[[i]]);setnames(z,c("rate","sd","rate_lower","rate_upper"));z[,`:=`(year=yrs,age_id=ages[i])];z}))
 pr<-melt(as.data.table(cbind(year=yrs,as.matrix(pm))),id.vars="year",variable.name="age_name",value.name="population")
 pr[,age_id:=ages[match(age_name,an)]]
 ar<-merge(ar,pr,by=c("year","age_id"));ar[,`:=`(number=rate*population,number_lower_approx=rate_lower*population,number_upper_approx=rate_upper*population)]
 n<-ar[,.(number=sum(number),number_lower_approx=sum(number_lower_approx),number_upper_approx=sum(number_upper_approx)),by=year]
 ar[,c("rate","sd","rate_lower","rate_upper"):=lapply(.SD,`*`,1e5),.SDcols=c("rate","sd","rate_lower","rate_upper")]
 list(total=merge(a,n,by="year"),age=ar)
}
TOT<-list();AGE<-list();VAL<-list();SENS<-list()
for(m in c(5,2))for(s in c(3,1,2)){
 key<-data.table(measure_id=m,sex_id=s);keyname<-paste(m,s,sep="_")
 d<-w[location_id==1 & measure_id==m & sex_id==s];cs<-mat(d,"val_1");cs<-round(cs)
 pm<-mat(pop[sex_id==s],"population")
 cs2<-rbind(cs,matrix(NA,29,9,dimnames=list(2022:2050,an)))
 cache<-file.path(O,"data/bapc_models",paste0(keyname,".rds"))
 if(file.exists(cache))res<-readRDS(cache)else{res<-run(cs2,pm,29);saveRDS(res,cache)}
 x<-extract(res,1990:2050,pm)
 truth<-ss[location_id==1 & measure_id==m & sex_id==s,.(year,observed_ASR=ASR,observed_number=number)]
 TOT[[keyname]]<-cbind(key,merge(x$total,truth,by="year",all.x=TRUE));AGE[[keyname]]<-cbind(key,x$age)
 vcache<-file.path(O,"data/bapc_models",paste0(keyname,"_validation.rds"))
 pmobs<-pm[as.character(1990:2021),]
 if(file.exists(vcache))rv<-readRDS(vcache)else{rv<-run(cs,pmobs,10);saveRDS(rv,vcache)}
 v<-merge(extract(rv,1990:2021,pmobs)$total,truth,by="year")[year>=2012]
 # Comparator uses fixed age-specific rates in 2011, with known held-out populations.
 last<-d[year==2011,.(age_id,last_rate=val_3)]
 naive<-merge(pop[sex_id==s & year %in% 2012:2021],last,by="age_id")
 naive<-merge(naive,aw[,.(age_id,weight)],by="age_id")[,.(naive_ASR=sum(weight*last_rate),naive_number=sum(last_rate*population/1e5)),by=year]
 VAL[[keyname]]<-cbind(key,merge(v,naive,by="year"))
 if(s==3){
  for(pr in c("rw1","half_prior_scale","double_prior_scale")){
   pc<-if(pr=="rw1")"rw1"else"rw2";b<-switch(pr,rw1=.00005,half_prior_scale=.000025,double_prior_scale=.0001)
   scache<-file.path(O,"data/bapc_models",paste0(keyname,"_",pr,".rds"))
   if(file.exists(scache))rs<-readRDS(scache)else{rs<-run(cs2,pm,29,pc,b);saveRDS(rs,scache)}
   se<-extract(rs,1990:2050,pm)$total[year==2050]
   SENS[[paste(m,pr)]]<-cbind(key,data.table(scenario=pr),se)
  }
 }
 cat("Completed measure",m,"sex",s,"\n");flush.console()
}
tot<-rbindlist(TOT);ag<-rbindlist(AGE);val<-rbindlist(VAL);sen<-rbindlist(SENS)
fwrite(tot,file.path(O,"tables/Table7_BAPC_global_annual_1990_2050.csv"))
fwrite(ag,file.path(O,"tables/TableS_BAPC_age_specific_1990_2050.csv"))
fwrite(val,file.path(O,"tables/TableS_temporal_validation_yearly_2012_2021.csv"))
metrics<-val[,.(n= .N,MAPE_ASR=100*mean(abs(ASR/observed_ASR-1)),bias_ASR=100*mean(ASR/observed_ASR-1),coverage_ASR=100*mean(observed_ASR>=ASR_lower & observed_ASR<=ASR_upper),MAPE_number=100*mean(abs(number/observed_number-1)),coverage_number_approx=100*mean(observed_number>=number_lower_approx & observed_number<=number_upper_approx),naive_MAPE_ASR=100*mean(abs(naive_ASR/observed_ASR-1)),naive_MAPE_number=100*mean(abs(naive_number/observed_number-1))),by=.(measure_id,sex_id)]
fwrite(metrics,file.path(O,"tables/Table8_temporal_validation_fit1990_2011_predict2012_2021.csv"))
base<-tot[year==2050 & sex_id==3];base[,scenario:="baseline_rw2"]
sen<-rbindlist(list(base,sen),fill=TRUE)
sen<-merge(sen,base[,.(measure_id,baseline_ASR=ASR,baseline_number=number)],by="measure_id")
sen[,`:=`(ASR_vs_baseline_pct=100*(ASR/baseline_ASR-1),number_vs_baseline_pct=100*(number/baseline_number-1))]
fwrite(sen,file.path(O,"tables/TableS_BAPC_prior_and_period_sensitivity_2050.csv"))
fwrite(tot[year==2050],file.path(O,"tables/Table9_projection_2050.csv"))
capture.output(sessionInfo(),file=file.path(O,"data/R_sessionInfo.txt"))
print(metrics);print(tot[year==2050])
