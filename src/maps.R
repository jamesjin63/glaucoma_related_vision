try(Sys.setlocale("LC_CTYPE","en_US.UTF-8"),silent=TRUE)
suppressPackageStartupMessages({library(data.table);library(sf);library(ggplot2);library(patchwork);library(scales)})
dest<-commandArgs(trailingOnly=TRUE)[1]
w<-st_read(file.path(Sys.getenv("GLAUCOMA_INPUT"),"support/world_map_support.geojson"),quiet=TRUE)
s<-fread(paste0(Sys.getenv("GLAUCOMA_OUTPUT"),"/tables/TableS_verified_summary_1990_2021.csv"))[level=="Country" & year==2021 & sex_id==3]
e<-fread(paste0(Sys.getenv("GLAUCOMA_OUTPUT"),"/tables/TableS_EAPC_ASR_vs_crude_and_HAC.csv"))[sex_id==3]
s<-merge(s,e,by=c("location_id","measure_id","sex_id"))
jobs<-data.table(m=c(5,5,2,2,5,2),v=c("number","ASR","number","ASR","EAPC","EAPC"),label=c("Prevalent cases, 2021","ASPR per 100,000, 2021","DALYs, 2021","ASDR per 100,000, 2021","ASPR EAPC (%), 1990–2021","ASDR EAPC (%), 1990–2021"))
plots<-list()
for(i in 1:6){
 z<-s[measure_id==jobs$m[i],.(location_id,value=get(jobs$v[i]))]
 d<-merge(w,z,by="location_id",all.x=TRUE)
 p<-ggplot(d)+geom_sf(aes(fill=value),colour="white",linewidth=.1)+coord_sf(ylim=c(-56,84),expand=FALSE)+theme_void(base_size=12)+labs(title=paste(LETTERS[i],jobs$label[i]),fill=NULL)+theme(plot.title=element_text(face="bold",size=12),legend.position="bottom",legend.key.width=unit(1.1,"cm"))
 if(jobs$v[i]=="EAPC")p<-p+scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,na.value="grey85")else p<-p+scale_fill_viridis_c(option="C",trans=if(jobs$v[i]=="number")"log10"else"identity",labels=label_number(scale_cut=cut_short_scale()),na.value="grey85")
 plots[[i]]<-p
 ggsave(file.path(dest,paste0("map_panel_",LETTERS[i],".pdf")),p,width=11,height=5.8,device=cairo_pdf,bg="white")
}
overview<-wrap_plots(plots,ncol=2)+plot_annotation(caption="Adults aged 55+, both sexes, GBD 2021. Number colour scales are logarithmic; grey denotes no mapped estimate. Large panels follow in the PDF.")
ggsave(file.path(dest,"map_overview.pdf"),overview,width=17,height=14,device=cairo_pdf,bg="white")
ggsave(file.path(dest,"Fig5_country_maps.png"),overview,width=17,height=14,dpi=300,device=grDevices::png,type="cairo",bg="white")
