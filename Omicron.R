library(tidyverse)
library(gtools)
library(broom)
library(gee)

#require(foreign)
require(nnet)
#require(ggplot2)
#require(reshape2)

setwd("/Users/cc19/Documents/GRUBAUGH_LAB/OMICRON")

## Function for plotting output from the regression model

plot.odds.ratio<-function(glm.model,log.coef=FALSE,as.is=FALSE){
  if(isTRUE(as.is)){
    input.model<-glm.model %>% dplyr::filter(term!="(Intercept)")
    log.coef<-FALSE
  }else{
    if(isTRUE(log.coef)){
      input.model<-glm.model %>% tidy(exponentiate=FALSE,conf.int=TRUE) %>%
        dplyr::filter(term!="(Intercept)")
    }else{
      input.model<-glm.model %>% tidy(exponentiate=TRUE,conf.int=TRUE) %>%
        dplyr::filter(term!="(Intercept)")
    }
  }
  
  if(isTRUE(log.coef)){
    vals.x<-pretty(c(min(as.double(input.model$conf.low))-abs(min(as.double(input.model$conf.high))-max(as.double(input.model$conf.low) ))/2,
                     max(as.double(input.model$conf.high) )+abs(min(as.double(input.model$conf.high))-max(as.double(input.model$conf.low) ))/2))
  }else{
    vals.x<-pretty(c(0,max(as.double(input.model$conf.high) )+abs(min(as.double(input.model$conf.high))-max(as.double(input.model$conf.low) ))/2))
  }
  
  vals.y<-c(0,seq(1,length(unique(input.model$term)) )-0.5,length(unique(input.model$term))-0.5+0.25)
  plot(0,0,ylim=c(0,max(vals.y)),col=rgb(1,1,1,alpha=0.01),lwd=1.5,
       xlim=c(min(vals.x),max(vals.x)),las=1,yaxt="n",xaxs="i",yaxs="i",bty="n",
       xlab=NA,ylab=NA,main=NA)
  axis(2,labels=c("",input.model$term,""),at=vals.y,las=1)
  if(isTRUE(log.coef)){
    mtext(bquote("Coefficient ("~beta~")"),side=1,line=2.5)
  }else{
    mtext(bquote("Odds ratio ("~beta~")"),side=1,line=2.5)
  }
  mtext(bquote("Covariate"),side=2,line=7)
  
  if(isTRUE(log.coef)){
    segments(0,0,0,max(length(unique(input.model$term)))+1.0,lty=4)
  }else{
    segments(1,0,1,max(length(unique(input.model$term)))+1.0,lty=4)
  }
  segments(min(vals.x),0,min(vals.x),max(length(unique(input.model$term)))+1.0,lwd=1.5)
  segments(min(vals.x),0,max(vals.x),0,lwd=1.5)
  
  for(i in 1:length(input.model$term)){
    j<-i-0.5
    points(as.double(input.model$estimate[i]),j,pch=20,cex=2)
    segments(as.double(input.model$conf.low[i]),j,as.double(input.model$conf.high[i]),j)
    #text(max(as.double(input.model$estimate)),i,bquote(italic("p")~"="~.(input.model$p.value[i])),cex=0.75)
    
    pval.text=ifelse(as.double(input.model$p.value[i])>=0.05, "",
                     ifelse(as.double(input.model$p.value[i])>0.01,"*",
                            ifelse(as.double(input.model$p.value[i])>0.001,"**",
                                   ifelse(as.double(input.model$p.value[i])<0.001,"***","")) ))
    
    if(pval.text!=""){
      text(max(vals.x)-1,j,bquote(.(pval.text)),cex=1.00)
      ##text(max(as.double(input.model$estimate)),i,bquote(italic("p")~"="~.(pval.text)),cex=0.75)
    }
  }
}

############ Read and prepare the data for the logistic regression analysis

Om.data.new<-as_tibble(data.table::fread("./DATA_FROM_ANDREAS/final_ynhh.tsv",header=TRUE)) %>%
  dplyr::select(hash_subject_id,sex,age)  %>% group_by(hash_subject_id) %>%
  add_count() %>% arrange(desc(n),hash_subject_id) %>% 
  mutate(age=ifelse(n>1,last(age),age)) %>% unique() %>% dplyr::select(-n)
  ##column_to_rownames("hash_subject_id")

## Load the cleaned data from N.G. and merge with the latest data from YNHH

Om.data<-as_tibble(data.table::fread("taqpath_12122021-12262021_2022-01-11.csv",header=TRUE)) %>%
  full_join(Om.data.new,by="hash_subject_id")

## Save the merged dataset

Om.data %>% write.table("Omicron_YNHH_merged.tsv",row.names=FALSE,sep="\t")
  
## Clean and prepared the merged dataset for the test positivity and regression analysis to show odds of detecting Omicron relative to Delta among infected individuals 

Om.data.corr<-Om.data %>% mutate(result=ifelse(ORF1ab_Ct>30,"negative",result)) %>% #rowwise() %>%
  dplyr::rename(S_Ct=`S-gene_Ct`,test_result=result) %>% 
  mutate(measurement_date=as.Date(measurement_date,"%m/%d/%y"),
         vax_2_date=as.Date(vax_2_date,"%m/%d/%y"),
         vax_1_date=as.Date(vax_1_date,"%m/%d/%y"),
         booster_1_date=as.Date(booster_1_date,"%m/%d/%y"),
         Doses=as.integer(Doses)) %>% #rowwise() %>%
  mutate(dose1=ifelse(Doses==1,1,0),
         age=ifelse(age==">=90","90",age),
         sex=ifelse(sex=="","Unknown",sex),
         days_from_2nd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_2_date,"%m/%d/%y")),
         days_from_1st_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_1_date,"%m/%d/%y")),
         days_from_3rd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(booster_1_date,"%m/%d/%y"))) %>%
  dplyr::filter(test_result!="negative") %>%
  mutate(vax_1_mfr=ifelse(vax_1_mfr=="J&J","1",
                          ifelse(vax_1_mfr=="Pfizer","2",
                          ifelse(vax_1_mfr=="Moderna","3","")))) %>%
  mutate(vax_2_mfr=ifelse(vax_2_mfr=="J&J","1",
                          ifelse(vax_2_mfr=="Pfizer","2",
                          ifelse(vax_2_mfr=="Moderna","3","")))) %>%
  mutate(booster_1_mfr=ifelse(booster_1_mfr=="J&J","1",
                              ifelse(booster_1_mfr=="Pfizer","2",
                              ifelse(booster_1_mfr=="Moderna","3","")))) %>%
  mutate(test_result=ifelse(test_result=="Omicron",1,0)) %>% 
  mutate(vax_status=ifelse(Doses==0,"0",
                           ifelse(Doses==1,"1",
                           ifelse(Doses==2 & days_from_2nd_dose/30<5,"2",
                           ifelse(Doses==2 & days_from_2nd_dose/30>=5,"3",
                           ifelse(Doses>=3,"4","Unknown")))))) %>%
  mutate(vax_status=fct_recode(vax_status,"unvax"="0","vax_1d"="1","vax_2d_5m"="2","vax_2d_5m_gt"="3","vax_3d"="4","Unknown"="Unknown")) %>%
  mutate(vax_status=fct_relevel(vax_status,"unvax",after=0)) %>% rowwise() %>%
  mutate(vax_status_=ifelse(Doses==0,"unvax",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="1","jj_1d_5m",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="1","jj_1d_5m_gt",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="2","pf_1d_5m",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="2","pf_1d_5m_gt",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="3","mod_1d_5m",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="3","mod_1d_5m_gt",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="1","jj_2d_5m",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="1","jj_2d_5m_gt",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m_gt",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m",
                            ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m_gt",
                            ifelse(Doses>=3,"vax_3d","unknown") )))))))))))))) %>% ungroup() %>%
  mutate(Doses=ifelse(Doses==3 & as.integer(measurement_date-booster_1_date)<=14,Doses-1,
                      ifelse(Doses==2 & as.integer(measurement_date-vax_2_date)<=14,Doses-1,
                             ifelse(Doses==1 & as.integer(measurement_date-vax_1_date)<=14,Doses-1,Doses)  ) )) %>%
  dplyr::filter(sex!="Unknown" & vax_status!="Unknown" & vax_status_!="Unknown") %>% 
  mutate(sex=as.factor(sex)) %>%
  group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
  mutate(vax_status_=fct_relevel(vax_status_,"unvax",after=0),
         sex=fct_relevel(sex,"Female",after=0),
         ORF1ab_Ct=as.double(ORF1ab_Ct),age=as.double(age)) %>%
  group_by(hash_subject_id) %>% arrange(measurement_date) %>% dplyr::filter(row_number()==1)

Om.data.corr$vax_status <- factor(Om.data.corr.ct$vax_status, levels = c("unvax", "vax_1d", "vax_2d_5m", "vax_2d_5m_gt", "vax_3d"))
Om.data.corr$vax_status_ <- factor(Om.data.corr.ct$vax_status_, levels = c("unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d"))


table(Om.data.corr$test_result,Om.data.corr$vax_status_) %>%
  write.table("Omicron_Delta_positives.tsv",sep="\t",row.names=FALSE)

############ Fit the logistic regression models

glm.out1<-glm(test_result~vax_status+age+sex,data=Om.data.corr,family=binomial(link="logit"),model=TRUE)
summary(glm.out1)
glm.out1 %>% tidy(exponentiate=TRUE,conf.int=TRUE)

glm.out2<-glm(test_result~vax_status_+age+sex,data=Om.data.corr,family=binomial(link="logit"),model=TRUE)
summary(glm.out2)
glm.out2 %>% tidy(exponentiate=TRUE,conf.int=TRUE)

##############
## Prepare data for the regression analysis of the PCR Ct values
  
Om.data.corr.ct<-Om.data %>% mutate(result=ifelse(ORF1ab_Ct>30,"negative",result)) %>% #rowwise() %>%
  dplyr::rename(S_Ct=`S-gene_Ct`,test_result=result) %>% 
  mutate(measurement_date=as.Date(measurement_date,"%m/%d/%y"),
         vax_2_date=as.Date(vax_2_date,"%m/%d/%y"),
         vax_1_date=as.Date(vax_1_date,"%m/%d/%y"),
         booster_1_date=as.Date(booster_1_date,"%m/%d/%y"),
         Doses=as.integer(Doses)) %>% #rowwise() %>%
  mutate(dose1=ifelse(Doses==1,1,0),
         age=ifelse(age==">=90","90",age),
         sex=ifelse(sex=="","Unknown",sex),
         days_from_2nd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_2_date,"%m/%d/%y")),
         days_from_1st_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_1_date,"%m/%d/%y")),
         days_from_3rd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(booster_1_date,"%m/%d/%y"))) %>%
  dplyr::filter(vax_1_mfr!="unknown" | vax_2_mfr!="unknown" | booster_1_mfr!="unknown") %>%
  dplyr::filter(test_result!="negative") %>%
  mutate(vax_1_mfr=ifelse(vax_1_mfr=="J&J","1",
                          ifelse(vax_1_mfr=="Pfizer","2",
                                 ifelse(vax_1_mfr=="Moderna","3","")))) %>%
  mutate(vax_2_mfr=ifelse(vax_2_mfr=="J&J","1",
                          ifelse(vax_2_mfr=="Pfizer","2",
                                 ifelse(vax_2_mfr=="Moderna","3","")))) %>%
  mutate(booster_1_mfr=ifelse(booster_1_mfr=="J&J","1",
                              ifelse(booster_1_mfr=="Pfizer","2",
                                     ifelse(booster_1_mfr=="Moderna","3","")))) %>%
  mutate(test_result=ifelse(test_result=="Omicron",1,0)) %>% 
  mutate(vax_status=ifelse(Doses==0,"0",
                           ifelse(Doses==1,"1",
                                  ifelse(Doses==2 & days_from_2nd_dose/30<5,"2",
                                         ifelse(Doses==2 & days_from_2nd_dose/30>=5,"3",
                                                ifelse(Doses>=3,"4","Unknown")))))) %>%
  mutate(vax_status=fct_recode(vax_status,"unvax"="0","vax_1d"="1","vax_2d_5m"="2","vax_2d_5m_gt"="3","vax_3d"="4","Unknown"="Unknown")) %>%
  mutate(vax_status=fct_relevel(vax_status,"unvax",after=0)) %>% rowwise() %>%
  ##dplyr::filter(vax_1_mfr==vax_2_mfr) %>%
  mutate(vax_status_=ifelse(Doses==0,"unvax",
                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="1","jj_1d_5m",
                                   ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="1","jj_1d_5m_gt",
                                          ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="2","pf_1d_5m",
                                                 ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="2","pf_1d_5m_gt",
                                                        ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="3","mod_1d_5m",
                                                               ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="3","mod_1d_5m_gt",
                                                                      ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="1","jj_2d_5m",
                                                                             ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="1","jj_2d_5m_gt",
                                                                                    ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m",
                                                                                           ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m_gt",
                                                                                                  ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m",
                                                                                                         ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m_gt",
                                                                                                                ifelse(Doses>=3,"vax_3d","unknown") )))))))))))))) %>% ungroup() %>%
  mutate(Doses=ifelse(Doses==3 & as.integer(measurement_date-booster_1_date)<=14,Doses-1,
                      ifelse(Doses==2 & as.integer(measurement_date-vax_2_date)<=14,Doses-1,
                             ifelse(Doses==1 & as.integer(measurement_date-vax_1_date)<=14,Doses-1,Doses)  ) )) %>%
  dplyr::filter(sex!="Unknown" & vax_status!="Unknown" & vax_status_!="Unknown") %>% 
  mutate(sex=as.factor(sex)) %>%
  dplyr::filter(!grepl(">",ORF1ab_Ct)) %>% 
  mutate(ORF1ab_Ct=as.double(ORF1ab_Ct)) %>%
  group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
  mutate(vax_status_=fct_relevel(vax_status_,"unvax",after=0),
         sex=fct_relevel(sex,"Female",after=0),
         ORF1ab_Ct=as.double(ORF1ab_Ct),age=as.double(age)) %>%
  group_by(hash_subject_id) %>% arrange(measurement_date) %>% dplyr::filter(row_number()==1)

Om.data.corr.ct$vax_status <- factor(Om.data.corr.ct$vax_status, levels = c("unvax", "vax_1d", "vax_2d_5m", "vax_2d_5m_gt", "vax_3d"))
Om.data.corr.ct$vax_status_ <- factor(Om.data.corr.ct$vax_status_, levels = c("unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d"))

## Plot a figure showing the PCR Ct values and regression output showing covariates associated with viral load

{
  pdf(file="PCR_Ct_figure.pdf",height=8.5,width=10.5)
  layout(matrix(c(1,2,3,4),byrow=TRUE,ncol=2))
  
  mai.mat<-c(0.751,1.751,0.35,0.25)
  
  {
    par(bty="n",mai=mai.mat)
    tmp.Om<-Om.data %>% mutate(test_result=ifelse(ORF1ab_Ct>30,"negative",result)) %>%
      dplyr::filter(result %in% c("Delta","Omicron")) %>%
      mutate(ORF1ab_Ct=as.double(ORF1ab_Ct)) %>% dplyr::select(ORF1ab_Ct,Doses,result)
    tmp.data<-list("ALL"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta"],
                   "0"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==0],
                   "1"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==1],
                   "2"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==2],
                   "3"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==3])
    
    tmp.data1<-list("ALL"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron"],
                    "0"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==0],
                    "1"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==1],
                    "2"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==2],
                    "3"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==3])
    
    y.vals=pretty(c(5,35))
    vioplot::vioplot(tmp.data,
                     las=1, plotCentre = "line", side = "left",lwd=0.001,
                     col=viridis::inferno(5)[c(4)], border="white", rectCol="palevioletred", lineCol="violetred", colMed="violet",
                     xlab="Vaccination status",ylab=bquote("Nasal swab PCR C"[t]~"-value"),yaxt="n",xaxt="n")
    axis(2,at=y.vals,labels=y.vals,las=1)
    axis(1,at=1:length(names(tmp.data)),labels=names(tmp.data),las=1)
    
    vioplot::vioplot(tmp.data1,
                     las=1, plotCentre = "line", side = "right",lwd=0.001,add = TRUE,
                     col=viridis::inferno(5)[c(3)], border="white", rectCol="palevioletred", lineCol="violetred", colMed="violet",
                     xlab="Vaccination status",ylab=bquote("Nasal swab PCR C"[t]~"-value"),yaxt="n",xaxt="n")
    
    legend("topleft",legend=c("Delta","Omicron"),#title="Doses",
           cex=0.75,
           fill=c(viridis::inferno(5)[c(4,3)]),bty="n")
  }
  
  {
    plot.new()
  }
  
  {
    par(mai=mai.mat)
    glm.out.ct.glm<-glm(ORF1ab_Ct~test_result+vax_status+age+sex,data=Om.data.corr.ct,family=gaussian,model=TRUE)
    glm.out.ct<-as.data.frame(coef(summary(glm.out.ct.glm)))

    glm.out<-glm.out.ct %>% rownames_to_column("term") %>% 
      mutate(conf.low=Estimate-1.96*`Std. Error`,
             conf.high=Estimate+1.96*`Std. Error`,
             std.error=`Std. Error`,
             p.value=`Pr(>|t|)`,statistic=`t value`) %>%
      mutate(estimate=exp(Estimate),conf.low=exp(conf.low),conf.high=exp(conf.high)) %>% #arrange(term) %>%
      dplyr::select(term,estimate,std.error,statistic,p.value,conf.low,conf.high)
    
    plot.odds.ratio(glm.out,as.is=TRUE)
    title(bquote("Nasal swab PCR C"[t]~"-value"),font.main=1,cex.main=1.0,line=1)
    
    glm.out %>% dplyr::filter(term!="(Intercept)") %>% 
      mutate(across(where(is.numeric), round, digits=2)) %>%
      mutate(estimate=paste0(estimate," (",conf.low,",",conf.high,")")) %>% #arrange(term) %>%
      dplyr::select(term,estimate,p.value) %>% write.table(file="Om_delta.OR_.tsv",sep="\t",row.names=FALSE)
  }
  
  {
    par(mai=mai.mat)
    glm.out.ct.glm<-glm(ORF1ab_Ct~test_result+vax_status_+age+sex,data=Om.data.corr.ct,family=gaussian,model=TRUE)
    glm.out.ct<-as.data.frame(coef(summary(glm.out.ct.glm)))

    glm.out<-glm.out.ct %>% rownames_to_column("term") %>% 
      mutate(conf.low=Estimate-1.96*`Std. Error`,
             conf.high=Estimate+1.96*`Std. Error`,
             std.error=`Std. Error`,
             p.value=`Pr(>|t|)`,statistic=`t value`) %>%
      mutate(estimate=exp(Estimate),conf.low=exp(conf.low),conf.high=exp(conf.high)) %>% #arrange(term) %>%
      dplyr::select(term,estimate,std.error,statistic,p.value,conf.low,conf.high)
    
    plot.odds.ratio(glm.out,as.is=TRUE)
    title(bquote("Nasal swab PCR C"[t]~"-value"),font.main=1,cex.main=1.0,line=1)
    
    glm.out %>% dplyr::filter(term!="(Intercept)") %>% 
      mutate(across(where(is.numeric), round, digits=4)) %>%
      mutate(estimate=paste0(estimate," (",conf.low,",",conf.high,")")) %>% #arrange(term) %>%
      dplyr::select(term,estimate,p.value) %>% write.table(file="Om_delta.OR__.tsv",sep="\t",row.names=FALSE)
  }
  dev.off()
}


# Plot the number of cases, test positivity rates, and regression estimates for the odds of infection with Omicron vs Delta in infected individuals

{
  pdf(file="Figure1_draft_.pdf",height=8.0,width=9.05)
  layout(matrix(c(1,1,2,3,4,5),ncol=2,byrow=TRUE),heights=c(1,1,1))
  mai.vals<-c(0.761,0.751,0.061,0.10)
  mai.vals1<-c(0.761,0.751,0.061,0.650)
  
  {
    par(mai=mai.vals1 )
    
    dat_freq_both_emerge<-as_tibble(readRDS("dat_freq_both_emerge.rds"))
    dat_freq_both_emerge$Date<-as.Date(dat_freq_both_emerge$Date,"%d/%m/%y")
    dat_freq_both_emerge
    
    dat_cases_long<-as_tibble(readRDS("dat_cases_long.rds")) %>% group_by(Variant) %>%
      mutate(Counter=1:n())
    dat_cases_long$Date<-as.Date(dat_cases_long$Date,"%d/%m/%y")
    
    variant.cases<-dat_cases_long %>% dplyr::select(Date,Variant,Num_Cases) %>%
      spread(key=Variant,value=Num_Cases) %>% drop_na() %>% column_to_rownames("Date") %>% as.matrix()
    
    
    variant.cases1<-variant.cases/rowSums(variant.cases)
    
    prop.Om<-dat_cases_long %>% dplyr::select(Date,Variant,prop_omicron_total) %>%
      drop_na() %>% dplyr::filter(Variant=="Omicron")
    
    ylim.vals<-pretty(c(0,max(variant.cases,na.rm=TRUE)+200))
    xlim.vals<-pretty(c(min(unique(dat_cases_long$Counter)),max(unique(dat_cases_long$Counter))),n=40)
    x <- barplot(t(variant.cases),beside=FALSE,las=1,col=c(viridis::inferno(5)[4],viridis::inferno(5)[3]),
                 ylim=c(0,max(ylim.vals)),ylab="Positive cases",xaxs="i",
                 xlab="Test date",xaxt="n")
    
    tmpW<-dat_cases_long$Date[dat_cases_long$Counter %in% xlim.vals] %>% as_tibble() %>% 
      mutate(XX=paste0(lubridate::month(value),"/",lubridate::day(value),"/",gsub("2022","21",gsub("2021","21",lubridate::year(value))))) %>% 
      dplyr::select(XX) %>% as.data.frame()
    text(cex=0.95,x=x[dat_cases_long$Counter %in% xlim.vals],y=-3,
         tmpW$XX,xpd=TRUE,srt=45,adj=1.10)
    par(new=TRUE)
    plot(prop.Om$Date,prop.Om$prop_omicron_total*max(ylim.vals),las=1,yaxs="i",xaxs="i",bty="n",
         ylim=c(0,max(ylim.vals)),ylab="Positive cases",xaxt="n",xlab=NA,type="p",lwd=2.0,pch=19,
         col=viridis::inferno(5)[3])
    lines(smooth.spline(x=prop.Om$Date, 
                        y=prop.Om$prop_omicron_total*max(ylim.vals),df=6),col=viridis::inferno(5)[3],lwd=2)
    axis(4,at=seq(0,1,length.out=6)*max(ylim.vals),labels=seq(0,1,length.out=6),las=1,xaxs="i")
    
    rect(as.Date("12/12/21","%d/%m/%y"),0,as.Date("26/12/21","%d/%m/%y"),1200)
    
    mtext("%SGTF",side=4,line=3,cex=0.85)
    legend("topleft",legend=c("Delta","Omicron"),
           fill=c(viridis::inferno(5)[4],viridis::inferno(5)[3]),
           bty="n",cex=0.80) #,title="Variant"
  }
  
  
  
  Om.data.tmp<-Om.data %>% mutate(result=ifelse(ORF1ab_Ct>30,"negative",result)) %>% #rowwise() %>%
    dplyr::rename(S_Ct=`S-gene_Ct`,test_result=result) %>% 
    mutate(measurement_date=as.Date(measurement_date,"%m/%d/%y"),
           vax_2_date=as.Date(vax_2_date,"%m/%d/%y"),
           vax_1_date=as.Date(vax_1_date,"%m/%d/%y"),
           booster_1_date=as.Date(booster_1_date,"%m/%d/%y"),
           Doses=as.integer(Doses)) %>% #rowwise() %>%
    mutate(dose1=ifelse(Doses==1,1,0),
           age=ifelse(age==">=90","90",age),
           sex=ifelse(sex=="","Unknown",sex),
           days_from_2nd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_2_date,"%m/%d/%y")),
           days_from_1st_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_1_date,"%m/%d/%y")),
           days_from_3rd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(booster_1_date,"%m/%d/%y"))) %>%
    #dplyr::filter(vax_1_mfr!="unknown" | vax_2_mfr!="unknown" | booster_1_mfr!="unknown") %>%
    ##dplyr::filter(test_result!="negative") %>%
    mutate(vax_1_mfr=ifelse(vax_1_mfr=="J&J","1",
                            ifelse(vax_1_mfr=="Pfizer","2",
                                   ifelse(vax_1_mfr=="Moderna","3","")))) %>%
    mutate(vax_2_mfr=ifelse(vax_2_mfr=="J&J","1",
                            ifelse(vax_2_mfr=="Pfizer","2",
                                   ifelse(vax_2_mfr=="Moderna","3","")))) %>%
    mutate(booster_1_mfr=ifelse(booster_1_mfr=="J&J","1",
                                ifelse(booster_1_mfr=="Pfizer","2",
                                       ifelse(booster_1_mfr=="Moderna","3","")))) %>%
    mutate(test_result=ifelse(test_result=="Omicron",1,0)) %>% 
    mutate(vax_status=ifelse(Doses==0,"0",
                             ifelse(Doses==1,"1",
                                    ifelse(Doses==2 & days_from_2nd_dose/30<5,"2",
                                           ifelse(Doses==2 & days_from_2nd_dose/30>=5,"3",
                                                  ifelse(Doses>=3,"4","Unknown")))))) %>%
    mutate(vax_status=fct_recode(vax_status,"unvax"="0","vax_1d"="1","vax_2d_5m"="2","vax_2d_5m_gt"="3","vax_3d"="4","Unknown"="Unknown")) %>%
    mutate(vax_status=fct_relevel(vax_status,"unvax",after=0)) %>% rowwise() %>%
    ##dplyr::filter(vax_1_mfr==vax_2_mfr) %>%
    mutate(vax_status_=ifelse(Doses==0,"unvax",
                              ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="1","jj_1d_5m",
                                     ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="1","jj_1d_5m_gt",
                                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="2","pf_1d_5m",
                                                   ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="2","pf_1d_5m_gt",
                                                          ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="3","mod_1d_5m",
                                                                 ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="3","mod_1d_5m_gt",
                                                                        ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="1","jj_2d_5m",
                                                                               ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="1","jj_2d_5m_gt",
                                                                                      ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m",
                                                                                             ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m_gt",
                                                                                                    ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m",
                                                                                                           ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m_gt",
                                                                                                                  ifelse(Doses>=3,"vax_3d","unknown") )))))))))))))) %>% ungroup() %>%
    mutate(Doses=ifelse(Doses==3 & as.integer(measurement_date-booster_1_date)<=14,Doses-1,
                        ifelse(Doses==2 & as.integer(measurement_date-vax_2_date)<=14,Doses-1,
                               ifelse(Doses==1 & as.integer(measurement_date-vax_1_date)<=14,Doses-1,Doses)  ) )) %>%
    #dplyr::filter(sex!="Unknown" & vax_status!="Unknown" & vax_status_!="Unknown")  %>% 
    mutate(sex=as.factor(sex)) %>% ungroup() %>%
    mutate(ORF1ab_Ct=as.double(gsub(">","",ORF1ab_Ct))) %>%
    group_by(vax_status_) %>% ######dplyr::filter(n()>=50 ) %>%
    mutate(vax_status_=fct_relevel(vax_status_,"unvax",after=0),
           sex=fct_relevel(sex,"Female",after=0),
           age=as.double(age)) %>%
    group_by(hash_subject_id) %>% arrange(measurement_date) %>% 
    dplyr::filter(row_number()==1) %>% ungroup() %>% 
    ####group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
    mutate(result=ifelse(ORF1ab_Ct>30,"negative",test_result)) %>% 
    dplyr::mutate(result=ifelse(result=="1","Omicron",ifelse(result=="0","Delta",result))) #%>% 
    #dplyr::filter(vax_status_!="Unknown")
  
  ####Om.data.tmp$vax_status <- factor(Om.data.tmp$vax_status, levels = c("unvax", "vax_1d", "vax_2d_5m", "vax_2d_5m_gt", "vax_3d"))
  ####Om.data.tmp$vax_status_ <- factor(Om.data.tmp$vax_status_, levels = c("unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d"))
  
  {
    {
      case.counts<-table(Om.data.tmp$result,Om.data.tmp$vax_status)
      case.counts<-cbind(case.counts,"ALL"=rowSums(case.counts))[,c("ALL","unvax","vax_1d","vax_2d_5m","vax_2d_5m_gt","vax_3d")]
      case.totals<-rbind(colSums(case.counts),colSums(case.counts),colSums(case.counts),colSums(case.counts))
      case.totals<-case.totals[,c("ALL","unvax","vax_1d","vax_2d_5m","vax_2d_5m_gt","vax_3d")]
      case.totals
      case.counts<-rbind(case.counts,"positive"=colSums(case.counts[-2,]))[,c("ALL","unvax","vax_1d","vax_2d_5m","vax_2d_5m_gt","vax_3d")]
      case.counts

      p1<-case.counts/case.totals; ser<-sqrt(p1*(1-p1)/case.totals)
      p<-round(p1,3)
      p.low<-round(p1-1.96*ser,3)
      p.high<-round(p1+1.96*ser,3)
      tmp.df<-as.data.frame(matrix(gsub("\\:","\\: ",gsub("\\("," \\(",gsub(" ","",paste(case.counts,"(",p,":",p.low,",",p.high,")")))),
                                   nrow=4,ncol=6,byrow=FALSE))
      rownames(tmp.df)<-rownames(case.counts)
      colnames(tmp.df)<-c("ALL","0","1","2d6m_gt","2d6m","3")
      prop.vals<-tmp.df#[!rownames(tmp.df) %in% c("negative"),] 
      prop.vals %>% dplyr::select("ALL","0","1","2d6m_gt","2d6m","3") %>%
        write.table(file="test.pos.rate_XX.tsv",sep="\t",row.names=FALSE)
    }
    
    p<-p1
    p.low<-p1-1.96*ser
    p.high<-p1+1.96*ser
    var.prop<-p[!rownames(p) %in% c("negative","positive"),] 
    var.prop.low<-p.low[!rownames(p.low) %in% c("negative","positive"),] 
    var.prop.high<-p.high[!rownames(p.high) %in% c("negative","positive"),] 

    par(mai=mai.vals )
    ylab.vals<-pretty(c(0,0.06))
    wdr<-barplot(var.prop,las=1,ylim=c(min(ylab.vals),max(ylab.vals)),
                 col=viridis::inferno(5)[c(4,3)],beside = TRUE,yaxt="n",xaxt="n",
                 xlab="Vaccine doses",ylab="Test positivity rate")
    text(colMeans(wdr),0,
         colnames(var.prop),xpd=TRUE,srt=45,adj=1.15)
    axis(2,at=ylab.vals,labels=ylab.vals,las=1,
         xlab="Vaccine doses")
    for(i in 1:dim(var.prop)[2]){
      segments(wdr[,i],var.prop.low[,i],
               wdr[,i],var.prop.high[,i],lwd=1.15)
    }
    legend("topleft",legend=rownames(var.prop.low),fill=viridis::inferno(5)[c(4,3)],
           cex=0.85,bty="n") 
  }
    
    #mtext(side=1,text="Vaccine doses",line=2.5)
    #mtext(side=2,text="Test positivity rate",line=3.0)
    
  Om.data.tmp<-Om.data %>% mutate(result=ifelse(ORF1ab_Ct>30,"negative",result)) %>% #rowwise() %>%
    dplyr::rename(S_Ct=`S-gene_Ct`,test_result=result) %>% 
    mutate(measurement_date=as.Date(measurement_date,"%m/%d/%y"),
           vax_2_date=as.Date(vax_2_date,"%m/%d/%y"),
           vax_1_date=as.Date(vax_1_date,"%m/%d/%y"),
           booster_1_date=as.Date(booster_1_date,"%m/%d/%y"),
           Doses=as.integer(Doses)) %>% #rowwise() %>%
    mutate(dose1=ifelse(Doses==1,1,0),
           age=ifelse(age==">=90","90",age),
           sex=ifelse(sex=="","Unknown",sex),
           days_from_2nd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_2_date,"%m/%d/%y")),
           days_from_1st_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_1_date,"%m/%d/%y")),
           days_from_3rd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(booster_1_date,"%m/%d/%y"))) %>%
    #dplyr::filter(vax_1_mfr!="unknown" | vax_2_mfr!="unknown" | booster_1_mfr!="unknown") %>%
    ##dplyr::filter(test_result!="negative") %>%
    mutate(vax_1_mfr=ifelse(vax_1_mfr=="J&J","1",
                            ifelse(vax_1_mfr=="Pfizer","2",
                                   ifelse(vax_1_mfr=="Moderna","3","")))) %>%
    mutate(vax_2_mfr=ifelse(vax_2_mfr=="J&J","1",
                            ifelse(vax_2_mfr=="Pfizer","2",
                                   ifelse(vax_2_mfr=="Moderna","3","")))) %>%
    mutate(booster_1_mfr=ifelse(booster_1_mfr=="J&J","1",
                                ifelse(booster_1_mfr=="Pfizer","2",
                                       ifelse(booster_1_mfr=="Moderna","3","")))) %>%
    mutate(test_result=ifelse(test_result=="Omicron",1,0)) %>% 
    mutate(vax_status=ifelse(Doses==0,"0",
                             ifelse(Doses==1,"1",
                                    ifelse(Doses==2 & days_from_2nd_dose/30<5,"2",
                                           ifelse(Doses==2 & days_from_2nd_dose/30>=5,"3",
                                                  ifelse(Doses>=3,"4","Unknown")))))) %>%
    mutate(vax_status=fct_recode(vax_status,"unvax"="0","vax_1d"="1","vax_2d_5m"="2","vax_2d_5m_gt"="3","vax_3d"="4","Unknown"="Unknown")) %>%
    mutate(vax_status=fct_relevel(vax_status,"unvax",after=0)) %>% rowwise() %>%
    ##dplyr::filter(vax_1_mfr==vax_2_mfr) %>%
    mutate(vax_status_=ifelse(Doses==0,"unvax",
                              ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="1","jj_1d_5m",
                                     ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="1","jj_1d_5m_gt",
                                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="2","pf_1d_5m",
                                                   ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="2","pf_1d_5m_gt",
                                                          ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="3","mod_1d_5m",
                                                                 ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="3","mod_1d_5m_gt",
                                                                        ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="1","jj_2d_5m",
                                                                               ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="1","jj_2d_5m_gt",
                                                                                      ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m",
                                                                                             ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m_gt",
                                                                                                    ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m",
                                                                                                           ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m_gt",
                                                                                                                  ifelse(Doses>=3,"vax_3d","unknown") )))))))))))))) %>% ungroup() %>%
    mutate(Doses=ifelse(Doses==3 & as.integer(measurement_date-booster_1_date)<=14,Doses-1,
                        ifelse(Doses==2 & as.integer(measurement_date-vax_2_date)<=14,Doses-1,
                               ifelse(Doses==1 & as.integer(measurement_date-vax_1_date)<=14,Doses-1,Doses)  ) )) %>%
    #dplyr::filter(sex!="Unknown" & vax_status!="Unknown" & vax_status_!="Unknown")  %>% 
    mutate(sex=as.factor(sex)) %>% ungroup() %>%
    mutate(ORF1ab_Ct=as.double(gsub(">","",ORF1ab_Ct))) %>%
    group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
    mutate(vax_status_=fct_relevel(vax_status_,"unvax",after=0),
           sex=fct_relevel(sex,"Female",after=0),
           age=as.double(age)) %>%
    group_by(hash_subject_id) %>% arrange(measurement_date) %>% 
    dplyr::filter(row_number()==1) %>% ungroup() %>% 
    ####group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
    mutate(result=ifelse(ORF1ab_Ct>30,"negative",test_result)) %>% 
    dplyr::mutate(result=ifelse(result=="1","Omicron",ifelse(result=="0","Delta",result))) #%>% 
  #dplyr::filter(vax_status_!="Unknown")
  
    {
      {
        case.counts<-table(Om.data.tmp$result,Om.data.tmp$vax_status_)
        case.counts<-cbind(case.counts,"ALL"=rowSums(case.counts))[,c("ALL","unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d")]
        case.totals<-rbind(colSums(case.counts),colSums(case.counts),colSums(case.counts),colSums(case.counts))
        case.totals<-case.totals[,c("ALL","unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d")]
        case.totals
        case.counts<-rbind(case.counts,"positive"=colSums(case.counts[-2,]))[,c("ALL","unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d")]
        case.counts
        
        p1<-case.counts/case.totals; ser<-sqrt(p1*(1-p1)/case.totals)
        p<-round(p1,3)
        p.low<-round(p1-1.96*ser,3)
        p.high<-round(p1+1.96*ser,3)
        tmp.df<-as.data.frame(matrix(gsub("\\:","\\: ",gsub("\\("," \\(",gsub(" ","",paste(case.counts,"(",p,":",p.low,",",p.high,")")))),
                                     nrow=4,ncol=7,byrow=FALSE))
        rownames(tmp.df)<-rownames(case.counts)
        colnames(tmp.df)<-c("ALL","unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d") ##c("ALL","0","1","2d6m_gt","2d6m","3")
        prop.vals<-tmp.df[!rownames(tmp.df) %in% c("negative"),] 
        prop.vals %>% dplyr::select("ALL","unvax", "jj_1d_5m_gt", "pf_2d_5m", "pf_2d_5m_gt", "mod_2d_5m_gt", "vax_3d") %>%
          write.table(file="test.pos.rate_.tsv",sep="\t",row.names=FALSE)
        }
      
      p<-p1
      p.low<-p1-1.96*ser
      p.high<-p1+1.96*ser
      var.prop<-p[!rownames(p) %in% c("negative","positive"),] 
      var.prop.low<-p.low[!rownames(p.low) %in% c("negative","positive"),] 
      var.prop.high<-p.high[!rownames(p.high) %in% c("negative","positive"),] 
      
      par(mai=mai.vals1 )
      ylab.vals<-pretty(c(0,0.08))
      wdr<-barplot(var.prop,las=1,ylim=c(min(ylab.vals),max(ylab.vals)),
                   col=viridis::inferno(5)[c(4,3)],beside = TRUE,yaxt="n",xaxt="n",
                   xlab="Vaccine doses",ylab="Test positivity rate")
      text(colMeans(wdr),0,
           colnames(var.prop),xpd=TRUE,srt=45,adj=1.15)
      axis(2,at=ylab.vals,labels=ylab.vals,las=1,
           xlab="Vaccine doses")
      for(i in 1:dim(var.prop)[2]){
        segments(wdr[,i],var.prop.low[,i],
                 wdr[,i],var.prop.high[,i],lwd=1.15)
      }
  }
  
  
  {
    par(mai=mai.vals)
    glm.out.var.glm<-glm(test_result~vax_status+age+sex,data=Om.data.corr,family=binomial)
    glm.out.var<-as.data.frame(coef(summary(glm.out.var.glm)))

    glm.out<-glm.out.var %>% rownames_to_column("term") %>% 
      mutate(conf.low=Estimate-1.96*`Std. Error`,
             conf.high=Estimate+1.96*`Std. Error`,
             std.error=`Std. Error`,
             p.value=`Pr(>|z|)`,statistic=`z value`) %>%
      mutate(estimate=exp(Estimate),conf.low=exp(conf.low),conf.high=exp(conf.high)) %>% arrange(term) %>%
      dplyr::select(term,estimate,std.error,statistic,p.value,conf.low,conf.high)
    
    plot.odds.ratio(glm.out,as.is=TRUE)
    title("Positive test: Omicron vs Delta",font.main=1,cex.main=1.0,line=1)
    
    glm.out %>% dplyr::filter(term!="(Intercept)") %>% 
      mutate(across(where(is.numeric), round, digits=4)) %>%
      mutate(estimate=paste0(estimate," (",conf.low,",",conf.high,")")) %>% #arrange(term) %>%
      dplyr::select(term,estimate,p.value) %>% write.table(file="Om_delta.reg__.tsv",sep="\t",row.names=FALSE)
    
  }
  
  {
    par(mai=mai.vals1)
    glm.out.var.glm<-glm(test_result~vax_status_+age+sex,data=Om.data.corr,family=binomial)
    glm.out.var<-as.data.frame(coef(summary(glm.out.var.glm)))

    glm.out<-glm.out.var %>% rownames_to_column("term") %>% 
      mutate(conf.low=Estimate-1.96*`Std. Error`,
             conf.high=Estimate+1.96*`Std. Error`,
             std.error=`Std. Error`,
             p.value=`Pr(>|z|)`,statistic=`z value`) %>%
      mutate(estimate=exp(Estimate),conf.low=exp(conf.low),conf.high=exp(conf.high)) %>% #arrange(term) %>%
      dplyr::select(term,estimate,std.error,statistic,p.value,conf.low,conf.high)
    
    plot.odds.ratio(glm.out,as.is=TRUE)
    title("Positive test: Omicron vs Delta",font.main=1,cex.main=1.0,line=1)
    
    glm.out %>% dplyr::filter(term!="(Intercept)") %>% 
      mutate(across(where(is.numeric), round, digits=4)) %>%
      mutate(estimate=paste0(estimate," (",conf.low,",",conf.high,")")) %>% #arrange(term) %>%
      dplyr::select(term,estimate,p.value) %>% write.table(file="Om_delta.reg.tsv",sep="\t",row.names=FALSE)
  }
  
  dev.off()
}



{
  pdf(file="Figure1_draft_v1.pdf",height=8.0/2.0,width=9.05)
  layout(matrix(1:2,ncol=2),widths=c(1,0.115))
  mai.vals1<-c(0.961,0.851,0.561,0.050)
  par(mai=mai.vals1 )
  
  #tmpW<-dat_cases_long$Date[dat_cases_long$Counter %in% xlim.vals] %>% as_tibble() %>% 
  #  mutate(XX=paste0(lubridate::month(value),"/",lubridate::day(value),"/",gsub("2022","21",gsub("2021","21",lubridate::year(value))))) %>% 
  #  dplyr::select(XX) %>% as.data.frame()
  
  tmpW<-colnames(t(variant.cases1)) %>% as_tibble() %>% 
    mutate(XX=paste0(lubridate::month(value),"/",lubridate::day(value),"/",gsub("2022","21",gsub("2021","21",lubridate::year(value))))) %>% 
    dplyr::select(XX) %>% as.data.frame()
  
  x <- barplot(t(variant.cases1),beside=FALSE,las=1,col=c(viridis::inferno(5)[4],viridis::inferno(5)[3]),
               ylim=c(0,1),ylab="Proportion",xaxs="i",
               names.arg=tmpW$XX,
               las=2,cex.names=0.5,
               xlab="Test date")
  
  mai.vals1<-c(0.961,0.051,0.561,0.050)
  par(mai=mai.vals1 )
  
  plot.new()
  legend("center",legend=c("Delta","Omicron"),xpd=TRUE,
         fill=c(viridis::inferno(5)[4],viridis::inferno(5)[3]),
         bty="n",cex=0.80) #,title="Variant"
  
  dev.off()
}

{
  pdf(file="PCR_figure.pdf",height=5.0,width=6.5)
  mai.vals<-c(0.861,0.851,0.161,0.150)
  par(mai=mai.vals,xpd=TRUE )
  
  {
    par(bty="n",mai=mai.vals)
    tmp.Om<-Om.data %>% mutate(result=ifelse(ORF1ab_Ct>30,"negative",result)) %>% #rowwise() %>%
      dplyr::rename(S_Ct=`S-gene_Ct`,test_result=result) %>% 
      mutate(measurement_date=as.Date(measurement_date,"%m/%d/%y"),
             vax_2_date=as.Date(vax_2_date,"%m/%d/%y"),
             vax_1_date=as.Date(vax_1_date,"%m/%d/%y"),
             booster_1_date=as.Date(booster_1_date,"%m/%d/%y"),
             Doses=as.integer(Doses)) %>% #rowwise() %>%
      mutate(dose1=ifelse(Doses==1,1,0),
             age=ifelse(age==">=90","90",age),
             sex=ifelse(sex=="","Unknown",sex),
             days_from_2nd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_2_date,"%m/%d/%y")),
             days_from_1st_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(vax_1_date,"%m/%d/%y")),
             days_from_3rd_dose=as.integer(as.Date(measurement_date,"%m/%d/%y")-as.Date(booster_1_date,"%m/%d/%y"))) %>%
      #dplyr::filter(vax_1_mfr!="unknown" | vax_2_mfr!="unknown" | booster_1_mfr!="unknown") %>%
      ##dplyr::filter(test_result!="negative") %>%
      mutate(vax_1_mfr=ifelse(vax_1_mfr=="J&J","1",
                              ifelse(vax_1_mfr=="Pfizer","2",
                                     ifelse(vax_1_mfr=="Moderna","3","")))) %>%
      mutate(vax_2_mfr=ifelse(vax_2_mfr=="J&J","1",
                              ifelse(vax_2_mfr=="Pfizer","2",
                                     ifelse(vax_2_mfr=="Moderna","3","")))) %>%
      mutate(booster_1_mfr=ifelse(booster_1_mfr=="J&J","1",
                                  ifelse(booster_1_mfr=="Pfizer","2",
                                         ifelse(booster_1_mfr=="Moderna","3","")))) %>%
      mutate(test_result=ifelse(test_result=="Omicron",1,0)) %>% 
      mutate(vax_status=ifelse(Doses==0,"0",
                               ifelse(Doses==1,"1",
                                      ifelse(Doses==2 & days_from_2nd_dose/30<5,"2",
                                             ifelse(Doses==2 & days_from_2nd_dose/30>=5,"3",
                                                    ifelse(Doses>=3,"4","Unknown")))))) %>%
      mutate(vax_status=fct_recode(vax_status,"unvax"="0","vax_1d"="1","vax_2d_5m"="2","vax_2d_5m_gt"="3","vax_3d"="4","Unknown"="Unknown")) %>%
      mutate(vax_status=fct_relevel(vax_status,"unvax",after=0)) %>% rowwise() %>%
      ##dplyr::filter(vax_1_mfr==vax_2_mfr) %>%
      mutate(vax_status_=ifelse(Doses==0,"unvax",
                                ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="1","jj_1d_5m",
                                       ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="1","jj_1d_5m_gt",
                                              ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="2","pf_1d_5m",
                                                     ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="2","pf_1d_5m_gt",
                                                            ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30<5 & vax_1_mfr=="3","mod_1d_5m",
                                                                   ifelse(Doses==1 & as.numeric(days_from_1st_dose)/30>=5 & vax_1_mfr=="3","mod_1d_5m_gt",
                                                                          ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="1","jj_2d_5m",
                                                                                 ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="1","jj_2d_5m_gt",
                                                                                        ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m",
                                                                                               ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="2" & vax_2_mfr=="2","pf_2d_5m_gt",
                                                                                                      ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30<5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m",
                                                                                                             ifelse(Doses==2 & as.numeric(days_from_2nd_dose)/30>=5 & vax_1_mfr=="3" & vax_2_mfr=="3","mod_2d_5m_gt",
                                                                                                                    ifelse(Doses>=3,"vax_3d","unknown") )))))))))))))) %>% ungroup() %>%
      mutate(Doses=ifelse(Doses==3 & as.integer(measurement_date-booster_1_date)<=14,Doses-1,
                          ifelse(Doses==2 & as.integer(measurement_date-vax_2_date)<=14,Doses-1,
                                 ifelse(Doses==1 & as.integer(measurement_date-vax_1_date)<=14,Doses-1,Doses)  ) )) %>%
      dplyr::filter(sex!="Unknown" & vax_status!="Unknown" & vax_status_!="Unknown")  %>% 
      mutate(sex=as.factor(sex)) %>% ungroup() %>%
      mutate(ORF1ab_Ct=as.double(gsub(">","",ORF1ab_Ct))) %>%
      group_by(vax_status_) %>% dplyr::filter(n()>=50 ) %>%
      mutate(vax_status_=fct_relevel(vax_status_,"unvax",after=0),
             sex=fct_relevel(sex,"Female",after=0),
             age=as.double(age)) %>%
      group_by(hash_subject_id) %>% arrange(measurement_date) %>% 
      dplyr::filter(row_number()==1) %>% ungroup() %>% 
      mutate(result=ifelse(ORF1ab_Ct>30,"negative",test_result)) %>% 
      dplyr::mutate(result=ifelse(result=="1","Omicron",ifelse(result=="0","Delta",result))) %>% 
      dplyr::filter(vax_status_!="Unknown")
    
    tmp.data<-list("ALL"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta"],
                   "0"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==0],
                   "1"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==1],
                   "2"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==2],
                   "3"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Delta" & tmp.Om$Doses==3])

    tmp.data1<-list("ALL"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron"],
                    "0"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==0],
                   "1"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==1],
                   "2"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==2],
                   "3"=tmp.Om$ORF1ab_Ct[tmp.Om$result=="Omicron" & tmp.Om$Doses==3]) 
    y.vals=pretty(c(5,35))
    vioplot::vioplot(tmp.data,xaxs="i",yaxs="i",ylim=c(min(y.vals),max(y.vals)),
                     las=1, plotCentre = "line", side = "left",lwd=0.001,
                     col=viridis::inferno(5)[c(4)], border="white", rectCol="palevioletred", lineCol="violetred", colMed="violet",
                     xlab="Vaccination status",ylab=bquote("Nasal swab PCR C"[t]~"-value"),yaxt="n",xaxt="n")
    axis(2,at=y.vals,labels=y.vals,las=1)
    axis(1,at=1:length(names(tmp.data)),labels=names(tmp.data),las=1)
    
    vioplot::vioplot(tmp.data1,
                     las=1, plotCentre = "line", side = "right",lwd=0.001,add = TRUE,
                     col=viridis::inferno(5)[c(3)], border="white", rectCol="palevioletred", lineCol="violetred", colMed="violet",
                     #xlab="Vaccination status",ylab=bquote("Nasal swab PCR C"[t]~"-value"),
                     yaxt="n",xaxt="n")
    
    legend("topleft",legend=c("Delta","Omicron"),#title="Doses",
           cex=0.75,
           fill=c(viridis::inferno(5)[c(4,3)]),bty="n")
    
  }
  dev.off()
  
  ##Testing for differences in the PCR Ct values for individuals infected with Omicron or Delta variants

  kruskal.test(list(tmp.data$`ALL`,tmp.data1$`ALL`))
  kruskal.test(list(tmp.data$`0`,tmp.data1$`0`))
  kruskal.test(list(tmp.data$`1`,tmp.data1$`1`))
  kruskal.test(list(tmp.data$`2`,tmp.data1$`2`))
  kruskal.test(list(tmp.data$`3`,tmp.data1$`3`))
  
  ##Mean and confidence intervals for the PCR Ct values for individuals infected with Omicron and Delta variants

  c(mean(tmp.data$ALL), mean(tmp.data$ALL)-1.96*sd((tmp.data$ALL)), mean(tmp.data$ALL)+1.96*sd((tmp.data$ALL)))
  c(mean(tmp.data1$ALL), mean(tmp.data1$ALL)-1.96*sd((tmp.data1$ALL)), mean(tmp.data1$ALL)+1.96*sd((tmp.data1$ALL)))
  
  c(mean(tmp.data$`0`), mean(tmp.data$`0`)-1.96*sd((tmp.data$`0`)), mean(tmp.data$`0`)+1.96*sd((tmp.data$`0`)))
  c(mean(tmp.data1$`0`), mean(tmp.data1$`0`)-1.96*sd((tmp.data1$`0`)), mean(tmp.data1$`0`)+1.96*sd((tmp.data1$`0`)))
  
  c(mean(tmp.data$`1`), mean(tmp.data$`1`)-1.96*sd((tmp.data$`1`)), mean(tmp.data$`1`)+1.96*sd((tmp.data$`1`)))
  c(mean(tmp.data1$`1`), mean(tmp.data1$`1`)-1.96*sd((tmp.data1$`1`)), mean(tmp.data1$`1`)+1.96*sd((tmp.data1$`1`)))
  
  c(mean(tmp.data$`2`), mean(tmp.data$`2`)-1.96*sd((tmp.data$`2`)), mean(tmp.data$`2`)+1.96*sd((tmp.data$`2`)))
  c(mean(tmp.data1$`2`), mean(tmp.data1$`2`)-1.96*sd((tmp.data1$`2`)), mean(tmp.data1$`2`)+1.96*sd((tmp.data1$`2`)))
  
  c(mean(tmp.data$`3`), mean(tmp.data$`3`)-1.96*sd((tmp.data$`3`)), mean(tmp.data$`3`)+1.96*sd((tmp.data$`3`)))
  c(mean(tmp.data1$`3`), mean(tmp.data1$`3`)-1.96*sd((tmp.data1$`3`)), mean(tmp.data1$`3`)+1.96*sd((tmp.data1$`3`)))
}
