setwd("/your_dir/")

## please make sure the following R library works well
library(survival)
library(ipred)
library(pec)

## "Phenotype_file.csv" each row is a sampleID, each column is a phenotype
Pheno=read.csv("Phenotype_file.csv") ## Pheno=read.table("",sep="\t",header=T)

dim(Pheno)##sample_ID, phenotypes
## please include the following columns in "Phenotype_file.csv":
## "days_to_event_num" -- Days_to_exam_death; "death" -- All_cause_mortality; "CVDDEATH_r" -- CVD_death; "Cancer_death";
## "age"; "sex"; "bmi"; "smoke"; "gramsperday" -- Alcohol Drinking; "EduYears"; "g2ex8phy" -- hysical Activity;
## "diabetes"-- prevalent diabetes; "prevCHD"-- prevalent CHD; "prevCHF" -- prevalent heart failure; "prevABI" -- prevalent stroke;
## "hrtnow" -- prevalent hypertension; "cancer" -- "prevalent cancer"

## Do not have NA in your data!
Pheno=Pheno[,c("sample_ID", "days_to_event_num","death", "CVDDEATH_r", "Cancer_death",
               "age", "sex", "gramsperday","smoke","bmi","diabetes",
               "prevCHD","hrtnow","cancer","g2ex8phy", "EduYears","prevCHF","prevABI")]
		 
	
## Read your DNA methylation data
## DNAm_Data=load("your_dna_methy_data.Rdata")
## beta values of DNAm, each row is a sampleID, each col is a CpG
load("DNAm_Data.Rdata")

ls()
dim(DNAm_Data) ##sample_ID,CpG_number
#DNAm_Data[1:3,1:3]

## We only need a subset of CpGs
CpGs_p05=read.csv("CpGs_EWAS.csv")

DNAm_Data_sub=DNAm_Data[,c("sample_ID",paste(CpGs_p05$CpG))]
dim(DNAm_Data_sub)

## Inverly-transformed the DNAm Beta values
DNAm_Data_s_inv=DNAm_Data_sub
row.names(DNAm_Data_s_inv)=DNAm_Data_sub[,1]
DNAm_Data_s_inv=DNAm_Data_s_inv[,-1]

for(i in 1:dim(DNAm_Data_s_inv)[2]) {
DNAm_Data_s_inv[,i]=qnorm(rank(DNAm_Data_s_inv[,i])/(length(DNAm_Data_s_inv[,i])+1),mean=0,sd=1)
}

## merge Pheno and DNAm_Data
Pheno_DNAm=merge(Pheno,DNAm_Data_s_inv, by.x="sample_ID", by.y="row.names")

## put all FHS models in a folder
## unzip the file:  tar -zxvf fhs_mods.tar.gz
FHS_mods_dir="/yourpath/fhs_mods/"
FHS_mods_names=dir(FHS_mods_dir)
length(FHS_mods_names)## 20

Pheno_riskScore=Pheno_DNAm[,c("sample_ID", "days_to_event_num","All_cause_mortality", "CVD_death", "Cancer_death",
               "age", "sex", "gramsperday","smoke","bmi","diabetes",
               "prevCHD","hrtnow","cancer","g2ex8phy", "EduYears","prevCHF","prevABI")]
	   
## read each model, calculate the riskScore, and merge riskScore with your phenotypes
for (i in 1:length(FHS_mods_names)){

c_mod=FHS_mods_names[i]
c_mod_dir=paste(FHS_mods_dir,c_mod, sep="")
FHS_mod = read.csv(c_mod_dir)

formula_riskScore=paste(FHS_mod[1:dim(FHS_mod)[1],1],"*Pheno_DNAm[,",'\"',FHS_mod[1:dim(FHS_mod)[1],2],'\"',"]",sep="")
formula_riskScore=toString(formula_riskScore,quote=F)
formula_riskScore=gsub('],','] + ', formula_riskScore)

riskScore = eval(parse(text=paste(formula_riskScore)))

DNAm_Mortality_RiskScore=as.matrix(riskScore)
DNAm_Mortality_RiskScore=as.numeric(qnorm(rank(DNAm_Mortality_RiskScore)/(length(DNAm_Mortality_RiskScore)+1),mean=0,sd=1))   

Pheno_riskScore[,paste(c_mod)]=DNAm_Mortality_RiskScore

}



#################################################
## Evaluate the DNAm_Mortality_RiskScore for predict all-cause mortality / CVD death / Cancer death

results_allcause=array(NA,c(0,13))
colnames(results_allcause)=c("Model", "Beta", "HR", "se", "z", "p", "CI.l", "CI.h", "C-index", "C-index.se", 
                    "chisq", "log-rank P", "Integrated_Brier_Score")
					
results_CVD=array(NA,c(0,13))
colnames(results_CVD)=c("Model", "Beta", "HR", "se", "z", "p", "CI.l", "CI.h", "C-index", "C-index.se", 
                    "chisq", "log-rank P", "Integrated_Brier_Score")

results_cancer=array(NA,c(0,13))
colnames(results_cancer)=c("Model", "Beta", "HR", "se", "z", "p", "CI.l", "CI.h", "C-index", "C-index.se", 
                    "chisq", "log-rank P", "Integrated_Brier_Score")

Pheno_riskScore$CVDDEATH_r[is.na(Pheno_riskScore$CVDDEATH_r)]=0
					
for(i in 1:20){

c_mod=FHS_mods_names[i]


Pheno_riskScore$risk.inv=Pheno_riskScore[,paste(c_mod)]
Pheno_riskScore$risk.inv_rank=0
median(Pheno_riskScore[,paste(c_mod)])
Pheno_riskScore$risk.inv_rank[Pheno_riskScore$risk.inv>=median(Pheno_riskScore[,paste(c_mod)])]=1

## all-cause
# coxph
res_cox<-coxph(Surv(days_to_event_num,death)~ risk.inv, data=Pheno_riskScore,na.action=na.omit,x = TRUE)
sum_cox=summary(res_cox)    
coef_cox=sum_cox$coefficients[1,]
conf=sum_cox$conf.int[1,c(3,4)]
C_index=sum_cox$"concordance"[1]
C_se=sum_cox$"concordance"[2]
# logRank
fit <- survfit(coxph(Surv(days_to_event_num,death)~ strata(risk.inv_rank), data=Pheno_riskScore,na.action=na.omit))
sdf=survdiff(Surv(days_to_event_num,death)~ risk.inv_rank, data=Pheno_riskScore, na.action=na.omit)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# Brier Score
smod <- Surv(Pheno_riskScore$days_to_event_num, Pheno_riskScore$death)
sbrier(smod,Pheno_riskScore$risk.inv,btime= max(smod[,1]))
	
perror=pec(list(CPH=res_cox),formula=Surv(days_to_event_num,death)~risk.inv,data=Pheno_riskScore)
crps(perror)[2]

c_res=c(c_mod,coef_cox,conf,C_index,C_se,sdf$chisq,p.val,crps(perror)[2])
results_allcause=rbind(results_allcause,c_res)


## CVD-death
# coxph
res_cox<-coxph(Surv(days_to_event_num,CVDDEATH_r)~ risk.inv, data=Pheno_riskScore,na.action=na.omit,x = TRUE)
sum_cox=summary(res_cox)    
coef_cox=sum_cox$coefficients[1,]
conf=sum_cox$conf.int[1,c(3,4)]
C_index=sum_cox$"concordance"[1]
C_se=sum_cox$"concordance"[2]
# logRank
fit <- survfit(coxph(Surv(days_to_event_num,CVDDEATH_r)~ strata(risk.inv_rank), data=Pheno_riskScore,na.action=na.omit))
sdf=survdiff(Surv(days_to_event_num,CVDDEATH_r)~ risk.inv_rank, data=Pheno_riskScore, na.action=na.omit)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# Brier Score
smod <- Surv(Pheno_riskScore$days_to_event_num, Pheno_riskScore$CVDDEATH_r)
sbrier(smod,Pheno_riskScore$risk.inv,btime= max(smod[,1]))
	
perror=pec(list(CPH=res_cox),formula=Surv(days_to_event_num,CVDDEATH_r)~risk.inv,data=Pheno_riskScore)
crps(perror)[2]

c_res=c(c_mod,coef_cox,conf,C_index,C_se,sdf$chisq,p.val,crps(perror)[2])
results_CVD=rbind(results_CVD,c_res)

## cancer-death
# coxph
res_cox<-coxph(Surv(days_to_event_num,Cancer_death)~ risk.inv, data=Pheno_riskScore,na.action=na.omit,x = TRUE)
sum_cox=summary(res_cox)    
coef_cox=sum_cox$coefficients[1,]
conf=sum_cox$conf.int[1,c(3,4)]
C_index=sum_cox$"concordance"[1]
C_se=sum_cox$"concordance"[2]
# logRank
fit <- survfit(coxph(Surv(days_to_event_num,Cancer_death)~ strata(risk.inv_rank), data=Pheno_riskScore,na.action=na.omit))
sdf=survdiff(Surv(days_to_event_num,Cancer_death)~ risk.inv_rank, data=Pheno_riskScore, na.action=na.omit)
p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
# Brier Score
smod <- Surv(Pheno_riskScore$days_to_event_num, Pheno_riskScore$Cancer_death)
sbrier(smod,Pheno_riskScore$risk.inv,btime= max(smod[,1]))
	
perror=pec(list(CPH=res_cox),formula=Surv(days_to_event_num,Cancer_death)~risk.inv,data=Pheno_riskScore)
crps(perror)[2]

c_res=c(c_mod,coef_cox,conf,C_index,C_se,sdf$chisq,p.val,crps(perror)[2])
results_cancer=rbind(results_cancer,c_res)

}


## save results, and please send me those results

write.table(results_allcause, "aric_EA_elasticNet_allcause.csv",row.names=F,quote=F, sep=",")
write.table(results_CVD, "aric_EA_elasticNet_cvdDeath.csv",row.names=F,quote=F, sep=",")
write.table(results_cancer, "aric_EA_elasticNet_cancerDeath.csv",row.names=F,quote=F, sep=",")






