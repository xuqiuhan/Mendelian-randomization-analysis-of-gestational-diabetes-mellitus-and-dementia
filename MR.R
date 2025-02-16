####多变量MR####
setwd("C:/Program Files (x86)/R")
library(TwoSampleMR)
library(data.table)
library(dplyr)
a<-fread(input="TWB.BFR.tsv",sep="\t")
a <-  a %>%dplyr::rename(SNP=rs_id,"beta.exposure"=beta,
                                                     "se.exposure"=standard_error,
                                                     "effect_allele.exposure"=effect_allele,
                                                     "other_allele.exposure"=other_allele,
                                                     "pval.exposure"=p_value,
                                                     "eaf.exposure"=effect_allele_frequency)%>%
  dplyr::select(SNP,"effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","eaf.exposure","pval.exposure")
a<-data.frame(a)
exp_a <- a[a$pval.exposure < 5e-7,]
#clump数据，去除连锁不平衡，人种设立为EAS
library("ieugwasr")
library(plinkbinr)
expo_clump <- ieugwasr::ld_clump(dplyr::tibble(rsid=exp_a$SNP, pval=exp_a$pval.exposure),
                                 clump_r2 = 0.001,
                                 clump_p = 1,
                                 clump_kb = 10000,
                                 plink_bin = plinkbinr::get_plink_exe(),
                                 bfile = "C:/Program Files (x86)/R/1kg.v3/1kg.v3/EAS",pop = "EAS")
exp_a <- exp_a%>%
  dplyr::filter(SNP%in%expo_clump$rsid)
exp_b <- TwoSampleMR::extract_instruments(outcomes = "ebi-a-GCST004904",
                                                          p1 = 5e-07,
                                                          clump = T,
                                                          r2 = 0.001,
                                                         kb = 10000)
Sys.setenv(OPENGWAS_JWT="eyJhbGciOiJSUzI1NiIsImtpZCI6ImFwaS1qd3QiLCJ0eXAiOiJKV1QifQ.eyJpc3MiOiJhcGkub3Blbmd3YXMuaW8iLCJhdWQiOiJhcGkub3Blbmd3YXMuaW8iLCJzdWIiOiJjaGVueXc3NTFAMTYzLmNvbSIsImlhdCI6MTczODg0OTAyNCwiZXhwIjoxNzQwMDU4NjI0fQ.OlixXGcBDkPHEVtoquYv-LO-wAQcQr-D7dlttjogrhD6CjpTs4yzu_M2KspWe1DcyJ2jdNVCyvgbxukLFIxdnsJkUHARBTfQ0EKbeng7KUlBDRlWrls2nnlqqV-Z_7dQmiEY_YUuqjrshXCF0UFVPJA_XJk8Rsq8D11s5uv7MVHlJiu-Fm6OGz2awszdgVVxLxjq_GZ8AclASeRPKKWMEnGpGhFJ-gd732Nc4uwixEbJ0QC_hQeKEXrHFPNWn3dhDH3E822mKvZ4c6OqudFqHYkbliEMojlSjux_05wbUuQxxH2quPt3PRIkYSGachoztDClTFLrn56Jtm3Qr37Cjg")
exp_a <- as.data.table(exp_a)
exp_b <- as.data.table(exp_b)
exp_a<-exp_a[,"SNP"]
exp_b<-exp_b[,"SNP"]
exp<-rbind(exp_a,exp_b)
exp<-unique(exp) 
exp_a<-merge(exp,a,by.x = "SNP",by.y = "SNP")
exp_b<- extract_outcome_data(snps = exp$SNP,outcomes = "ebi-a-GCST004904",proxies = F)
exp_c<-extract_outcome_data(snps = exp$SNP,outcomes = "ebi-a-GCST90018696",proxies = F)
exp_a <- exp_a %>%
  distinct(SNP,.keep_all = T)
exp_b <- exp_b %>%
  distinct(SNP,.keep_all = T)
exp_c <- exp_c %>%
  distinct(SNP,.keep_all = T)
merge_a1 <- intersect(exp_a$SNP,exp_b$SNP )
merge_a2 <- intersect(merge_a1,exp_c$SNP)

exp_a <- exp_a[exp_a$SNP %in% merge_a2]
exp_b <- exp_b[exp_b$SNP %in% merge_a2,]
exp_c <- exp_c[exp_c$SNP %in% merge_a2,]
exp_a$beta.exposure <- as.numeric(exp_a$beta.exposure)
exp_a$se.exposure <- as.numeric(exp_a$se.exposure)
exp_b$beta.outcome <- as.numeric(exp_b$beta.outcome)
exp_b$se.outcome <- as.numeric(exp_b$se.outcome)
library(MendelianRandomization)
MRMVInputObject <- mr_mvinput(bx = cbind(exp_a$beta.exposure,exp_b$beta.outcome),
                              bxse = cbind(exp_a$se.exposure,exp_b$se.outcome),
                              by = exp_c$beta.outcome,                                 byse = exp_c$se.outcome)
MRMVInputObject
MRMVObject_IVW <- mr_mvivw(MRMVInputObject, 
                           model = "default",
                           correl = FALSE,
                           distribution = "normal",
                           alpha = 0.05)

MRMVObject_IVW 





# 示例系数 (Beta) 值
beta_values <- c(0.744,-0.251)

# 示例标准误 (SE) 值
se_values <- c(0.299,0.233)
#将p值输入进去
p_values <- c("0.013","0.282")
# 创建一个数据框来存储结果
results <- data.frame(
  Beta = beta_values,
  SE = se_values,
  OR = exp(beta_values),
  LCI = beta_values - 1.96 * se_values,
  UCI = beta_values + 1.96 * se_values,
  Lower_CI = exp(beta_values - 1.96 * se_values),
  Upper_CI = exp(beta_values + 1.96 * se_values),
  p = p_values
)
results <- results %>% 
  dplyr::mutate("Exposure" = c("Body fat rate","Body mass index"),"Outcome" = "Obstructive sleep apena")
write.csv(results ,file = "D:/R results/mvmr_ivw.csv",row.names = F)
