# Portable recorded-OS sensitivity. Existing R survival package; no installation/network.
args <- commandArgs(trailingOnly=TRUE)
input_dir <- if (length(args)>=1) args[[1]] else 'data'
output_dir <- if (length(args)>=2) args[[2]] else 'results'
dir.create(output_dir,showWarnings=FALSE,recursive=TRUE)
library(survival)
d <- read.delim(file.path(input_dir,'ICGC_survival_inputs.tsv'),check.names=FALSE)
rows <- list(); tests <- list()
for (cohort in c('HCC_203','mixed_231_sensitivity')) {
 x <- if(cohort=='HCC_203') d[d$tumour_histological_type=='8170/3',] else d
 x$score_z <- as.numeric(scale(x$MISP_RNA_suppression_score,center=TRUE,scale=TRUE))
 x$event <- as.integer(x$donor_vital_status=='deceased')
 x$male <- as.integer(x$donor_sex=='male')
 for (model in c('univariable','age_sex_stage')) {
  formula <- if(model=='univariable') Surv(donor_survival_time,event)~score_z else Surv(donor_survival_time,event)~score_z+donor_age_at_diagnosis+male+tumour_stage
  fit <- coxph(formula,data=x,ties='efron',x=TRUE,model=TRUE)
  s <- summary(fit);ph <- cox.zph(fit,transform='km')
  for (term in rownames(s$coefficients)) {
   k <- match(term,rownames(s$coefficients));ci <- s$conf.int[k,]
   rows[[length(rows)+1]] <- data.frame(cohort=cohort,model=model,term=term,n=fit$n,events=fit$nevent,beta=s$coefficients[k,'coef'],se=s$coefficients[k,'se(coef)'],HR=ci['exp(coef)'],lower=ci['lower .95'],upper=ci['upper .95'],p=s$coefficients[k,'Pr(>|z|)'],PH_term_p=ph$table[term,'p'],PH_global_p=ph$table['GLOBAL','p'])
  }
 }
}
out <- do.call(rbind,rows)
# Five coefficient tests within each cohort, matching the stated exploratory family.
out$FDR_within_cohort <- ave(out$p,out$cohort,FUN=function(p) p.adjust(p,method='BH'))
out$FDR_four_MISP_models <- NA_real_
score_rows <- out$term=='score_z'
out$FDR_four_MISP_models[score_rows] <- p.adjust(out$p[score_rows],method='BH')
write.table(out,file.path(output_dir,'ICGC_HCC_survival_sensitivity.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
print(out[out$term=='score_z',]);cat('survival package',as.character(packageVersion('survival')),'\n')
