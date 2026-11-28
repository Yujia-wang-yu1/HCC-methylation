args <- commandArgs(trailingOnly=TRUE)
D <- if(length(args)>=1) args[[1]] else 'data'
O <- if(length(args)>=2) args[[2]] else 'results'
dir.create(O,showWarnings=FALSE,recursive=TRUE)
library(MASS)
if (requireNamespace('EpiDISH', quietly=TRUE) && packageVersion('EpiDISH') == '2.28.0') {
  DoRPC <- getFromNamespace('DoRPC', 'EpiDISH')
} else {
  source(file.path(D, 'EpiDISH_2.28.0_epidish.R'))
}
for(cohort in 'TCGA') {
 f <- list()
 for(nm in c('centEpiFibIC.m','centBloodSub.m')) {
  original_n <- if(nm=='centEpiFibIC.m') 716 else 188
  x <- as.matrix(read.delim(file.path(D,paste0(cohort,'_',nm,'_QC_noMISP.tsv')),row.names=1,check.names=FALSE))
  ref <- as.matrix(read.delim(file.path(D,paste0(cohort,'_',nm,'_reference.tsv')),row.names=1,check.names=FALSE))
  est <- matrix(NA_real_,ncol(x),ncol(ref),dimnames=list(colnames(x),colnames(ref)))
  for(j in seq_len(ncol(x))) {
   keep <- is.finite(x[,j]) & apply(ref,1,function(v)all(is.finite(v)))
   if(sum(keep)>=ceiling(2*original_n/3)) est[j,] <- DoRPC(x[keep,j,drop=FALSE],ref[keep,,drop=FALSE],maxit=50)$estF[1,]
  }
  f[[nm]] <- est
 }
 broad <- f[['centEpiFibIC.m']];blood <- f[['centBloodSub.m']]
 expected <- read.delim(file.path(D,paste0(cohort,'_DNA_composition_proxies.tsv')),check.names=FALSE)
 stopifnot(identical(rownames(broad),expected$sample_id))
 err <- max(abs(broad-as.matrix(expected[,c('Epi','Fib','IC')])),na.rm=TRUE)
 stopifnot(err<1e-9)
 write.table(data.frame(sample_id=rownames(broad),broad,check.names=FALSE),file.path(O,paste0(cohort,'_reproduced_DNA_proxies.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
 cat(cohort,'maximum broad proxy error',err,'eligible',sum(complete.cases(broad)),'\n')
}
