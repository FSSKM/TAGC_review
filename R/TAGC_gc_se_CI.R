

TAGC_gc_se_CI = function(TAGC_bias, covar_path=NULL, pheno_path, pheno_col_num, prs_path, prs_col_num, number_bs_samples=500, bs_sample_size=1500, CI_conf_level=0.95) {

  # load covar (if given), phenotype, and prs
  pheno_input = as.data.frame(data.table::fread(pheno_path))
  prs_input = as.data.frame(data.table::fread(prs_path))
  if (!is.null(covar_path)) {
    covar_input = as.data.frame(data.table::fread(covar_path))

    # take intersection
    intersect_ids = intersect(intersect(covar_input[,1], pheno_input[,1]), prs_input[,1])
    if (length(intersect_ids) == 0) {
      stop(paste0('There is no overlap of subjects in the first column of ', covar_path, ', ', pheno_path, ', and ', prs_path, '.\n Please make sure the first column contains subject IDs.'))
    }
    covar_input = covar_input[covar_input[,1] %in% intersect_ids,]
    pheno_input = pheno_input[pheno_input[,1] %in% intersect_ids,]
    prs_input = prs_input[prs_input[,1] %in% intersect_ids,]
    # sort
    covar_input = covar_input[order(covar_input[,1]),]
    pheno_input = pheno_input[order(pheno_input[,1]),]
    prs_input = prs_input[order(prs_input[,1]),]
  } else {
    # take intersection
    intersect_ids = intersect(pheno_input[,1], prs_input[,1])
    if (length(intersect_ids) == 0) {
      stop(paste0('There is no overlap of subjects in the first column of ', covar_path, ', ', pheno_path, ', and ', prs_path, '.\n Please make sure the first column contains subject IDs.'))
    }
    pheno_input = pheno_input[pheno_input[,1] %in% intersect_ids,]
    prs_input = prs_input[prs_input[,1] %in% intersect_ids,]
    # sort
    pheno_input = pheno_input[order(pheno_input[,1]),]
    prs_input = prs_input[order(prs_input[,1]),]
  }

  summary2 = matrix(NA, 1, 8)
  tagc = NA
  tagc_se = NA
  tagc_bs_se = NA
  if (!all(is.na(pheno_input[,pheno_col_num]))) {
    if (!is.null(covar_path)) {
      l1<-lm(scale(pheno_input[,pheno_col_num])~scale(prs_input[,prs_col_num])+scale(as.matrix(covar_input[,-1])))
      orig_gc = summary(l1)$coef[2,1]
    } else {
      l1<-lm(scale(pheno_input[,pheno_col_num])~scale(prs_input[,prs_col_num]))
      orig_gc = cor(scale(pheno_input[,pheno_col_num]), scale(prs_input[,prs_col_num]))
    }

    summary2[1,1] <- summary(l1)$coef[2,4] # p-value
    summary2[1,2] <- orig_gc # original gc
    summary2[1,3] <- summary(l1)$coef[2,2] # original SE
    summary2[1,4] <- tagc <- orig_gc*TAGC_bias # TAGC gc
    tagc_se <- summary(l1)$coef[2,2]*TAGC_bias # TAGC SE, being used only if bootstrap SE does not work
    summary2[1,8]<-nobs(l1) # number of samples in lm

    ### if TAGC estimate is beyond [-1,+1]
    if (abs(tagc) > 1) {
      warning('The estimated genetic correlation exceeds the bounds [-1,+1]. This could be attributed to two potential reasons. Firstly, please check if the input h2 is very small. If h2 is extremely low or not significantly different from zero, our estimator may be unstable due to the limited genetic variance/influence in the trait. However, if h2 falls within a typical range (e.g., >0.05), a large genetic correlation may indicate a high level of genetic correlation but does not imply that the genetic correlation exceeds 1. It is recommended to consider the uncertainty associated with the point estimate arising from sample variance and to use confidence intervals when interpreting the results.')
    }

    #####
    ###bootstrap to get the SE and CI
    TEMP<-rep(NA, number_bs_samples)
    bs_sample_size_ub = min(1500, length(pheno_input[,pheno_col_num]))
    if (bs_sample_size > bs_sample_size_ub) {
      bs_sample_size = bs_sample_size_ub
    }
    set.seed(1000)
    for(ii in 1:length(TEMP)){
      temp_data<-sample(which(!is.na(pheno_input[,pheno_col_num])),size=bs_sample_size,replace = T)

      if (!is.null(covar_path)) {
        l1t<-try(lm(scale(pheno_input[temp_data,pheno_col_num])~scale(prs_input[,prs_col_num])[temp_data]+scale(as.matrix(covar_input[temp_data,-1]))))
      } else {
        l1t<-try(lm(scale(pheno_input[temp_data,pheno_col_num])~scale(prs_input[,prs_col_num])[temp_data]))
      }
      if (class(l1t)== "try-error"){
        cat(paste0("Error in bootstrap for iter", ii," ...skip. \n"))
        next
      }
      if (class(prs_path)!= "try-error"){
        TEMP[ii]<-summary(l1t)$coef[2,1]*TAGC_bias
      }
    }
    ###Boot SE
    (summary2[1,5] <- tagc_bs_se <- sd(TEMP,na.rm = T))
    if (is.na(summary2[1,4])) {
      ##If boot not worksing, e.g., for trait 14
      cat(paste0("Error in bootstrap, use regression standard error for bootstrap standard error.\n"))
      summary2[1,5] <- tagc_bs_se <- tagc_se
    }
    ###Boot CI
    summary2[1,6]<-quantile(TEMP,prob=(1-CI_conf_level)/2,na.rm=T)
    summary2[1,7]<-quantile(TEMP,prob=1-(1-CI_conf_level)/2,na.rm=T)
    ###Boot CI with Normality
    bs_CI_normal_LB <- tagc + qnorm((1-CI_conf_level)/2)*tagc_bs_se
    bs_CI_normal_UB <- tagc + qnorm(1-(1-CI_conf_level)/2)*tagc_bs_se
    if (is.na(tagc_bs_se)) {
      ## If boot not worksing, e.g., for trait 14 (note this happens if TEMP contains only NA)
      warning('None of the bootstrap samples have a non-NA standard error. Using regression stndard error for confidence interval.\n')
      summary2[1,5]<-bs_CI_normal_LB
      summary2[1,6]<-bs_CI_normal_UB
    }
  }
  summary2_df = data.frame(summary2)
  colnames(summary2_df) = c("P_value","GC0","SE0",
                            "GC","GC_SE",
                            "GC_CI1","GC_CI2",
                            "N_obs")
  return(summary2_df)
}


