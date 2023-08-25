

TAGC_gc_se_CI_old = function(TAGC_bias, covar_path=NULL, pheno_path, pheno_col_num, prs_path, prs_col_num, number_bs_samples=500, bs_sample_size=1500, CI_conf_level=0.95) {

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

  summary2 = matrix(NA, 1, 14)
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

    print(summary2[1,1] <- summary(l1)$coef[2,4]) # p-value
    print(summary2[1,2] <- tagc <- orig_gc*TAGC_bias) # TAGC gc
    print(summary2[1,3] <- tagc_se <- summary(l1)$coef[2,2]*TAGC_bias) # regression se
    (summary2[1,14]<-nobs(l1)) # number of samples in lm
    #####
    ###bootstrap to get the SE and CI
    TEMP<-rep(NA, number_bs_samples)
    set.seed(1000)
    for(ii in 1:length(TEMP)){
      temp_data<-sample(which(!is.na(pheno_input[,pheno_col_num])),size=bs_sample_size,replace = T)
      l1t<-try(lm(scale(pheno_input[temp_data,pheno_col_num])~scale(prs_input[,prs_col_num])[temp_data]+scale(as.matrix(covar_input[temp_data,-1]))))
      if (class(l1t)== "try-error"){
        cat(paste0("Error in bootstrap for iter", ii," ...skip. \n"))
        next
      }
      if (class(prs_path)!= "try-error"){
        (TEMP[ii]<-summary(l1t)$coef[2,1]*TAGC_bias)
      }
    }
    ###Boot SE
    (summary2[1,4] <- tagc_bs_se <- sd(TEMP,na.rm = T))
    if (is.na(summary2[1,4])) {
      ##If boot not worksing, e.g., for trait 14
      cat(paste0("Error in bootstrap, use regression standard error for bootstrap standard error.\n"))
      summary2[1,4] <- tagc_bs_se <- tagc_se
    }
    ###Boot CI
    (summary2[1,5]<-quantile(TEMP,prob=(1-CI_conf_level)/2,na.rm=T))
    (summary2[1,6]<-quantile(TEMP,prob=1-(1-CI_conf_level)/2,na.rm=T))
    ###Boot CI with Normality
    (summary2[1,7] <- tagc + qnorm((1-CI_conf_level)/2)*tagc_bs_se)
    (summary2[1,8] <- tagc + qnorm(1-(1-CI_conf_level)/2)*tagc_bs_se)
    if (is.na(summary2[1,4])) {
      ## If boot not worksing, e.g., for trait 14 (note this happens if TEMP contains only NA)
      warning('None of the bootstrap samples have a non-NA standard error. Using regression stndard error for confidence interval.\n')
      (summary2[1,5]<-summary2[1,7])
      (summary2[1,6]<-summary2[1,8])
    }
    ####FZ Boot SE
    (summary2[1,9]<-sd(DescTools::FisherZ(TEMP)[!is.na(DescTools::FisherZ(TEMP))]))
    ####FZ CI with bootstrap n
    temp_ci<-CorCI_fz(tagc, n=dim(pheno_input)[1],
                       sigma0 = summary2[1,9],
                       conf.level = CI_conf_level, alternative = c("two.sided"))
    (summary2[1,10]<-temp_ci[2])
    (summary2[1,11]<-temp_ci[3])
    ####FZ CI with the theoretical n
    temp_ci2<-DescTools::CorCI(tagc, n=max(10,(1/tagc_bs_se)^2),
                    conf.level = CI_conf_level, alternative = c("two.sided"))
    (summary2[1,12]<-temp_ci2[2])
    (summary2[1,13]<-temp_ci2[3])
  }
  summary2_df = data.frame(summary2)
  colnames(summary2_df) = c("P_value","GC","SE",
                            "Boot_SE","Boot_CI1","Boot_CI2","Boot_N_CI1","Boot_N_CI2",
                            "FZ_Boot_SE","FZ_Boot_CI1","FZ_Boot_CI2","FZ_Boot_N_CI1","FZ_Boot_N_CI2","N_obs")
  return(summary2_df)
}


