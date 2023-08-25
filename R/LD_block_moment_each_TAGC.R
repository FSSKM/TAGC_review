
# code/TAGC/LD_block_moment/LD_block_moment_each_TAGC.R



### compute LD matrix moment for each LD block, TAGC version
TAGC_LD_block_moment = function(bfile_path_popn1, bfile_path_popn2=NULL, bfile_name_popn1, bfile_name_popn2=NULL, 
  popn1_list_path=NULL, popn2_list_path=NULL, output_path=NULL, output_filename=NULL) {

  if (is.null(bfile_path_popn1)) {
    stop(paste0('bfile_path_popn1 is not provided. Please provide the path to the folder containing LD block bfiles.'))
  }
  if (is.null(bfile_path_popn2)) {
    warning(paste0('bfile_path_popn2 is not provided. Will try using bfile_path_popn1, ', bfile_path_popn1, ', if the bfiles for both populations are in ', bfile_path_popn1, '.'))
  }
  if (is.null(bfile_name_popn1)) {
    stop(paste0('bfile_name_popn1 is not provided. Please provide the name of the LD block bfile without extension.'))
  }
  if (is.null(bfile_name_popn2)) {
    warning(paste0('bfile_path_popn2 is not provided. Will try using bfile_name_popn1, ', bfile_name_popn1, ', if the subjects from both populations are in ', bfile_name_popn1, '.'))
  }
  if (is.null(popn1_list_path)) {
    if (is.null(bfile_name_popn2)) {
      stop(paste0('Please either provide the popn1_list_path or ensure two distinct bfiles (bfile_name_popn1 and bfile_name_popn2) are provided.'))
    } else {
      warning(paste0('popn1_list_path is not provided. Will use all subjects in bfile_name_popn1, ', bfile_name_popn1, '.'))
    }
  }
  if (is.null(popn2_list_path)) {
    warning(paste0('popn2_list_path is not provided. Will use all subjects in bfile_name_popn2, ', bfile_name_popn2, ' if bfile_name_popn2 is provided.'))
    if (is.null(bfile_name_popn2)) {
      stop(paste0('Neither popn2_list_path nor bfile_name_popn2 is provided. If you intend to use one set of subjects for LD block moment computation, you should reconsider if it is a within-ancestry setup and use WAGC_LD_block_moment instead.'))
    }
  }
  if (is.null(output_path)) {
    warning(paste0('output_path is not provided. Will use the current working directory as the output path.'))
    output_path = '.'
  }

  ### helper code
  Input_by_mean<-function (Y){
    check<-rep(NA,dim(Y)[2])
    for ( i in 1:dim(Y)[2]){
      check[i]<-sum(is.na(Y[,i])*1)
      if(check[i]!=0){
        temp1<-mean(Y[,i],na.rm = T)
        Y[which(is.na(Y[,i])==T),i]<-temp1
      }
    }
    names(check)<-colnames(Y)
    out<-list(check=check,Y=Y)
    return(out)
  }

  # create a sub-directory for the output
  output_path = paste0(output_path, '/', 'TAGC_LD_block_moment')
  dir.create(output_path, recursive=T, showWarnings=F)

  ### load subjects
  if (!is.null(popn1_list_path)) {
    popn1_list = read.table(popn1_list_path, header=F)
  } else {
    popn1_list = read.table(paste0(bfile_path_popn1, '/', bfile_name_popn1, '.fam'), header=F)
  }
  if (!is.null(popn2_list_path)) {
    popn2_list = read.table(popn2_list_path, header=F)
  } else {
    if (is.null(bfile_path_popn2)) {
      popn2_list = read.table(paste0(bfile_path_popn1, '/', bfile_name_popn2, '.fam'), header=F)
    } else {
      popn2_list = read.table(paste0(bfile_path_popn2, '/', bfile_name_popn2, '.fam'), header=F)
    }
  }
  
  # a quick check if in fact popn1 and popn2 are identical (in this case, should not use TAGC, but WAGC)
  if (all(popn1_list[,1] %in% popn2_list[,1])) {
    warning(paste0('All of popn 1 subjects are in popn 2. Please consider applying WAGC (WAGC_LD_block_moment) if popn 1 and popn 2 are identical.\n'))
  } else if (all(popn2_list[,1] %in% popn1_list[,1])) {
    warning(paste0('All of popn 2 subjects are in popn 1. Please consider applying WAGC (WAGC_LD_block_moment) if popn 1 and popn 2 are identical.\n'))
  }

  # if user only provided one bfile (bfile_name_popn1)
  if (is.null(bfile_name_popn2)) {
    # name for the output file
    if (is.null(output_filename)) {
      output_filename = paste0(bfile_name_popn1, '_TAGC_traces')
    }

    # load LD block of both populations from bfile_name_popn1
    setwd(bfile_path_popn1)
    cat(paste0("\nProcessing Population 1 and 2 LD Block ", bfile_name_popn1, '\n'))
    data<-try(BEDMatrix::BEDMatrix(bfile_name_popn1))

    if (class(data) != "try-error") {
      data_matrix<-as.matrix(data)
      data_matrix2<-data_matrix
      dim(data_matrix2)
      data_matrix3<-data_matrix2[rownames(data_matrix2)%in%paste(popn1_list$V1,popn1_list$V2,sep="_"),]
      data_matrix3<-Input_by_mean(data_matrix3)$Y
      if (nrow(data_matrix3) == 0) {
        stop(paste0('Zero subjects remaining in the bfile after matching with the subjects in popn1_list_path and LD moment cannot be computed.'))
      }
      data_matrix3a<-data_matrix2[rownames(data_matrix2)%in%paste(popn2_list$V1,popn2_list$V2,sep="_"),]
      data_matrix3a<-Input_by_mean(data_matrix3a)$Y
      if (nrow(data_matrix3a) == 0) {
        stop(paste0('Zero subjects remaining in the bfile after matching with the subjects in popn2_list_path and LD moment cannot be computed.'))
      }
    } else {
      stop(paste0('Cannot load the bfile ', bfile_path_popn1, '/', bfile_name_popn1, '. Please make sure the bfile exists.\n'))
    }
  } else { # if user provided two bfiles
    # name for the output file
    if (is.null(output_filename)) {
      output_filename = paste0(bfile_name_popn1, '_', bfile_name_popn2, '_TAGC_traces')
    }

    # load LD block of population 1 from bfile_name_popn1
    setwd(bfile_path_popn1)
    cat(paste0("\nProcessing Population 1 LD Block ", bfile_name_popn1, '\n'))
    data<-try(BEDMatrix::BEDMatrix(bfile_name_popn1))
    if (class(data) != "try-error") {
      data_matrix<-as.matrix(data)
      data_matrix2<-data_matrix
      dim(data_matrix2)
      data_matrix3<-data_matrix2[rownames(data_matrix2)%in%paste(popn1_list$V1,popn1_list$V2,sep="_"),]
      data_matrix3<-Input_by_mean(data_matrix3)$Y
      if (nrow(data_matrix3) == 0) {
        stop(paste0('Zero subjects remaining in the bfile after matching with the subjects in popn1_list_path and LD moment cannot be computed.'))
      }
    } else {
      stop(paste0('Cannot load the bfile ', bfile_path_popn1, '/', bfile_name_popn1, '. Please make sure the bfile exists.\n'))
    }

    # load LD block of population 2 from bfile_name_popn2
    if (!is.null(bfile_path_popn2)) {
      setwd(bfile_path_popn2)
    }
    cat(paste0("\nProcessing Population 2 LD Block ", bfile_name_popn2, '\n'))
    data<-try(BEDMatrix::BEDMatrix(bfile_name_popn2))
    if (class(data) != "try-error") {
      data_matrix<-as.matrix(data)
      data_matrix2<-data_matrix
      dim(data_matrix2)
      data_matrix3a<-data_matrix2[rownames(data_matrix2)%in%paste(popn2_list$V1,popn2_list$V2,sep="_"),]
      data_matrix3a<-Input_by_mean(data_matrix3a)$Y
      if (nrow(data_matrix3a) == 0) {
        stop(paste0('Zero subjects remaining in the bfile after matching with the subjects in popn2_list_path and LD moment cannot be computed.'))
      }
    } else {
      stop(paste0('Cannot load the bfile ', getwd(), '/', bfile_name_popn2, '. Please make sure the bfile exists.\n'))
    }
  }

  # compute LD matrix and LD moment
  list00<-unique(c(which(apply(data_matrix3,2,sd)==0), which(apply(data_matrix3a,2,sd)==0)))
  if(length(list00)>0){
    # data_matrix2<-data_matrix2[,-list00]
    data_matrix3 = data_matrix3[,-list00]
    data_matrix3a = data_matrix3a[,-list00]
    cat(paste0(length(list00),"SNPs...removed due to zero sd.\n"))
  }
  output<-matrix(NA,4,3)
  ######### popn1
  dim(data_matrix3)
  ######### compute LD matrix of popn1
  cor_temp<-cor(data_matrix3,use = "complete")
  cor_temp2<-cor_temp%*%cor_temp
  cor_temp3<-cor_temp2%*%cor_temp
  output[1,1]<-dim(data_matrix3)[2]
  output[2,1]<-psych::tr(cor_temp)
  output[3,1]<-psych::tr(cor_temp2)
  output[4,1]<-psych::tr(cor_temp3)
  #########
  ######### popn2
  dim(data_matrix3a)
  ######### compute LD matrix of popn2
  cor_tempa<-cor(data_matrix3a,use = "complete")
  cor_temp2a<-cor_tempa%*%cor_tempa
  cor_temp3a<-cor_temp2a%*%cor_tempa
  output[1,2]<-dim(data_matrix3a)[2]
  output[2,2]<-psych::tr(cor_tempa)
  output[3,2]<-psych::tr(cor_temp2a)
  output[4,2]<-psych::tr(cor_temp3a)
  ######### popn1 to popn2
  cor_tempc1<-cor_temp%*%cor_tempa
  cor_tempc2<-cor_temp2%*%cor_tempa
  output[1,3]<-dim(data_matrix3)[2]
  output[2,3]<-psych::tr(cor_tempc1)
  output[3,3]<-psych::tr(cor_tempc2)
  #########
  cat(paste(dim(data_matrix3)[2],"SNPs...done trace\n"))
  #########
  write(t(output),
        paste0(output_path, '/', output_filename, '.txt'),
        ncolumns=ncol(output))
  cat(paste0('LD block moments written to ', paste0(output_path, '/', output_filename, '.txt'), '.\n'))
  return(paste0(output_path, '/', output_filename, '.txt'))
}


