
# code/TAGC/LD_block_moment/LD_block_moment_each_WAGC.R



### compute LD matrix moment for each LD block, WAGC version
WAGC_LD_block_moment = function(bfile_path, bfile_name, popn_list_path=NULL, output_path=NULL) {

  if (is.null(bfile_path)) {
    stop(paste0('bfile_path is not provided. Please provide the path to the folder containing LD block bfiles.'))
  }
  if (is.null(bfile_name)) {
    stop(paste0('bfile_name is not provided. Please provide the name of the LD block bfile without extension.'))
  }
  if (is.null(popn_list_path)) {
    warning(paste0('popn_list_path is not provided. Will use all subjects in bfile_name, ', bfile_name, '.'))
  }
  if (is.null(output_path)) {
    warning(paste0('output_path is not provided. Will use the current working directory as the output path.'))
    output_path = '.'
  }

  # create a sub-directory for the output
  output_path = paste0(output_path, '/', 'WAGC_LD_block_moment')
  dir.create(output_path, recursive=T, showWarnings=F)
  # name for the output file
  output_filename = paste0(bfile_name, '_traces')

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

  ### load subjects
  # popn_list = read.table(popn_list_path, header=F)
  if (!is.null(popn_list_path)) {
    popn_list = read.table(popn_list_path, header=F)
  } else {
    popn_list = read.table(paste0(bfile_path, '/', bfile_name, '.fam'), header=F)
  }

  ### load LD block
  setwd(bfile_path)
  cat(paste0("\nProcessing LD Block ", bfile_name, '\n'))
  data<-try(BEDMatrix::BEDMatrix(bfile_name))

  if (class(data) != "try-error") {
    data_matrix<-as.matrix(data)
    data_matrix2<-data_matrix
    dim(data_matrix2)
    data_matrix3<-data_matrix2[rownames(data_matrix2)%in%paste(popn_list$V1,popn_list$V2,sep="_"),]
    data_matrix3<-Input_by_mean(data_matrix3)$Y
    list00<-unique(c(which(apply(data_matrix3,2,sd)==0)))
    if(length(list00)>0){
      data_matrix2<-data_matrix2[,-list00]
      cat(paste0(length(list00),"SNPs...removed due to zero sd.\n"))
    }
    output<-matrix(NA,4,1)
    ######### popn1
    data_matrix3<-data_matrix2[rownames(data_matrix2)%in%paste(popn_list$V1,popn_list$V2,sep="_"),]
    data_matrix3<-Input_by_mean(data_matrix3)$Y
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
    cat(paste(dim(data_matrix3)[2],"SNPs...done trace\n"))
    #########
    write(t(output),
          paste0(output_path, '/', output_filename, '.txt'),
          ncolumns=ncol(output))
    cat(paste0('LD block moments written to ', paste0(output_path, '/', output_filename, '.txt'), '.\n'))
    return(paste0(output_path, '/', output_filename, '.txt'))
  } else {
    stop(paste0('Cannot load the bfile ', bfile_path, '/', bfile_name, '. Please make sure the bfile exists.\n'))
  }
}



