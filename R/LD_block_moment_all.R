

### combine the moments of each LD block matrix, TAGC version
# n1 is the sample size of population I
LD_block_moment_all_TAGC = function(LD_block_mom_path, n1, output_path='.', output_filename=NULL) {
  if (is.null(output_filename)) {
    output_filename = paste0('TAGC_LD_block_mom1_mom2_mom3')
  }
  output_path = paste0(output_path, '/', 'TAGC_moment')
  dir.create(output_path, recursive=T, showWarnings=F)

  # a list of file names that are output from TAGC_LD_block_moment or WAGC_LD_block_moment
  file_lst = list.files(LD_block_mom_path)
  df_lst = lapply(paste0(LD_block_mom_path, '/', file_lst), FUN=read.table)

  if (ncol(df_lst[[1]]) != 3) {
    stop('Number of columns in the txt files in the provided LD_block_mom_path should be 3. Please check if LD_block_mom_path is correct. It is also possible you should call LD_block_moment_all_WAGC instead.')
  }

  temp1 = sum(sapply(df_lst, FUN=function(x){x[1,3]}))
  temp2 = sum(sapply(df_lst, FUN=function(x){x[2,3]}))
  temp3 = sum(sapply(df_lst, FUN=function(x){x[3,3] - 1/n1*x[1,3]*x[2,3]}))

  # output un-normalized version so that 'temp1' gives the value of 'p'
  write(c(temp1, temp2, temp3), file=paste0(output_path, '/', output_filename, '.txt'))
  cat(paste0('Output written to ', paste0(output_path, '/', output_filename, '.txt'), '\n'))
  return(c(temp1, temp2, temp3))
}


### combine the moments of each LD block matrix, TAGC version
# n1 is the sample size of the population
LD_block_moment_all_WAGC = function(LD_block_mom_path, n1, output_path='.', output_filename=NULL) {
  if (is.null(output_filename)) {
    output_filename = paste0('WAGC_LD_block_mom1_mom2_mom3')
  }
  output_path = paste0(output_path, '/', 'WAGC_moment')
  dir.create(output_path, recursive=T, showWarnings=F)

  # a list of file names that are output from TAGC_LD_block_moment or WAGC_LD_block_moment
  file_lst = list.files(LD_block_mom_path)
  df_lst = lapply(paste0(LD_block_mom_path, '/', file_lst), FUN=read.table)

  if (ncol(df_lst[[1]]) != 1) {
    stop('Number of columns in the txt files in the provided LD_block_mom_path should be 1. Please check if LD_block_mom_path is correct. It is also possible you should call LD_block_moment_all_TAGC instead.')
  }
  temp1 = sum(sapply(df_lst, FUN=function(x){x[2,1]}))
  temp2 = sum(sapply(df_lst, FUN=function(x){x[3,1] - x[2,1]^2/n1}))
  temp3 = sum(sapply(df_lst, FUN=function(x){x[4,1] - x[2,1]^3/n1^2 - 3*x[2,1]*(x[3,1] - x[2,1]^2/n1)/n1 }))

  # output un-normalized version so that 'temp1' gives the value of 'p'
  write(c(temp1, temp2, temp3), file=paste0(output_path, '/', output_filename, '.txt'))
  cat(paste0('Output written to ', paste0(output_path, '/', output_filename, '.txt'), '\n'))
  return(c(temp1, temp2, temp3))
}


