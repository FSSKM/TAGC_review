
# code/TAGC/TAGC/TAGC_bias_factor_v2.R

# mom1, mom2, mom3 are unnormalized moments, output from LD_block_moment.R. mom1 = p
bias_factor = function(mom1, mom2, mom3, hb_2, ha_2, n) {
  b1<-mom1/mom1 #1
  b2<-mom2/mom1 #3.176597
  b3<-mom3/mom1 #82.82365
  if (b1 < 0 | b2 < 0 | b3 < 0) {
    stop('The given LD moments are negative, please check correctness of LD moments.')
  }
  if (any(hb_2 <= 0) | any(ha_2 <= 0)) {
    warning('The given heritablities are negative, please check correctness of heritabilities.')
  }
  bias<-sqrt((b3/b2^2+(mom1/n)/(b2*hb_2))/ha_2) # formula in the eq at the end of section 3.1
  return(bias)
}






