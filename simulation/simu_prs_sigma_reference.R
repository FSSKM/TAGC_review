library(psych)
library(MASS)
library(bdsmatrix)
library(MatrixEQTL)
library(glmnet)
library(e1071)
library(data.table)
data0<-fread("TAGC/LD-EUR-CHR22-hapmap3.csv")
dim(data0)
data0<-data.frame(data0)
data0<-data0[1:10000,1:10000]
head(data0[,1:10])
##
data00<-fread("TAGC/LD-EAS-CHR22-hapmap3.csv")
dim(data00)
data00<-data.frame(data00)
data00<-data00[1:10000,1:10000]
head(data00[,1:10])
#########################################################################################################
##################################the SNP data###########################################################
#########################################################################################################
n<-dim(data0)[1]
n.run<-100

args=(commandArgs(TRUE))
prop_signal<-as.numeric(args[1])
prop_h2<-as.numeric(args[2])
gc_corr<-as.numeric(args[3])
iter<-as.numeric(args[4])
E2<-R2<-matrix(NA,20,4)
nw<-5000
nz<-500
n.run<-100
####################
size<-1
for(size in c(1)){
  Sigma<-as.matrix(data0)
  Sigma2<-as.matrix(data00)
  dim(Sigma)
  dim(Sigma2)
  p<-dim(Sigma)[1]
  (omega<-p/n)
  head(Sigma[,1:20])
  ###############
  load("TAGC/LD-EUR-CHR22-hapmap3-moments-v1.RData")
  gmoments<-read.csv("TAGC/LD-EUR-CHR22-hapmap3-moments-v1.csv")
  gmoments<-gmoments[,-1]
  ###############

  (beta1<-gmoments[1,1])
  (beta2<-gmoments[1,2])
  (beta3<-gmoments[1,3])
  
  #######Use the B-N model########
  set.seed(123456+100000*size+2000*iter+200) 
  X0<-Z0<-W0a<-W0b<-matrix(NA,n,p)
  for(ii in 1:p){
    mfa1<-runif(n=1, min = 0.05, max = 0.45)
    ###
    fst1<-runif(n=1, min =-0.04, max = 0.04)
    fst2<-runif(n=1, min =-0.04, max = 0.04)
    ###
    mfa2<-mfa1+fst1
    mfa3<-mfa1+fst2
    ###
    X0[1:n,ii]<-rbinom(n,2,mfa2)
    Z0[1:n,ii]<-rbinom(n,2,mfa3)
    W0a[1:n,ii]<-rbinom(n,2,mfa2)
    W0b[1:n,ii]<-rbinom(n,2,mfa3)
    ###
  }
  X<-X0%*%Sigmah
  Z<-Z0%*%Sigmah2
  Wa<-W0a%*%Sigmah
  Wb<-W0b%*%Sigmah2
  ##handle some rare cases (if any)
  X[,which((apply(X, 2, var))==0)]<-rbinom(n,2,0.5) 
  X<-scale(X)
  Z[,which((apply(Z, 2, var))==0)]<-rbinom(n,2,0.5) 
  Z<-scale(Z)
  Z<-Z[1:nz,]
  Wa[,which((apply(Wa, 2, var))==0)]<-rbinom(n,2,0.5) 
  Wa<-scale(Wa)
  Wa<-Wa[1:nw,]
  Wb[,which((apply(Wb, 2, var))==0)]<-rbinom(n,2,0.5) 
  Wb<-scale(Wb)
  Wb<-Wb[1:nw,]
  #######projection matrix########
  dim(X)
  dim(Z)
  dim(Wa)
  dim(Wb)
  
  PWa<-t(Wa)%*%Wa
  PWb<-t(Wb)%*%Wb
  PWc<-(t(Wa[1:floor(nw/2),])%*%Wa[1:floor(nw/2),]+
          t(Wb[1:floor(nw/2),])%*%Wb[1:floor(nw/2),])
  #####
  sigma<-c(1,(1-prop_h2[1])/prop_h2[1])
  lambda0<-as.numeric(omega*(1-prop_h2[1])/prop_h2[1])
  #######dataset##########
  p.causal<-p*(prop_signal[1])
  #####
  error.train<-rnorm(n,mean=0,sd=sqrt(sigma[2]))
  error.test<-rnorm(nz,mean=0,sd=sqrt(sigma[2]))
  ##### genetic correlation
  Sigma.beta<- matrix(c(sigma[1],gc_corr,
                        gc_corr,sigma[1]),ncol=2)
  coef<-mvrnorm(n=p.causal,
                mu=c(0,0),
                Sigma=1*Sigma.beta)
  cor(coef)
  dim(coef)
  ####
  beta1<-as.matrix(coef[,1])
  beta2<-as.matrix(coef[,2])
  index<-sample((1:p),size=(p.causal))
  ##################################
  signal.train<-scale(X[,index]%*%beta1)
  signal.test<-scale(Z[,index]%*%beta2)
  ##################################
  y.train<-signal.train+error.train
  y.test<-signal.test+error.test
  #####
  dim(y.train)
  dim(y.test)
  ###################################
  ###################################
  (h2<-var(signal.train)/var(y.train))
  var(signal.test)/var(y.test)
  
  #####
  ################################################################
  ###################1.Marginal###################################
  ################################################################
  ####Create three SlicedData objects for the analysis
  snps1 <- SlicedData$new(t(X));
  gene1 <- SlicedData$new(t(y.train));
  cvrt1 <- SlicedData$new();
  ####Produce no output files
  filename <- NULL; 
  ####Call the main analysis function
  m.temp1 <- Matrix_eQTL_main(
    snps = snps1, 
    gene = gene1, 
    cvrt = cvrt1, 
    output_file_name = filename, 
    pvOutputThreshold = 1, 
    useModel = modelLINEAR, 
    errorCovariance = numeric(), 
    verbose = F,
    pvalue.hist = FALSE,
    min.pv.by.genesnp=F);
  ####Pull Matrix eQTL results
  pvalue.temp<-m.temp1$all$eqtls$pvalue
  snps.temp0<-m.temp1$all$eqtls$snps
  snps.temp<-as.numeric(substring(snps.temp0, 4))
  beta.temp<-m.temp1$all$eqtls$beta
  eqtl.temp<-data.frame(pvalue.temp,snps.temp,beta.temp)
  head(eqtl.temp)
  ######out-of-sample######
  test.pred1<-Z[,eqtl.temp$snps.temp]%*%as.matrix(eqtl.temp$beta.temp)
  (R2[1:4,size]<-as.vector(summary(lm(scale(y.test)~scale(test.pred1)))$coef[2,]))
  print(R2[5,size]<-(R2[1,size])*sqrt((gmoments[3,2]/gmoments[2,2]^2+omega/(prop_h2*gmoments[2,2]))/prop_h2))
  ################################################
  ################Reference panel:v1##############
  ################################################
  lambda<-as.numeric(omega*(1-h2)/h2)
  Ridge.P<-solve(PWa+diag(lambda*nw,p))
  ridge.coef<-Ridge.P%*%t(X)%*%(y.train)
  ########out-of-sample########
  test.pred2<-Z%*%as.matrix(ridge.coef)
  (R2[6:9,size]<-as.vector(summary(lm(scale(y.test)~scale(test.pred2)))$coef[2,]))
  ########
  Rw0<-Ridge.P*nw
  Rw1<-prop_h2*(tr((Sigma2)%*%Rw0%*%(Sigma)))^2/p
  Rw2<-prop_h2*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma)%*%(Sigma))
  Rw3<-omega*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma))
  print(R2[10,size]<-(R2[6,size])/sqrt(prop_h2*Rw1/(Rw2+Rw3)))
  ################################################
  ################Reference panel:v2##############
  ################################################
  Ridge.P<-solve(PWb+diag(lambda*nw,p))
  ridge.coef<-Ridge.P%*%t(X)%*%(y.train)
  ########out-of-sample########
  test.pred3<-Z%*%as.matrix(ridge.coef)
  (R2[11:14,size]<-as.vector(summary(lm(scale(y.test)~scale(test.pred3)))$coef[2,]))
  Rw0<-Ridge.P*nw
  Rw1<-prop_h2*(tr((Sigma2)%*%Rw0%*%(Sigma)))^2/p
  Rw2<-prop_h2*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma)%*%(Sigma))
  Rw3<-omega*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma))
  print(R2[15,size]<-(R2[11,size])/sqrt(prop_h2*Rw1/(Rw2+Rw3)))
  ################################################
  ################Reference panel:v3##############
  ################################################
  Ridge.P<-solve(PWc+diag(lambda*nw,p))
  ridge.coef<-Ridge.P%*%t(X)%*%(y.train)
  ########out-of-sample########
  test.pred4<-Z%*%as.matrix(ridge.coef)
  (R2[16:19,size]<-as.vector(summary(lm(scale(y.test)~scale(test.pred4)))$coef[2,]))
  Rw0<-Ridge.P*nw
  Rw1<-prop_h2*(tr((Sigma2)%*%Rw0%*%(Sigma)))^2/p
  Rw2<-prop_h2*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma)%*%(Sigma))
  Rw3<-omega*tr(Rw0%*%(Sigma2)%*%Rw0%*%(Sigma))
  print(R2[20,size]<-(R2[16,size])/sqrt(prop_h2*Rw1/(Rw2+Rw3)))
  ################################
}
####################################
####################################
write.csv(R2,
          paste0("TAGC/results/R2/TAGC-Simu-PRS-June17-2023-1-ref-sigma-m",prop_signal[1],"-h",prop_h2[1],"-corr",gc_corr[1],"-iter",iter,".csv"),
          row.names = F)
####################################
####################################
save(eqtl.temp,Z,y.test,y.train,test.pred1,
     file = paste0("TAGC/results/data/TAGC-Simu-PRS-June17-2023-1-ref-sigma-m",prop_signal[1],"-h",prop_h2[1],"-corr",gc_corr[1],"-iter",iter,".RData"))






##################################
##################################
##################################
##################################
rm(list=ls())
##4 n/p ratio
##6, h2
##6, m/p
##10 methods
##100 replicates
CORR_ALL0<-CORR_ALL<-array(NA,c(3,2,2,4,200))
l<-1;
j<-1;
k<-1;
h<-1;
jlist<-c(0.05, 0.3)
hlist<-c(0.3, 0.6)
klist<-c(0, 0.45, 0.9)
for(l in 1:200){print(l)
  for(j in 1:2){
    for(h in 1:2){
      for(k in 1:3){
        ###############
        truth<-try(read.csv(paste0("TAGC/results/R2/TAGC-Simu-PRS-June17-2023-1-ref-sigma-m",jlist[j],"-h",hlist[h],"-corr",klist[k],"-iter",l,".csv"),header = T))
        ###############
        if (class(truth) != "try-error"){
          CORR_ALL[k,j,h,,l]<-as.matrix(truth[c(5,10,15,20),1])
          CORR_ALL0[k,j,h,,l]<-as.matrix(truth[c(1,6,11,16),1])
        }
        ###############
      }
    }
  }
}
#############
apply(CORR_ALL[3,2,2,,],1,mean,na.rm=T)
apply(CORR_ALL0[3,2,2,,],1,mean,na.rm=T)
apply(CORR_ALL[2,2,2,,],1,mean,na.rm=T)
apply(CORR_ALL0[2,2,2,,],1,mean,na.rm=T)
apply(CORR_ALL[1,2,2,,],1,mean,na.rm=T)
apply(CORR_ALL0[2,2,2,,],1,mean,na.rm=T)
mean(CORR_ALL[3,,,,],na.rm=T)
mean(CORR_ALL0[3,,,,],na.rm=T)
mean(CORR_ALL[2,,,,],na.rm=T)
mean(CORR_ALL0[2,,,,],na.rm=T)
mean(CORR_ALL[1,,,,],na.rm=T)
mean(CORR_ALL0[1,,,,],na.rm=T)
#############
save(CORR_ALL,CORR_ALL0,
     file="TAGC/results/Corr-Simu-PRS-June17-2023-1-ref-sigma.RData")
####################################
####################################

load("Pure_simulation/results/Corr-Simu-PRS-June17-2023-1-ref-sigma.RData")
R2<-(CORR_ALL)
dim(R2)
###################################
###################################
cols<-c('sienna','palevioletred1','royalblue2','darkturquoise',
                'darkviolet','brown','midnightblue','dimgrey','darkorange1',
                "forestgreen","blueviolet",'orange','yellow','green','red')
                
names<-c("Marginal","Ref-X","Ref-Z","Ref-Mixed")              
n<-10000
pall<-c(0, 0.45, 0.9)
mp<-c(0.05,0.3)
sigma<-c(0.3, 0.6)
np<-c(1,1,1,1)
##################
pdf("Pure_simulation/results/Corr-Simu-PRS-June17-2023-1-ref-sigma-1.pdf",width=12,height=4.5)
par(mar=c(6.5,4,2,2), xpd=F)
par(mfrow=c(1,3))
par(cex.axis=1.25)
#options(scipen=100000)
l<-1
j<-1
R22<-array(NA,c(3,4,(2*2*200)))
for(l in c(1:3)){
    for(j in 1:4){
      R22[l,j,]<-as.vector(R2[l,,,j,])
      ####
    }
}
dim(R22)
for(l in 1:3){
  boxplot(t(R22[l,,]), las = 2,
          col =cols[1:4],ylim=c(-2,2),
          names =names,cex=1.05,
          main=paste())
  mtext('Genetic orrelation', side = 2, line = 2.5,cex = 1.05)
  mtext(paste("True genetic correlation=", pall[l],"",sep = ""), side = 3, line = 0.5,cex=1.05)
  grid(col = "lightgray", lty = "dotted",lwd=2,equilogs = TRUE)
  abline(h=pall[l],lty=2.5,lwd=2,col="2")
  points(y=apply(t(R22[l,,]),2,mean,na.rm=T),x=c(1:length(names)),pch=20)
}
dev.off()
##################################### 

####################################
####################################
load("Pure_simulation/results/Corr-Simu-PRS-June17-2023-1-ref-sigma.RData")
R2<-(CORR_ALL0)
dim(R2)
###################################
###################################
cols<-c('sienna','palevioletred1','royalblue2','darkturquoise',
                'darkviolet','brown','midnightblue','dimgrey','darkorange1',
                "forestgreen","blueviolet",'orange','yellow','green','red')
                
names<-c("Marginal","Ref-X","Ref-Z","Ref-Mixed")              

n<-10000
pall<-c(0, 0.45, 0.9)
mp<-c(0.05,0.3)
sigma<-c(0.3, 0.6)
np<-c(1,1,1,1)
##################
pdf("Pure_simulation/results/Corr-Simu-PRS-June17-2023-1-ref-sigma-2.pdf",width=12,height=4.5)
par(mar=c(6.5,5,2,2), xpd=F)
par(mfrow=c(1,3))
par(cex.axis=1.25)

l<-1
j<-1
R22<-array(NA,c(3,4,(2*2*200)))
for(l in c(1:3)){
  for(j in 1:4){
    R22[l,j,]<-as.vector(R2[l,,,j,])
    ####
  }
}
for(l in 1:3){
  boxplot(t(R22[l,,]), las = 2,
          col =cols[1:4],ylim=c(-1,1),
          names =names,cex=1.05,
          main=paste())
  mtext('Genetic orrelation', side = 2, line = 3.5,cex = 1.05)
  mtext(paste("True genetic correlation=", pall[l],"",sep = ""), side = 3, line = 0.5,cex=1.05)
  grid(col = "lightgray", lty = "dotted",lwd=2,equilogs = TRUE)
  abline(h=pall[l],lty=2.5,lwd=2,col="2")
  points(y=apply(t(R22[l,,]),2,mean,na.rm=T),x=c(1:length(names)),pch=20)
}
dev.off()
##################################### 
