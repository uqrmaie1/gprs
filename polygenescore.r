# polygenescore: power, area under curve and correlation for estimated gene scores
# Frank Dudbridge, December 2012
#
#
# INPUTS
# n1 = discovery sample size
# nsnp = number of independent SNPs in the gene score
# vg1 = proportion of total variance that is explained by genetic effects in discovery sample
# n2 = replication sample size
# vg2 = proportion of total variance that is explained by genetic effects in replication sample
# corr = correlation between genetic effect sizes in the two populations
# plower = lower bound on p-value for selection from discovery sample
# pupper = upper bound on p-value for selection from discovery sample
# weighted = T if effect sizes used as weights in forming gene score
# alpha = type-1 error for testing association in replication sample
# binary = T for binary traits
# prevalence1 = disease prevalence in discovery sample
# prevalence2 = disease prevalence in replication sample
# sampling1 = case/control sampling fraction in discovery sample
# sampling2 = case/control sampling fraction in replication sample
# lambdaS1 = sibling relative recurrence risk in discovery sample, can be specified instead of vg1
# lambdaS2 = sibling relative recurrence risk in replication sample, can be specified instead of vg2
# nullfraction = proportion of SNPs with no effects
# shrinkage = T if effect sizes are to be shrunk to BLUPs
# logrisk = T if binary trait arises from log-risk model, not liability threshold

# OUTPUTS
# R2 = squared correlation between estimated gene score and replication trait
# NCP = non-centrality parameter of association test between score and replication trait
# p = expected p-value of association test
# power = power of association
# AUC = for binary traits, area under ROC curve
# MSE = for quantitative traits, mean square error between replication trait and estimated gene score

polygenescore=function(n1,nsnp,vg1=1,n2=n1,vg2=vg1,corr=1,plower=0,pupper=1,weighted=T,alpha=0.05,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {

# variance of the phenotype
varY1 = 1
varY2 = 1
if (binary) {
  varY1 = sampling1*(1-sampling1)
  varY2 = sampling2*(1-sampling2)
}

# sampling variance of each beta
samplingVar = varY1/n1

# conversion from lambdaS to vg
if (!is.na(lambdaS1)) {
  if (logrisk) {
    vg1=2*log(lambdaS1)
  }
  else {
    t=-qnorm(prevalence1)
    t1=qnorm(1-lambdaS1*prevalence1)
    vg1 = 2*(t-t1*sqrt(1-(t^2-t1^2)*(1-t*prevalence1/dnorm(t))))/(dnorm(t)/prevalence1+t1^2*(dnorm(t)/prevalence1-t))
    if (vg1>1 | is.na(vg1)) vg1=1
  }
}
if (!is.na(lambdaS2)) {
  if (logrisk) {
    vg2=2*log(lambdaS2)
  }
  else {
    t=-qnorm(prevalence2)
    t1=qnorm(1-lambdaS2*prevalence2)
    vg2 = 2*(t-t1*sqrt(1-(t^2-t1^2)*(1-t*prevalence2/dnorm(t))))/(dnorm(t)/prevalence2+t1^2*(dnorm(t)/prevalence2-t))
    if (vg2>1 | is.na(vg2)) vg2=1
  }
}
# variance of true betas
betaVar = vg1/(nsnp*(1-nullfraction))
betaVar2 = vg2/(nsnp*(1-nullfraction))

# transform from liability scale to observed scale
if (logrisk) {
  liab2obs1=prevalence1*sampling1*(1-sampling1)/prevalence1/(1-prevalence1)
  liab2obs2=prevalence2*sampling2*(1-sampling2)/prevalence2/(1-prevalence2)
}
else {
  liab2obs1=dnorm(qnorm(prevalence1))*sampling1*(1-sampling1)/prevalence1/(1-prevalence1)
  liab2obs2=dnorm(qnorm(prevalence2))*sampling2*(1-sampling2)/prevalence2/(1-prevalence2)
}
if (binary) betaVar = betaVar*liab2obs1^2
if (binary) betaVar2 = betaVar2*liab2obs2^2

shrink=1
if (shrinkage) {
  shrink = 1-samplingVar/(betaVar*(1-nullfraction)+samplingVar)
#  betaVar = betaVar*shrink^2
#  betaVar2 = betaVar2*shrink^2
#  samplingVar = samplingVar*shrink^2
}

# threshold on betahat based on its p-value
betaHatThreshLo = -qnorm(plower/2)*sqrt(samplingVar)
betaHatThreshHi = -qnorm(pupper/2)*sqrt(samplingVar)

# expected number of selected SNPs
betaHatSD = sqrt(betaVar+samplingVar)
probTruncBeta = 2*nsnp*(1-nullfraction)*abs(pnorm(-betaHatThreshHi,sd=betaHatSD)-
                                            pnorm(-betaHatThreshLo,sd=betaHatSD))
nullHatSD = sqrt(samplingVar)
probTruncNull = 2*nsnp*nullfraction*abs(pnorm(-betaHatThreshHi,sd=nullHatSD)-
                                        pnorm(-betaHatThreshLo,sd=nullHatSD))

# variance of the estimated gene score
if (weighted) {
  if (plower==0) term1=0 else term1=betaHatThreshLo/betaHatSD*dnorm(betaHatThreshLo/betaHatSD)
  if (pupper==0) term2=0 else term2=betaHatThreshHi/betaHatSD*dnorm(betaHatThreshHi/betaHatSD)
  varBetaHat = betaHatSD^2*(1+(term1-term2)/(pnorm(betaHatThreshHi/betaHatSD)-pnorm(betaHatThreshLo/betaHatSD)))
  if (plower==0) term1=0 else term1=betaHatThreshLo/nullHatSD*dnorm(betaHatThreshLo/nullHatSD)
  if (pupper==0) term2=0 else term2=betaHatThreshHi/nullHatSD*dnorm(betaHatThreshHi/nullHatSD)
  varNullHat = samplingVar*(1+(term1-term2)/(pnorm(betaHatThreshHi/nullHatSD)-pnorm(betaHatThreshLo/nullHatSD)))
  varGeneScoreHat = varBetaHat*probTruncBeta+varNullHat*probTruncNull
  #browser()
}
else {
  varGeneScoreHat = probTruncBeta+probTruncNull
}

# covariance between Y2 and estimated gene score
if (weighted) {
# coefficient in SNPs with effects
  scoreCovariance = corr*sqrt(betaVar*betaVar2)/(betaVar+samplingVar)
# covariance in SNPs with effects
  scoreCovariance = scoreCovariance*varBetaHat*probTruncBeta
}
else {
  scoreCovariance = 2*sqrt(betaVar2/betaVar)*corr*(1-nullfraction)*nsnp*
    integrate(discordantSign,0,Inf,sqrt(betaVar),betaHatThreshLo,betaHatThreshHi,sqrt(samplingVar),abs.tol=1e-12)$value
}

# Coefficient of determination!
R2 = scoreCovariance^2/varGeneScoreHat/varY2
# Non-centrality parameter!
NCP=n2*R2/(1-R2)
# Power!
power=pchisq(qchisq(1-alpha,1),1,lower=F,ncp=NCP)

thresholdDensity = dnorm(qnorm(prevalence2))/prevalence2
caseMean = thresholdDensity*R2*varY2/liab2obs2^2
caseVariance = R2*varY2/liab2obs2^2*(1-caseMean*(thresholdDensity+qnorm(prevalence2)))
thresholdDensity = dnorm(qnorm(prevalence2))/(1-prevalence2)
controlMean = -thresholdDensity*R2*varY2/liab2obs2^2
controlVariance = R2*varY2/liab2obs2^2*(1+controlMean*(thresholdDensity-qnorm(prevalence2)))

# debugging
#print(c(probTruncBeta,probTruncNull,varGeneScoreHat,scoreCovariance,caseMean,controlMean,caseVariance,controlVariance))
#print(varGeneScoreHat)

# area under ROC curve!
if (binary) {
  if (logrisk) {
    AUC=pnorm(sqrt(R2*(1-prevalence2)^2/sampling2/(1-sampling2)/2))
  }
  else {
    AUC = pnorm((caseMean-controlMean)/sqrt(caseVariance+controlVariance))
  }
  MSE=NULL
}
else {
  AUC = NULL
  MSE = 1+shrink^2*varGeneScoreHat-2*shrink*scoreCovariance
}

# R2 on liability scale for binary traits
if (binary) R2=R2/liab2obs2^2*sampling2*(1-sampling2)

return(list(R2=R2,NCP=NCP,p=pchisq(NCP+1,1,lower=F),power=power,AUC=AUC,MSE=MSE))
}

discordantSign=function(x,xsigma,threshLo,threshHi,asigma) {
  x*dnorm(x,sd=xsigma)*
  (pnorm(threshLo,mean=x,sd=asigma)-pnorm(threshHi,mean=x,sd=asigma)-pnorm(-threshHi,mean=x,sd=asigma)+pnorm(-threshLo,mean=x,sd=asigma))
  
}

#########
# sampleSizeForAUC: what size of training sample would lead to a given AUC
# parameters as above
#########
sampleSizeForAUC=function(AUC,nsnp,vg1=1,vg2=vg1,corr=1,weighted=T,binary=T,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F,maxN=1e10) {
  obj2=function(p,n1) {
    polygenescore(n1=n1,nsnp=nsnp,vg1=vg1,pupper=p,vg2=vg2,corr=corr,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$AUC
  }
  obj1=function(n1) {
    (optimise(obj2,c(0,1),n1,maximum=T)$objective-AUC)^2
  }
  fit=optimise(obj1,c(0,maxN))
  return(list(n=fit$minimum,p=optimise(obj2,c(0,1),fit$minimum,maximum=T)$maximum))
}

#########
# sampleSizeForCorrelation: what size of training sample would lead to a given correlation
# parameters as above
#########
sampleSizeForCorrelation=function(rho,nsnp,vg1=1,vg2=vg1,corr=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F,maxN=1e10) {
  obj2=function(p,n1) {
    sqrt(polygenescore(n1=n1,nsnp=nsnp,vg1=vg1,pupper=p,vg2=vg2,corr=corr,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$R2)
  }
  obj1=function(n1) {
    (optimise(obj2,c(0,1),n1,maximum=T)$objective-rho)^2
  }
  fit=optimise(obj1,c(0,maxN))
  return(list(n=fit$minimum,p=optimise(obj2,c(0,1),fit$minimum,maximum=T)$maximum))
}

#########
# sampleSizeForPower: what size of sample would lead to a given power, assuming training and sample sizes are equal
# parameters as above
#########
sampleSizeForPower=function(power,nsnp,vg1=1,vg2=vg1,corr=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F,maxN=1e4) {
  obj2=function(p,n1) {
    polygenescore(n1=n1,nsnp=nsnp,vg1=vg1,pupper=p,vg2=vg2,corr=corr,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$power
  }
  obj1=function(n1) {
    (optimise(obj2,c(0,1),n1,maximum=T)$objective-power)^2
  }
  fit=optimise(obj1,c(0,maxN))
  return(list(n=fit$minimum,obj=fit$objective,p=optimise(obj2,c(0,1),fit$minimum,maximum=T)$maximum))
}

#########
# estimateVg1FromP: estimate genetic variance explained by marker panel in discovery sample, given p-value for polygenic score test
# parameters as above
#########
estimateVg1FromP=function(p,n1,nsnp,n2=n1,vg2=0,corr=1,plower=0,pupper=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {
  obj1=function(vg1) {
    if (vg2==0) vg2here=vg1
    else vg2here=vg2
    (sqrt(polygenescore(n1=n1,nsnp=nsnp,vg1=vg1,n2=n2,vg2=vg2here,corr=corr,plower=plower,pupper=pupper,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$NCP)-ncp)^2
  }
#  ncp=qt(p/2,399,lower=F)
ncp=qnorm(p/2,lower=F)
  vg=optimise(obj1,c(0,1))$minimum
#  ncp=qt(.025,399,ncp=qt(p/2,399,lower=F))
  ncp=qnorm(.025,mean=qnorm(p/2,lower=F))
  vgLo=optimise(obj1,c(0,1))$minimum
#  ncp=qt(.975,399,ncp=qt(p/2,399,lower=F))
  ncp=qnorm(.975,mean=qnorm(p/2,lower=F))
  vgHi=optimise(obj1,c(0,1))$minimum
  list(vg=vg,vgLo=vgLo,vgHi=vgHi)
}

#########
# estimateVg2FromP: estimate genetic variance explained by marker panel in replication sample, given p-value for polygenic score test
# parameters as above
#########
estimateVg2FromP=function(p,n1,nsnp,vg1=0,n2=n1,corr=1,plower=0,pupper=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {
  obj1=function(vg2) {
    if (vg1==0) vg1here=vg2
    else vg1here=vg1
    (sqrt(polygenescore(n1=n1,nsnp=nsnp,vg1=vg1here,n2=n2,vg2=vg2,corr=corr,plower=plower,pupper=pupper,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$NCP)-ncp)^2
  }
  ncp=qnorm(p/2,lower=F)
  vg=optimise(obj1,c(0,1))$minimum
  ncp=qnorm(.025,mean=qnorm(p/2,lower=F))
  vgLo=optimise(obj1,c(0,1))$minimum
  ncp=qnorm(.975,mean=qnorm(p/2,lower=F))
  vgHi=optimise(obj1,c(0,1))$minimum
  list(vg=vg,vgLo=vgLo,vgHi=vgHi)
}

#########
# estimateCorrFromP: estimate genetic correlation between two traits explained by marker panel, given p-value for polygenic score test
# parameters as above
#########
estimateCorrFromP=function(p,n1,nsnp,vg1=1,n2=n1,vg2=vg1,plower=0,pupper=1,weighted=T,binary=F,prevalence1=0.1,prevalence2=prevalence1,sampling1=prevalence1,sampling2=prevalence2,lambdaS1=NA,lambdaS2=NA,nullfraction=0,shrinkage=F,logrisk=F) {
  obj1=function(corr) {
    (sqrt(polygenescore(n1=n1,nsnp=nsnp,vg1=vg1,n2=n2,vg2=vg2,corr=corr,plower=plower,pupper=pupper,weighted=weighted,binary=binary,prevalence1=prevalence1,prevalence2=prevalence2,sampling1=sampling1,sampling2=sampling2,lambdaS1=lambdaS1,lambdaS2=lambdaS2,nullfraction=nullfraction,shrinkage=shrinkage,logrisk=logrisk)$NCP)-ncp)^2
  }
  ncp=qnorm(p/2,lower=F)
  corr=optimise(obj1,c(0,1))$minimum
  ncp=qnorm(.025,mean=qnorm(p/2,lower=F))
  corrLo=optimise(obj1,c(0,1))$minimum
  ncp=qnorm(.975,mean=qnorm(p/2,lower=F))
  corrHi=optimise(obj1,c(0,1))$minimum
  list(corr=corr,corrLo=corrLo,corrHi=corrHi)
}

