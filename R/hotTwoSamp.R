hotTwosamp<-function(x1,x2,delta){
  n1<-nrow(x1)
  n2<-nrow(x2)
  p<-ncol(x1)
  x1bar<-as.matrix(colMeans(x1),nrow=2)
  x2bar<-as.matrix(colMeans(x2),nrow=2)
  spooled<-((n1-1)*cov(x1)+(n2-1)*cov(x2))/(n1+n2-2)
  sp<-((1/n1)+(1/n2))*spooled
  T2<-t((x1bar-x2bar)-delta)%*%solve(sp)%*%((x1bar-x2bar)-delta)
  F_stat<-(n1+n2-p-1)*T2/((n1+n2-2)*p)
  pvalue<-pf(F_stat,p,n1+n2-p-1,lower.tail = F)
  return(list(T2=T2,Fs=F_stat,p=pvalue))
}
#' hotelling two sample T^2 test
#' 
#' takes multivariate data, with two independent vector. 
#' @param x A matrix
#' @param y A matrix
#' @param delta Difference
#' @returns T^2 statistic, F statistic and the p-value.

