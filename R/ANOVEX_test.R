ANOVEX_test=function(X,kn=trunc(min(sapply(X,length))/10),L=2,alpha=NULL){

  if(is.list(X)==F){
    stop("X must be a list.")
  }

  J=length(X)

  if(J<2){
    stop("X must contain at least 2 samples.")
  }

  if (length(kn) > 1) {
    stop("k must be of length 1.")
  }

  if (kn > min(sapply(X,length))-1 || kn < 1) {
    stop("kn must be greater than 1, and not exceed the smallest sample size.")
  }

  n=trunc(mean(sapply(X,length)))

  if(is.null(alpha)==T && is.null(L)==T){
    stop("L or alpha missing.")
  }

  if(is.null(alpha)==F && is.null(L)==F && length(alpha)!=L){
    stop("L and alpha must have the same length.")
  }

  if(is.null(alpha)==T && is.null(L)==F){
    alpha=1-(1:L)/n
  }

  if(is.null(alpha)==F && is.null(L)==T){
    L=length(alpha)
  }

  if(any(alpha<=0)==T || any(alpha>=1)==T){
    stop("alpha must be strictly between 0 and 1.")
  }

  if(L<2){
    stop("L must be strictly greater than 1.")
  }

  weisquant=function(Y){
    qhat=quantile(Y,1-kn/length(Y),type=1)
    gammahat=mean(log(quantile(Y,1-(0:(kn-1))/length(Y),type=1)/qhat))
    return(qhat*(kn/(length(Y)*(1-alpha)))^gammahat)
  }
  matquant=sapply(X,weisquant)
  mu=mean(log(matquant))
  mul=apply(log(matquant),1,mean)
  Delta2=mean((mul-mu)^2)
  Delta1l=apply((log(matquant)-matrix(mul,L,J))^2,1,mean)
  Delta1=mean(Delta1l)
  statistic=J*var(log(1-alpha))*kn/(mean(log(kn/(n*(1-alpha)))^2))*Delta1/Delta2
  pvalue=1-pchisq(statistic,df=J-1)
  names(statistic)="T"
  output <- list(statistic = statistic,
               p.value = pvalue,
               method = "ANOVEX test",
               data.name = deparse(substitute(X)))
  class(output) <- "htest"
  return(output)
}
