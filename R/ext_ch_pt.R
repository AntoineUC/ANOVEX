ext_ch_pt=function(Y,X=1:length(Y),kn=trunc(length(Y)/20),method="ANOVEX",grid=seq(quantile(X,0.15,type=1),quantile(X,0.85,type=1),length=sum((as.numeric(names(table(X)))<=quantile(X,0.85))*(as.numeric(names(table(X)))>=quantile(X,0.15)))),plot=FALSE){

  n=length(Y)

  if (n != length(X)) {
    stop("X and Y must have the same length !")
  }

  if (length(kn) > 1) {
    stop("kn must be of length 1.")
  }

  if (max(grid) > max(X) || min(grid) < min(X)) {
    stop("grid must be between min(X) and max(X).")
  }

  if (kn > n - 1 || kn < 1) {
    stop("kn must be between 1 and n-1.")
  }

  if(method!="ANOVEX" && method!="KL" && method!="LR-GPD" && method!="LR-Pareto" && method!="Quantile"){
    stop("method must be ANOVEX, KL, LR-GPD, LR-Pareto or Quantile.")
  }

  if(!is.logical(plot)){
    stop("plot must be boolean.")
  }

  if(method=="ANOVEX"){

    testat=function(s){

      Y1=Y[which(X<=s)]
      Y2=Y[which(X>s)]

      n1=length(Y1)
      n2=length(Y2)

      L=min(trunc(kn*n/n1),trunc(kn*n/n2))
      alphan=1-(1:L)/n

      q1=quantile(Y1,1-kn/n1,type=1)
      q2=quantile(Y2,1-kn/n2,type=1)

      xi1=mean(log(quantile(Y1,1-(0:(kn-1))/n1,type=1)/quantile(Y1,1-kn/n1,type=1)))
      xi2=mean(log(quantile(Y2,1-(0:(kn-1))/n2,type=1)/quantile(Y2,1-kn/n2,type=1)))

      q1ex=rep(q1,length(alphan))*(kn/(n1*(1-alphan)))^rep(xi1,length(alphan))
      q2ex=rep(q2,length(alphan))*(kn/(n2*(1-alphan)))^rep(xi2,length(alphan))

      mu_alpha_n=(sum(log(q1ex))+sum(log(q2ex)))/(2*L)
      mu_alpha_ln=(log(q1ex)+log(q2ex))/2
      Delta_1l=((log(q1ex)-mu_alpha_ln)^2+(log(q2ex)-mu_alpha_ln)^2)/2
      Delta_1=sum(Delta_1l)/L
      Delta_2=sum((mu_alpha_ln-mu_alpha_n)^2)/L
      return(Delta_1/Delta_2*2*kn*(L-1)*var(log(1:L))/sum(log(kn/(1:L))^2))
    }
    title="ANOVEX"
  }

  if(method=="KL"){

    testat=function(s){

      u=quantile(Y,1-kn/n,type=1)
      Mt=sum(log((Y[which(X<=s)])[which(Y[which(X<=s)]>u)]/u))
      MT=sum(log(Y[which(Y>u)]/u))

      return(abs(Mt-s*MT/max(X)))
    }
    title="KL"
  }

  if(method=="Quantile"){

    testat=function(s){

      u=quantile(Y,1-kn/n,type=1)
      St=sum((1-kn/n)-(Y[which(X<=s)]<=u))
      ST=sum((1-kn/n)-(Y<=u))

      return(abs(St-s*ST))
    }
    title="Quantile"
  }

  if(method=="LR-Pareto"){

    testat=function(s){

      Y1=Y[which(X<=s)]
      Y2=Y[which(X>s)]
      n1=length(Y1)
      n2=length(Y2)
      k1=trunc((sum(X<=s)/n)*kn)
      Hill_1=mean(log(quantile(Y1,1-(0:(k1-1))/n1,type=1)/quantile(Y1,1-k1/n1,type=1)))
      Hill_2=mean(log(quantile(Y2,1-(0:(kn-k1-1))/n2,type=1)/quantile(Y2,1-(kn-k1)/n2,type=1)))
      Hill=mean(log(quantile(Y,1-(0:(kn-1))/n,type=1)/quantile(Y,1-kn/n,type=1)))

      return(-2*(k1*log(Hill_1)+(kn-k1)*log(Hill_2)-kn*log(Hill))-(k1*Hill_1+(kn-k1)*Hill_2-kn*Hill)/Hill)
    }
    title="LR-Pareto"
  }

  if(method=="LR-GPD"){

    testat=function(s){

      k1=trunc((sum(X<=s)/n)*kn)
      m=sum(X<=s)

      Y1=Y[which(X<=s)]
      Y2=Y[which(X>s)]

      u=quantile(Y,1-kn/n,type=1)
      u1=quantile(Y1,1-k1/m,type=1)
      u2=quantile(Y2,1-(kn-k1)/(n-m),type=1)

      gpd_est1=mev::fit.gpd(Y1[which(Y1>u1)]-u1,threshold = 0)$estimate
      gpd_est2=mev::fit.gpd(Y2[which(Y2>u2)]-u2,threshold = 0)$estimate
      gpd_est=mev::fit.gpd(Y[which(Y>u)]-u,threshold = 0)$estimate

      xi1=gpd_est1[2]
      sigma1=gpd_est1[1]
      xi2=gpd_est2[2]
      sigma2=gpd_est2[1]
      xi_all=gpd_est[2]
      sigma_all=gpd_est[1]

      L1=-k1*log(sigma1)-(1/xi1+1)*sum(log(1+xi1*(Y1[which(Y1>u1)]-u1)/sigma1))
      L2=-(kn-k1)*log(sigma2)-(1/xi2+1)*sum(log(1+xi2*(Y2[which(Y2>u2)]-u2)/sigma2))
      L3=-kn*log(sigma_all)-(1/xi_all+1)*sum(log(1+xi_all*(Y[which(Y>u)]-u)/sigma_all))

      return(2*(L1+L2-L3))

    }
    title="LR-GPD"
  }

  vec_testat=mapply(testat,grid)

  shat=grid[which.max(vec_testat)]

  if(plot==TRUE){
    plot(y=Y,x=X,pch=16,ylab="Y",xlab="X",main=paste(title,"split"))
    par(new=T)
    plot(y=vec_testat,x=grid,type="l",lwd=2,col="blue",xlab="",ylab="",axes=F,xlim=c(min(X),max(X)))
    abline(v=shat,lwd=2,col="red",lty=5)
  }

  return(shat)

}
