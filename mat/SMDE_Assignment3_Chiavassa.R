---
  title: "SMDE_Assignment3_Mathieu_Chiavassa"
output: html_document
---
  
  
#Creation m with my E[tau]=78 and my w=exp(2^2) :
createM=function(p,E=78,w=exp(2^2)) { (p*E)/w^0.5 }

#Creation fct to calculate mean and var and C2 of lognormal distrib 
MeanVarLognorm=function(p, sigma=2) {
  m=createM(p)
  w=exp(sigma^2)
  mean=m*(w^0.5)
  var=(m^2)*w*(w-1)
  C2=var/mean^2
  return(list(mean=mean, variance=var, coefVariation=C2))
}

#Creation fct to calculate theorical mean and var of Erlang distrib
MeanVarErlang=function(E_tau){
  mean=E_tau
  lambda = 1/E_tau
  var=mean^2
  return(list(mean=mean, variance=var))
}



#Creation fct to generate data using the recurrent relations :
GenerateData=function(p,N) {
  m=createM(p)
  #generation of service times with sdlog=sigma from my assigned data
  X=rlnorm(N, meanlog=log(m), sdlog=2)
  #generation of arrival times with rate=1/E(tau) from my assigned data
  tau=rexp(N, rate=1/78)
  #Initialization of statistics variables
  L=0; W=0; Lq=0; Wq=0;
  t=vector(mode='numeric', length=N)
  t[1]=0
  ts=vector(mode='numeric', length=N)
  ts[1]=0
  teta=vector(mode='numeric', length=N)
  teta[1]=-1/0
  wq=vector(mode='numeric', length=N)
  wq[1]=ts[1]-t[1]
  w=vector(mode='numeric', length=N)
  w[1]=wq[1]+X[1]
  Li=vector(mode='numeric', length=N)
  Li[1]=w[1]
  L=L+Li[1]
  LT=vector(mode='numeric', length=N)
  LT[1]=0
  lq=vector(mode='numeric', length=N)
  lq[1]=wq[1]
  #Recurrent calculation of the statistics
  for(i in 2:N){
    t[i]=t[i-1]+tau[i-1]
    ts[i]=max(teta[i-1],t[i])
    teta[i]=ts[i]+X[i]
    wq[i]=ts[i]-t[i]
    w[i]=wq[i]+X[i]
    Li[i]=w[i]
    L=L+Li[i]
    LT[i]=L/(t[i]-t[1])
    lq[i]=wq[i]
    Lq=Lq+lq[i]
  }
  
  #After N clients :
  W=sum(w)/N
  Wq=sum(wq)/N
  L=L/t[N]
  Lq=Lq/t[N]
  
  return(list(W=W,Wq=Wq,L=L,Lq=Lq,t=t,LT=LT,lq=lq, wq=wq,N=N,tn=t[N]))
}

#Function to calculate the confience intervals for Wq and Lq:
confidence_inter=function(simu) {
  print("Interval and mean for Wq")
  for(i in 1:length(simu)){
    conf_inter=t.test(simu[[i]]$wq)
    print(conf_inter$conf.int)
    print(conf_inter$estimate)
  }
  print("Interval and mean for Lq")
  for(i in 1:length(simu)){
    conf_inter=t.test(simu[[i]]$wq)
    print((conf_inter$conf.int*simu[[i]]$N)/simu[[i]]$tn)
    print((conf_inter$estimate*simu[[i]]$N)/simu[[i]]$tn)
    
  }
}
#Function to analyze the service time on N instances and compare it to the theoretical values + creating histograms :
analysisServTime=function(p,N) {
  m=createM(p)
  X=rlnorm(N,log(m),sdlog=2)
  meanX=mean(X)
  varX=var(X)
  C2=varX/meanX^2
  theoricValues=MeanVarLognorm(p)
  par(mfrow=c(2,1), mar=c(4,4,4,4))
  
  #histogram on all the sample :
  hist(X,main="Service Time Histogram",xlab="Service Time", ylab="Frequency", ylim=c(0,length(X)), breaks=100, col='red')
  
  #histogram zoom on the revelant part :
  hist(X,main="Zoom",xlab="Service Time", ylab="Frequency",xlim=c(0,300), ylim=c(0,length(X)), breaks=200, col='red')
  
  return(list(meanNum=meanX,varNum=varX, coefVar=C2,meanTheoric=theoricValues$mean,varTheoric=theoricValues$variance, coefVarTheoric=theoricValues$coefVariation))
}

#Allen Cuneen simplified formula with E[tau]=78 and sigma=2 from my assigned data (simplification can be found in the report)
Allen_Cuneen=function(p,E=78,w=exp(2^2)){
  Approximation=E*(p^2/(1-p))*(w/2)
  return(Approximation)
}

#Main function calling almost all the previous ones. You just need to run it to obtain 10 simulations, their graphs, the approximated value of Wq and W with Allen Cuneen formula, the confidence intervals, and others revelant statistics in the "return" :
Main=function(p,N){
  #Generation of 10 simulations
  simu=list(GenerateData(p,N))
  simuWq=vector(mode='numeric', length=10)
  simuWq[1]=simu[[1]]$Wq
  simuLq=vector(mode='numeric', length=10)
  simuLq[1]=simu[[1]]$Lq
  simuW=vector(mode='numeric', length=10)
  simuW[1]=simu[[1]]$W
  simuL=vector(mode='numeric', length=10)
  simuL[1]=simu[[1]]$L
  for(i in 2:10){
    simu=c(simu, list(GenerateData(p,N)))
    simuWq[i]=simu[[i]]$Wq
    simuLq[i]=simu[[i]]$Lq
    simuW[i]=simu[[i]]$W
    simuL[i]=simu[[i]]$L
  }
  
  #plot the LT parameter versus t for the 10 generations 
  file_name=sprintf("./LT_vs_t.jpeg")
  jpeg(filename=file_name)
  par(mfrow=c(5,2), mar=c(1,1,1,1))
  
  for(i in 1:10){
    plot(simu[[i]]$t,simu[[i]]$LT,type='l',xlab='t',ylab='LT')
  }
  dev.off()
  
  #Calculate theorical Allen Cuneen values:
  Approx=Allen_Cuneen(p)
  
  #Calcul of theorical W,Lq,L values with Allen Cuneen approximation
  TheoricW=Approx+p*78
  
  #Confidence interval of the 10 simulations
  confidence_inter(simu)
  
  return(list(Allen_Cuneen=Approx,TheoreticalW=TheoricW,simuWq=simuWq,simuLq=simuLq,simuW=simuW,simuL=simuL,meanWq=mean(simuWq),meanLq=mean(simuLq),meanW=mean(simuW),meanL=mean(simuL)))
}









