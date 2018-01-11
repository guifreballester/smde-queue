# Function to analyze the service time on N instances and compare it to the theoretical values 
# + creating histograms :
analysisServTime=function(rho) {
  a = 0.6521
  E_tau = 66
  
  N=10000
  
  b = create_b(rho)
  
    
  samples = rweibull(N, shape = a, scale = b )
  
  mean_service_samples = mean(samples)
  var_service_samples = var(samples)
  coeff_var_service_samples = sqrt(var_service_samples)/mean_service_samples
  
  mean_service_theoretical = getMeanVarianceWeibull(rho)$mean
  var_service_theoretical = getMeanVarianceWeibull(rho)$var
  coeff_var_service_theoretical = getMeanVarianceWeibull(rho)$coeff_var
  
  par(mfrow=c(2,1), mar=c(4,4,4,4))
  
  #histogram on all the sample :
  hist(samples,main="Service Time Histogram",xlab="Service Time", 
       ylab="Frequency", breaks = 100, col='red')
  
  #histogram zoom on the revelant part :
  hist(samples,main="Zoom",xlab="Service Time", ylab="Frequency",
       xlim=c(0,200), breaks=200, col='red')
  
  return(list(mean_service_samples=mean_service_samples, 
              mean_service_theoretical=mean_service_theoretical, 
              var_service_samples=var_service_samples, 
              var_service_theoretical=var_service_theoretical,
              coeff_var_service_samples=coeff_var_service_samples,
              coeff_var_service_theoretical=coeff_var_service_theoretical))
}

#Creation b with my E[tau]=66:
create_b=function(p,E_tau=66,a=0.6521) {
  return (p * E_tau / gamma((a+1)/a))
}

#Creation fct to calculate theorical mean and var of Erlang distrib
getMeanVarianceErlang=function(k=3, E_tau = 66){
  
  mean=E_tau
  lambda = k/E_tau
  var=k/(lambda)^2
  return(list(mean=mean, variance=var))
}


getMeanVarianceWeibull <- function(p, a=0.6521) {
  b = create_b(p)
  
  mean = b * gamma(1+1/a)
  variance = b^2 * (gamma(1+2/a) - (gamma(1+1/a))^2 )
  
  coeff_var = sqrt(variance)/mean
  
  return(list=list(mean = mean, var = variance, coeff_var = coeff_var))
}

# Return Wq approximation using Allen Cuneen's formula
allen_cunnen_approx <- function(p) {
  #Obtain mean and variance for arrival times and service times for better accuracy
  #Using above given formulas
  E_tau <- getMeanVarianceErlang()$mean
  E_x <- getMeanVarianceWeibull(p)$mean
  var_tau <- getMeanVarianceErlang()$var
  var_x <- getMeanVarianceWeibull(p)$var
  
  lambda <- 1/E_tau
  mu <- 1/E_x
  omega <-lambda/mu
  C <- omega/(1-p+omega)
  
  Wq_approx <- C*(lambda^2*var_tau + mu^2*var_x)/(2*mu*(1-p))
  
  Lq_approx <- (lambda^2*var_tau+p^2*mu^2*var_x)/(2*(1-p))
  return(list(Wq_approx=Wq_approx, Lq_approx=Lq_approx))
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


#Creation fct to generate data using the recurrent relations :
GenerateData=function(p,N=100000) {
  a = 0.6521
  b = create_b(p)
  
  #generation of service times
  X = rweibull(N, shape = a, scale = b )
  
  E_tau = 66
  E_stage = 22
  K = 3 #shape
  
  
  #generation of arrival times with rate=1/E(stage)
  tau = rgamma(N, shape = K, rate = 1/E_stage)
  
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

main=function(p,N){
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
  Approx=allen_cunnen_approx(p)
  
  #Calcul of theorical W,Lq,L values with Allen Cuneen approximation
  #TheoricW=Approx+p*78
  
  #Confidence interval of the 10 simulations
  confidence_inter(simu)
  
  return(list(Allen_Cuneen=Approx,
              simuWq=simuWq,simuLq=simuLq,simuW=simuW,
              simuL=simuL,meanWq=mean(simuWq),meanLq=mean(simuLq),
              meanW=mean(simuW),meanL=mean(simuL)))
}



############1###########
analysisServTime(rho=0.4)



############2###########
p=c(0.4, 0.7, 0.85, 0.925)
N = 100000
  
main(p[4],N)
  
  
  
  