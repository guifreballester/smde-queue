# Individual parameters from .pdf
std <- 1.6333 # (Service)
a <- 2; b = 72 # Weibull distribution for interarrival times. b will vary for service time
loadingFactors <- c(0.4, 0.7, 0.85, 0.925) # loading factors
N <- 100000 # number of clients

# Function to obtain mean and variance of lognormal distribution,
# given free parameter m and sigma
getMeanVarianceLognor <- function(m, sigma = std) {
  omega <- exp(sigma^2)
  mean <- m*sqrt(omega)
  variance <- m^2*omega*(omega-1)
  
  return(list=list(mean = mean, var = variance))
}

# Since a = 2, we can simplify values of gamma functions
# Function to obtain mean and variance of Weibull distribution
# Default parameters are given in .pdf file
getMeanVarianceWeibull <- function(aa = a, bb = b) {
  mean <- bb*sqrt(pi)/2
  variance <- bb^2*(1 - pi/4)
  
  return(list=list(mean = mean, var = variance))
}

# Function to generate N clients of Waiting System
# Interarrival times follow Weibull distribution
# Service times follow Lognormal distribution
# Function will return N arrival, service, and entrance times,
# LT, sojourn times, and sojourn times in queue. And print statistics:
# Average occupancy in W.S. and Queue, and Average waiting time in W.S. and Queue

generateData <- function(N, a, b, m, sigma = std) {
  tau <- vector(mode = 'numeric', length = N) # interarrival times for N clients
  x <- vector(mode = 'numeric', length = N) # service times for N clients 
  
  t <- vector(mode='numeric',length = N) # entrance time instant for the W.S. 
  omega = vector(mode = 'numeric', length = N) # exit time instant from W.S. for N clients
  omega[1] = -1/0 # infinity
  
  ts = vector(mode = 'numeric', length = N) # arrival time instant to the service system
  L = vector(mode = 'numeric', length = N) # contribution of i-th client to occupancy 
  LT = vector(mode = 'numeric', length = N) 
  w = vector(mode = 'numeric', length = N) # sojourn time of client i in W.S.
  wq = vector(mode = 'numeric', length = N) # sojourn time of client i in queue
  lq = vector(mode = 'numeric', length = N)
  
  # Generating data using recurrent relations
  LL = 0; W = 0; Lq = 0; Wq = 0; # Average occupancy, sojourn time, occupancy in queue, sojourn time in queue
  for(i in 1:N) {
    tau[i] <- rweibull(1, a, b)
    x[i] <- rlnorm(1, meanlog = log(m), sdlog = sigma)
    if(i > 1) {
      ts[i] = max(omega[i-1], t[i]) 
    } else {
      ts[i] = t[i]
    }
    omega[i] = ts[i] + x[i]
    t[i+1] = t[i] + tau[i]
    
    L[i] = omega[i] - t[i]
    w[i] = L[i]
    LL = LL + L[i]
    LT[i] = LL/(t[i] - t[1])
    W = W + w[i]
    
    lq[i] = ts[i] - t[i]
    wq[i] = lq[i]
    Lq = Lq + lq[i]
    Wq = Wq + wq[i]
  }
  
  # After N clients report, calculate final printing statistics
  W = W/N
  Wq = Wq/N
  LL = LL/(t[N] - t[1])
  Lq = Lq/(t[N] - t[1])
  t = t[1:N]
  
  # Return all necessary values as a list
  return(list(arrivalTime=tau,serviceTime=x,entranceTime=t, LT=LT, sojTime=w, sojTimeQ=wq,
              occQ = lq, avOcc=LL,avSojTime=W,avOccQ=Lq,avSojTimeQ=Wq,length=N))
}

# Function to find m of Lognormal distribution. expArrTime = E[tau]=64, sigma=1.6333 from pdf

findM <- function(loadFactor, expArrTime = 64, sigma=std) {
  omega = exp(sigma*sigma)
  m = loadFactor*expArrTime/sqrt(omega)
  return(m)
}

# Function that performs 10 simulations in parallel, plots them in big figure,
# and later saves in the same directory where file was called
# Function also returns list of 10 simulations, that contain output of function generateData(...)
data <- NULL
library(parallel)
plot_LT <- function(rpos, dat = NULL) {
  # Generate data using parallel processing
  if(is.null(dat)) {
    nCores <- detectCores() - 1
    m <<- mclapply(1:4, function(i) findM(loadingFactors[i]), mc.silent = TRUE, mc.cores = nCores)
    data <<- mclapply(1:10, function(i) generateData(N,a,b,m[[rpos]]), mc.silent = TRUE, mc.cores = nCores)
  } else {
    data <<- dat
  }
  
  # Saving plot to file, without plotting
  k <- 0
  while(TRUE) {
    file_name <- sprintf("./LT_t_%d_%d.png",rpos,k)
    if(!file.exists(file_name)) {
      break
    }
    k <- k + 1
  }
  png(filename = file_name, width = 1900, height = 1200)
  par(mfrow=c(5,2))
  for(i in 1:10) {
    d <- data[[i]]
    title <- sprintf("p (Loading Factor) = %.3f",loadingFactors[rpos])
    plot(d$entranceTime, d$LT, type='l', xlab = 't', ylab = 'LT', main = title)
    abline(d$avOcc, 0 , col = 'red')
  }
  dev.off()
  return(list(data=data, M=m))
}

# Save image of histograms of Service Time distributions in the fun
plot_histograms <- function(rpos, data) {
  #file_name <- sprintf("~/Documents/tex/SMDE_Lab3/images/histogram_%d.png",rpos)
  k <- 0
  while(TRUE) {
    file_name <- sprintf("./histogram_%d_%d.png",rpos, k)
    if(!file.exists(file_name)) {
      break
    }
    k <- k + 1
  }
  
  png(filename = file_name, width = 1900, height = 1200)
  par(mfrow=c(5,2))
  for(i in 1:10) {
    d <- data[[i]]
    title <- sprintf("Service Times Histrogram for p (Loading Factor) = %.3f",loadingFactors[rpos])
    hist(d$serviceTime, xlab = "Service Time",ylab = "Frequency", main = title, breaks = 200, col = 'blue')
  }
  dev.off()
}

# Return Wq approximation using Allen Cuneen's formula
allen_cunnen_approx <- function(p, m) {
  #Obtain mean and variance for arrival times and service times for better accuracy
  #Using above given formulas
  E_x <- getMeanVarianceLognor(m)$mean
  E_tau <- getMeanVarianceWeibull()$mean
  var_tau <- getMeanVarianceWeibull()$var
  var_x <- getMeanVarianceLognor(m)$var
  
  lambda <- 1/E_tau
  mu <- 1/E_x
  omega <-lambda/mu
  C <- omega/(1-p+omega)
  
  approx <- C*(lambda^2*var_tau + mu^2*var_x)/(2*mu*(1-p))
  return(approx)
}

# Print confidence 95% intervals using t-statistic

print_confidence_intervals <- function(data) {
  for(i in 1:length(data)) {
    conf_interval <- t.test(data[[i]]$sojTimeQ)
    print(conf_interval)
  }
}