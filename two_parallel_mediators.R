calc_2mediator_power <- function(powReps = 10, mcmcReps = 10, seed = 2, conf=95, N =100,
                                 stdcoef=matrix(rep(0,ntimes=6)),SD=matrix(rep(0,ntimes=4))){
  
  
  #create coefficient matrix
  a1 <- as.numeric(stdcoef[1])
  a2 <- as.numeric(stdcoef[2])
  b1 <- as.numeric(stdcoef[3])
  b2 <- as.numeric(stdcoef[4])
  cprime <- as.numeric(stdcoef[5])
  #COR32=r FOUND HERE: https://github.com/schoam4/mc_power_med/blob/master/code/two_parallel_mediators_stdcoef_ui.R
  cor32 <- as.numeric(stdcoef[6])
  
  
  #stdcoef = a1, a2, b1, b2, cprime, cor32
  #stdcoef = c(0.14,0.24,0.3,0.33, 0.15,.14)
  #combine coefficients into matrix
  corMat <- diag(4)
  corMat[2,1] <- a1
  corMat[1,2] <- a1
  corMat[3,1] <- a2
  corMat[1,3] <- a2
  corMat[2,3] <- cor32
  corMat[3,2] <- cor32
  corMat[4,1] <- cprime + a1*b1 + a2*b2
  corMat[1,4] <- cprime + a1*b1 + a2*b2
  corMat[2,4] <- a1*cprime + b1 + b2*cor32
  corMat[4,2] <- a1*cprime + b1 + b2*cor32
  corMat[3,4] <- a2*cprime + b2 + b1*cor32
  corMat[4,3] <- a2*cprime + b2 + b1*cor32
  
  #standard deviations for each variable
  #SD <- x, m1,m2,y
  #SD <- c(1,1,1,1)
  SDX <- as.numeric(SD[1])
  SDM1 <- as.numeric(SD[2])
  SDM2 <- as.numeric(SD[3])
  SDY <- as.numeric(SD[4])
  # Get diagonal matrix of SDs
  SDs <- diag(c(SDX, SDM1, SDM2, SDY))
  
  # Convert to covariance matrix
  covMat <- SDs%*%corMat%*%SDs
  
  
  # Start power simulation
  # Create function for 1 rep
  seed = 1234
  powRep <- function(seed = 1234, Ns = N, covMatp = covMat){
    set.seed(seed)
    require(MASS)
    
    dat <- mvrnorm(Ns, mu = c(0,0,0,0), Sigma = covMatp)
    # Run regressions
    m1 <- lm(dat[,2] ~ dat[,1])
    m2 <- lm(dat[,3] ~ dat[,1])
    m3 <- lm(dat[,4] ~ dat[,2] + dat[,3] + dat[,1])
    
    # Output parameter estimates and standard errors
    a1 <- rnorm(mcmcReps, coef(m1)[2], sqrt(vcov(m1)[2,2]))
    a2 <- rnorm(mcmcReps, coef(m2)[2], sqrt(vcov(m2)[2,2]))
    b1 <- rnorm(mcmcReps, coef(m3)[2], sqrt(vcov(m3)[2,2]))
    b2 <- rnorm(mcmcReps, coef(m3)[3], sqrt(vcov(m3)[3,3]))
    
    a1b1 <- a1*b1
    a2b2 <- a2*b2
    diff <- a1*b1 - a2*b2
    
    # Calculate confidence intervals
    low <- (1 - (conf / 100)) / 2
    upp <- ((1 - conf / 100) / 2) + (conf / 100)
    LL1 <- quantile(a1b1, low)
    UL1 <- quantile(a1b1, upp)
    LL2 <- quantile(a2b2, low)
    UL2 <- quantile(a2b2, upp)        
    LLd <- quantile(diff, low)
    ULd <- quantile(diff, upp)  
    
    # Is rep significant?
    c(LL1*UL1 > 0, LL2*UL2 > 0, LLd*ULd > 0)
    
  }
  set.seed(1234)
  pow <- lapply(sample(1:50000, powReps), powRep)
  
  #Turn results into a matrix and then compute power for each effect
  df <- data.frame("Parameter" = c("a1b1", "a2b2", "difference"),
                   "N" = rep(N, 3),
                   "Power" = colSums(matrix(unlist(pow), nrow = powReps, byrow = TRUE)) / powReps)
  
  return(df)
}