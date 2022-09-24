#2 parallel mediators

require(MASS)

#--- IMPORT USER-SPECIFIED VALUES --------------------------------------------#

powReps <- input$replicationInput
mcmcReps <- input$mcdrawInput
seed <- input$seedInput
conf <- input$ciInput
samplesize <- input$samplesizeInput
input_method <- input$input_method
numIterations <- input$numIterations

#check for user-input value errors:
# CHECK: Is the number of replications > 5 and an integer?
if (powReps < 5 | !abs(powReps - round(powReps)) < .Machine$double.eps ^ 0.5) {
  stop("\"# of Replications\" must be an integer greater than 5. Please change this value.")
}

# CHECK: Is the number of MC replications > 5 and an integer?
if (mcmcReps < 5 | !abs(mcmcReps - round(mcmcReps)) < .Machine$double.eps ^ 0.5) {
  stop("\"Monte Carlo Draws per Rep\" must be an integer greater than 5. Please change this value.")
}

# CHECK: Is the seed > 5 and an integer?
if (seed < 5 | !abs(seed - round(seed)) < .Machine$double.eps ^ 0.5) {
  stop("\"Seed\" must be an integer greater than 5. Please change this value.")
}

# CHECK: Is the confidence level (%) between 0 and 100?
if (conf < 0 | conf > 100) {
  stop("\"Confidence Level (%)\" must be a number between 0 and 100. Please change this value.")
}

calc_2mediator_power <- function(powReps = 10, mcmcReps = 10, seed = 2, 
                                 conf=95, N =100,
                                 input_type = c("sc", "corr"),
                                 stdcoef=matrix(rep(0,ntimes=6)),
                                 corr = matrix(rep(0,ntimes=6)),
                                 SD=matrix(rep(0,ntimes=5))){
  
  if(input_type=="sc"){
    a1 <- as.numeric(stdcoef[1])
    a2 <- as.numeric(stdcoef[2])
    b1 <- as.numeric(stdcoef[3])
    b2 <- as.numeric(stdcoef[4])
    cprime <- as.numeric(stdcoef[5])
    #variables defined here
    #https://github.com/schoam4/mc_power_med/blob/master/code/three_parallel_mediators_stdcoef_ui.R
    #core32 = rM1M2 in GUI for package
    cor32 <- as.numeric(stdcoef[6])
    
    
    if(abs(a1)> .999 | abs(a2)> .999  |
       abs(b1)> .999 | abs(b2)> .999 |
       abs(cprime)> .999 | abs(cor32) > .999
    ) {
      stop("One or more standardized coefficients are out of range (greater than 1 or less than -1)
         check your inputs and try again")
    }
    
    corMat <- diag(4)
    corMat[2,1] <- corMat[1,2] <- a1
    corMat[3,1] <- corMat[1,3] <- a2
    corMat[2,3] <- corMat[3,2] <- cor32
    
    corMat[4,1] <- cprime + a1*b1 + a2*b2
    corMat[1,4] <- cprime + a1*b1 + a2*b2
    corMat[2,4] <- a1*cprime + b1 + b2*cor32
    corMat[4,2] <- a1*cprime + b1 + b2*cor32
    corMat[3,4] <- a2*cprime + b2 + b1*cor32
    corMat[4,3] <- a2*cprime + b2 + b1*cor32
  }else{
    cor21 <- as.numeric(corr[1])
    cor31 <- as.numeric(corr[2])
    cor32 <- as.numeric(corr[3])
    cor41 <- as.numeric(corr[4])
    cor42 <- as.numeric(corr[5])
    cor43 <- as.numeric(corr[6])
    cor51 <- as.numeric(corr[7])
    cor52 <- as.numeric(corr[8])
    cor53 <- as.numeric(corr[9])
    cor54 <- as.numeric(corr[10])
    
    
    # Create correlation / covariance matrix
    corMat <- diag(5)
    corMat[2,1] <- corMat[1,2] <- cor21
    corMat[3,1] <- corMat[1,3] <- cor31
    corMat[2,3] <- corMat[3,2] <- cor32
    corMat[4,1] <- corMat[1,4] <- cor41
    corMat[2,4] <- corMat[4,2] <- cor42
    corMat[3,4] <- corMat[4,3] <- cor43
    corMat[5,1] <- corMat[1,5] <- cor51
    corMat[5,2] <- corMat[2,5] <- cor52
    corMat[5,3] <- corMat[3,5] <- cor53
    corMat[5,4] <- corMat[4,5] <- cor54
    
  }
  
  
  SDX <- as.numeric(SD[1])
  SDM1 <- as.numeric(SD[2])
  SDM2 <- as.numeric(SD[3])
  SDY <- as.numeric(SD[4])
  # Get diagonal matrix of SDs
  SDs <- diag(c(SDX, SDM1, SDM2, SDY))
  
  
  # Convert to covariance matrix
  covMat <- SDs%*%corMat%*%SDs
  
  
  
  
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
  set.seed(seed)
  pow <- lapply(sample(1:50000, powReps), powRep)
  
  #Turn results into a matrix and then compute power for each effect
  df <- data.frame("Parameter" = c("a1b1", "a2b2", "difference12"),
                   "N" = rep(N, 3),
                   "Power" = colSums(matrix(unlist(pow), nrow = powReps, byrow = TRUE)) / powReps)
  
}


#convert 
if (input_method == "Correlations") {
  # Import model input values
  cor21 <- as.numeric(input$cor21)
  cor31 <- as.numeric(input$cor31)
  cor32 <- as.numeric(input$cor32)
  cor41 <- as.numeric(input$cor41)
  cor42 <- as.numeric(input$cor42)
  cor43 <- as.numeric(input$cor43)
  cor51 <- as.numeric(input$cor51)
  cor52 <- as.numeric(input$cor52)
  cor53 <- as.numeric(input$cor53)
  cor54 <- as.numeric(input$cor54)
  
  
  if(any(abs(corMat) > 1)) {
    stop("One or more correlations are out of range (greater than 1 or less than -1)
         check your inputs and try again")
  }
} else {
  a1_low <- as.numeric(input$STa1_low)
  a1_high <- as.numeric(input$STa1_high)
  a2_low <- as.numeric(input$STa2_low)
  a2_high <- as.numeric(input$STa2_high)
  b1_low <- as.numeric(input$STb1_low)
  b1_high <- as.numeric(input$STb1_high)
  b2_low <- as.numeric(input$STb2_low)
  b2_high <- as.numeric(input$STb2_high)
  cor32 <- as.numeric(input$cor32)
  cprime_low <- as.numeric(input$STc_low)
  cprime_high <- as.numeric(input$STc_high)
  
  #create data from inputs
  input_data <- data.frame(a1=runif(numIterations,min= a1_low,max = a1_high),
                           a2=runif(numIterations,min= a2_low,max = a2_high),
                           b1=runif(numIterations,min= b1_low,max = b1_high),
                           b2=runif(numIterations,min= b2_low,max = b2_high),
                           cprime=runif(numIterations,min=cprime_low,max=cprime_high),
                           rm1m2=rep(cor32,numIterations))
}


power_estimates_2mediators <- data.frame(a1b1 = rep(0,numIterations), a2b2 = rep(0,numIterations),
                                         difference12 =rep(0,numIterations))

for(i in 1:numIterations){
#  withProgress(message = 'Running Replications', value = 0, {
    #define standard deviation
    SD <- c(1,1,1,1)
    #define standard coefficients
    SC <- input_data[i,]
    power_estimates_2mediators[i,] <- calc_2mediator_power(powReps=powReps, mcmcReps = mcmcReps, seed = seed,
                                                           conf = conf, N = samplesize, input_type=input_method, stdcoef = SC, SD = SD)[,3]
    print(i)
  #})
}

power_means <- data.frame(a1b1_mean = mean(power_estimates_2mediators[,1]),
                          a2b2_mean = mean(power_estimates_2mediators[,2]),
                          difference12 = mean(power_estimates_2mediators[,3]))
power_means












