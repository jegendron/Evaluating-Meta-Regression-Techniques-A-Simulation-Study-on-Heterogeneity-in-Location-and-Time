
### To allow audio feedback
library(beepr)

##############################################################################
###                     DATA SETUP - (Find True Parameters)                ###
##############################################################################

### Start by setting up variance covariance matrix (presumed the same for all studies)
library(mvtnorm)
library(Matrix)

# matrix syntax: (TL, BL, TR, BR)
S00 <- matrix(c(1, 0.5, 0.5, 1), nrow = 2, ncol = 2)      # top left
S01 <- matrix(c(0.25, 0, 0.2, 0.1), nrow = 2, ncol = 2)   # top right
S10 <- t(S01)                                             # bottom left
S11 <- matrix(c(1, 0.1, 0.1, 1), nrow = 2, ncol = 2)      # bottom right

SIGMA <- rbind(cbind(S00, S01), cbind(S10, S11))          # VarCov matrix

# Elements used to calculate the true parameters
varY <- S00[1,1]
covY <- SIGMA[2:length(SIGMA[1,]),1:1]
covX <- SIGMA[2:length(SIGMA[1,]),2:length(SIGMA[1,])]   
  # the length() extracts the length of one row
trueBetas <- solve(covX) %*% t(t(covY)) # t(t()) to transpose once, likely because of datatype

### The true parameters
sig2True <- varY - covY %*% trueBetas                     # Omega
b1True <- trueBetas[1,1]
b2True <- trueBetas[2,1]
b3True <- trueBetas[3,1]

R2True <- 1 - sig2True/varY #Not used due to issues with how R2 is measured for fixed vs random vs mixed effects



##############################################################################
###                       DATA SETUP - (Generate the data)                 ###
##############################################################################

#HOW THE DATA IS GENERATED
  #-Each (country) has 1 study ran per year for 5 years
  #-Heterogeneity in countries is captured by the different mean value for each country
  #-Heterogeneity in time is captured by the trend in y=y+0.5*i

library(tidyverse)

# Means
mu_x1 <- 1
mu_x2 <- 0.5
mu_x3 <- 1.5

numYears <- 5

#####################################
### SPECIFY SIMULATION PARAMETERS ###
#####################################

# Monte Carlo iterations
#
N <- 10000
#N <- 1000
#N <- 100
#N <- 1
  
# Sample size per study
#
n <- 50
#n <- 100
#n <- 150

case <- 1 # Options: 1-12

model1 <- TRUE # RE s
model2 <- TRUE # RE l
model3 <- TRUE # FE s
model4 <- TRUE # FE l
model5 <- TRUE # FE t
model6 <- TRUE # FE lt
model7 <- TRUE # FE s (with Trend)
model8 <- TRUE # FE l (with Trend)

model9 <- TRUE # ME s
model10 <- TRUE# ME l

if(case==1){
  numCountries <-5
  timeHet <-0.1 #small
  spreadMuY <- -2 #small
}else if(case==2){
  numCountries <-5
  timeHet <-0.1 #small
  spreadMuY <- -10 #large
}else if(case==3){
  numCountries <-9
  timeHet <-0.1 #small
  spreadMuY <- -2 #small
}else if(case==4){
  numCountries <-9
  timeHet <-0.1 #small
  spreadMuY <- -10 #large
}else if(case==7){
  numCountries <-5
  timeHet <-0.5 #large
  spreadMuY <- -2 #small
}else if(case==8){
  numCountries <-5
  timeHet <-0.5 #large
  spreadMuY <- -10 #large
}else if(case==9){
  numCountries <-9
  timeHet <-0.5 #large
  spreadMuY <- -2 #small
}else if(case==10){
  numCountries <-9
  timeHet <-0.5 #large
  spreadMuY <- -10 #large
}else if(case==5){
  numCountries <-15
  timeHet <-0.1 #small
  spreadMuY <- -2 #small
}else if(case==6){
  numCountries <-15
  timeHet <-0.1 #small
  spreadMuY <- -10 #large
}else if(case==11){
  numCountries <-15
  timeHet <-0.5 #large
  spreadMuY <- -2 #small
}else if(case==12){
  numCountries <-15
  timeHet <-0.5 #large
  spreadMuY <- -10 #large
}
# timeHet <- Will time heterog be small or large?
# spreadMuY <- Will muY have a large variation (-10 to 10), or small variation (-2 to 2)?



# This is how we account for the Large N problem, following Good (1988) "The Interface between Statistics and Philosophy of Science"
if(numYears==5 & numCountries==5){
  if(n==50){            # N = 1,250
    alphaLevel<-0.001
  }
  else if(n==100){      # N = 2,250
    alphaLevel<-0.0001
  }
  else if(n==150){      # N = 3,750
    alphaLevel<-0.00005   # May need to adjust if overpowered (as long as >= 0.000001)
  }
}
if(numYears==5 & numCountries==9){
  if(n==50){            # N = 2,500
    alphaLevel<-0.0001
  }
  else if(n==100){      # N = 4,500
    alphaLevel<-0.00001   # May need to adjust if overpowered (as long as >= 0.000001)
  }
  else if(n==150){      # N = 7,500
    alphaLevel<-0.000005  # May need to adjust if overpowered (as long as >= 0.000001)
  }
}
if(numYears==5 & numCountries==15){
  if(n==50){            # N = 3,750
    alphaLevel<-0.00005   # May need to adjust if overpowered (as long as >= 0.000001)
  }
  else if(n==100){      # N = 6,750
    alphaLevel<-0.000005  # May need to adjust if overpowered (as long as >= 0.000001)
  }
  else if(n==150){      # N = 11,250
    alphaLevel<-0.000001   
  }
}



##############################################################################
###                                 MONTE CARLO                            ###
##############################################################################

library(tidyverse)
library(lmtest)
library(plm)
library(lme4)
library(sandwich)

# All of these empty vectors will store various results for each model
  # (i.e. REs, FEs, FE_lt, ...)

REb0 <- numeric(N)
FEb0 <- numeric(N)
FE2b0 <- numeric(N)
REl_b0 <- numeric(N)
FEl_b0 <- numeric(N)
FEt_b0 <- numeric(N)
FEsT_b0 <- numeric(N)
FElT_b0 <- numeric(N)
MEb0 <- numeric(N)
MEl_b0 <- numeric(N)

REb1 <- numeric(N) # x1
FEb1 <- numeric(N) # x1
FE2b1 <- numeric(N) # x1
REl_b1 <- numeric(N) # x1
FEl_b1 <- numeric(N) # x1
FEt_b1 <- numeric(N) # x1
FEsT_b1 <- numeric(N) # x1
FElT_b1 <- numeric(N) # x1
MEb1 <- numeric(N)
MEl_b1 <- numeric(N)

REb2 <- numeric(N) # x2
FEb2 <- numeric(N) # x2
FE2b2 <- numeric(N) # x2
REl_b2 <- numeric(N) # x2
FEl_b2 <- numeric(N) # x2
FEt_b2 <- numeric(N) # x2
FEsT_b2 <- numeric(N) # x2
FElT_b2 <- numeric(N) # x2
MEb2 <- numeric(N)
MEl_b2 <- numeric(N)

REb3 <- numeric(N) # x3
FEb3 <- numeric(N) # x3
FE2b3 <- numeric(N) # x3
REl_b3 <- numeric(N) # x3
FEl_b3 <- numeric(N) # x3
FEt_b3 <- numeric(N) # x3
FEsT_b3 <- numeric(N) # x3
FElT_b3 <- numeric(N) # x3
MEb3 <- numeric(N)
MEl_b3 <- numeric(N)

FEsT_bTrend <- numeric(N)
FElT_bTrend <- numeric(N)

REsig2 <- numeric(N)
FEsig2 <- numeric(N)
FE2sig2 <- numeric(N) 
REl_sig2 <- numeric(N)
FEl_sig2 <- numeric(N)
FEt_sig2 <- numeric(N)
FEsT_sig2 <- numeric(N)
FElT_sig2 <- numeric(N) 
MEsig2 <- numeric(N)
MEl_sig2 <- numeric(N)

REb0SE <- numeric(N)
FEb0SE <- numeric(N)
FE2b0SE <- numeric(N)
REl_b0SE <- numeric(N)
FEl_b0SE <- numeric(N)
FEt_b0SE <- numeric(N)
FEsT_b0SE <- numeric(N)
FElT_b0SE <- numeric(N)
MEb0SE <- numeric(N)
MEl_b0SE <- numeric(N)

REb1SE <- numeric(N) # x1
FEb1SE <- numeric(N) # x1
FE2b1SE <- numeric(N) # x1
REl_b1SE <- numeric(N) # x1
FEl_b1SE <- numeric(N) # x1
FEt_b1SE <- numeric(N) # x1
FEsT_b1SE <- numeric(N) # x1
FElT_b1SE <- numeric(N) # x1
MEb1SE <- numeric(N)
MEl_b1SE <- numeric(N)

REb2SE <- numeric(N) # x2
FEb2SE <- numeric(N) # x2
FE2b2SE <- numeric(N) # x2
REl_b2SE <- numeric(N) # x2
FEl_b2SE <- numeric(N) # x2
FEt_b2SE <- numeric(N) # x2
FEsT_b2SE <- numeric(N) # x2
FElT_b2SE <- numeric(N) # x2
MEb2SE <- numeric(N)
MEl_b2SE <- numeric(N)

REb3SE <- numeric(N) # x3
FEb3SE <- numeric(N) # x3
FE2b3SE <- numeric(N) # x3
REl_b3SE <- numeric(N) # x3
FEl_b3SE <- numeric(N) # x3
FEt_b3SE <- numeric(N) # x3
FEsT_b3SE <- numeric(N) # x3
FElT_b3SE <- numeric(N) # x3
MEb3SE <- numeric(N)
MEl_b3SE <- numeric(N)

FEsT_bTrendSE <- numeric(N)
FElT_bTrendSE <- numeric(N)

RE_i2 <-numeric(N)
FE_i2 <- numeric(N)
FE2_i2 <- numeric(N)
REl_i2 <- numeric(N)
FEl_i2 <- numeric(N)
FEt_i2 <- numeric(N)
FEsT_i2 <- numeric(N)
FElT_i2 <- numeric(N)
ME_i2 <- numeric(N)
MEl_i2 <- numeric(N)

RE_h2 <-numeric(N)
FE_h2 <- numeric(N)
FE2_h2 <- numeric(N)
REl_h2 <- numeric(N)
FEl_h2 <- numeric(N)
FEt_h2 <- numeric(N)
FEsT_h2 <- numeric(N)
FElT_h2 <- numeric(N)
ME_h2 <- numeric(N)
MEl_h2 <- numeric(N)

RE_Q <-numeric(N)
FE_Q <- numeric(N)
FE2_Q <- numeric(N)
REl_Q <- numeric(N)
FEl_Q <- numeric(N)
FEt_Q <- numeric(N)
FEsT_Q <- numeric(N)
FElT_Q <- numeric(N)
ME_Q <- numeric(N)
MEl_Q <- numeric(N)

REr2 <- numeric(N)
FEr2 <- numeric(N)
FE2r2 <- numeric(N)
REl_r2 <- numeric(N)
FEl_r2 <- numeric(N)
FEt_r2 <- numeric(N)
FEsT_r2 <- numeric(N)
FElT_r2 <- numeric(N)
ME_r2 <- numeric(N)
MEl_r2 <- numeric(N)

REa_r2 <- numeric(N)
FEa_r2 <- numeric(N)
FE2a_r2 <- numeric(N)
REl_a_r2 <- numeric(N)
FEl_a_r2 <- numeric(N)
FEt_a_r2 <- numeric(N)
FEsT_a_r2 <- numeric(N)
FElT_a_r2 <- numeric(N)
ME_a_r2 <- numeric(N)
MEl_a_r2 <- numeric(N)

REaic <- numeric(N)
FEaic <- numeric(N)
FE2aic <- numeric(N)
REl_aic <- numeric(N)
FEl_aic <- numeric(N)
FEt_aic <- numeric(N)
FEsT_aic <- numeric(N)
FElT_aic <- numeric(N)
ME_aic <- numeric(N)
MEl_aic <- numeric(N)

REmse <- numeric(N)
FEmse <- numeric(N)
FE2mse <- numeric(N)
REl_mse <- numeric(N)
FEl_mse <- numeric(N)
FEt_mse <- numeric(N)
FEsT_mse <- numeric(N)
FElT_mse <- numeric(N)
MEmse <- numeric(N)
MEl_mse <- numeric(N)

REmae <- numeric(N)
FEmae <- numeric(N)
FE2mae <- numeric(N)
REl_mae <- numeric(N)
FEl_mae <- numeric(N)
FEt_mae <- numeric(N)
FEsT_mae <- numeric(N)
FElT_mae <- numeric(N)
MEmae <- numeric(N)
MEl_mae <- numeric(N)

REmpe <- numeric(N)
FEmpe <- numeric(N)
FE2mpe <- numeric(N)
REl_mpe <- numeric(N)
FEl_mpe <- numeric(N)
FEt_mpe <- numeric(N)
FEsT_mpe <- numeric(N)
FElT_mpe <- numeric(N)
MEmpe <- numeric(N)
MEl_mpe <- numeric(N)

REmape <- numeric(N)
FEmape <- numeric(N)
FE2mape <- numeric(N)
REl_mape <- numeric(N)
FEl_mape <- numeric(N)
FEt_mape <- numeric(N)
FEsT_mape <- numeric(N)
FElT_mape <- numeric(N)
MEmape <- numeric(N)
MEl_mape <- numeric(N)

REmse_x1 <- numeric(N)
FEmse_x1 <- numeric(N)
FE2mse_x1 <- numeric(N)
REl_mse_x1 <- numeric(N)
FEl_mse_x1 <- numeric(N)
FEt_mse_x1 <- numeric(N)
FEsT_mse_x1 <- numeric(N)
FElT_mse_x1 <- numeric(N)
ME_mse_x1 <- numeric(N)
MEl_mse_x1 <- numeric(N)

REmae_x1 <- numeric(N)
FEmae_x1 <- numeric(N)
FE2mae_x1 <- numeric(N)
REl_mae_x1 <- numeric(N)
FEl_mae_x1 <- numeric(N)
FEt_mae_x1 <- numeric(N)
FEsT_mae_x1 <- numeric(N)
FElT_mae_x1 <- numeric(N)
ME_mae_x1 <- numeric(N)
MEl_mae_x1 <- numeric(N)

REmpe_x1 <- numeric(N)
FEmpe_x1 <- numeric(N)
FE2mpe_x1 <- numeric(N)
REl_mpe_x1 <- numeric(N)
FEl_mpe_x1 <- numeric(N)
FEt_mpe_x1 <- numeric(N)
FEsT_mpe_x1 <- numeric(N)
FElT_mpe_x1 <- numeric(N)
ME_mpe_x1 <- numeric(N)
MEl_mpe_x1 <- numeric(N)

REmape_x1 <- numeric(N)
FEmape_x1 <- numeric(N)
FE2mape_x1 <- numeric(N)
REl_mape_x1 <- numeric(N)
FEl_mape_x1 <- numeric(N)
FEt_mape_x1 <- numeric(N)
FEsT_mape_x1 <- numeric(N)
FElT_mape_x1 <- numeric(N)
ME_mape_x1 <- numeric(N)
MEl_mape_x1 <- numeric(N)

REmse_x2 <- numeric(N)
FEmse_x2 <- numeric(N)
FE2mse_x2 <- numeric(N)
REl_mse_x2 <- numeric(N)
FEl_mse_x2 <- numeric(N)
FEt_mse_x2 <- numeric(N)
FEsT_mse_x2 <- numeric(N)
FElT_mse_x2 <- numeric(N)
ME_mse_x2 <- numeric(N)
MEl_mse_x2 <- numeric(N)

REmae_x2 <- numeric(N)
FEmae_x2 <- numeric(N)
FE2mae_x2 <- numeric(N)
REl_mae_x2 <- numeric(N)
FEl_mae_x2 <- numeric(N)
FEt_mae_x2 <- numeric(N)
FEsT_mae_x2 <- numeric(N)
FElT_mae_x2 <- numeric(N)
ME_mae_x2 <- numeric(N)
MEl_mae_x2 <- numeric(N)

REmpe_x2 <- numeric(N)
FEmpe_x2 <- numeric(N)
FE2mpe_x2 <- numeric(N)
REl_mpe_x2 <- numeric(N)
FEl_mpe_x2 <- numeric(N)
FEt_mpe_x2 <- numeric(N)
FEsT_mpe_x2 <- numeric(N)
FElT_mpe_x2 <- numeric(N)
ME_mpe_x2 <- numeric(N)
MEl_mpe_x2 <- numeric(N)

REmape_x2 <- numeric(N)
FEmape_x2 <- numeric(N)
FE2mape_x2 <- numeric(N)
REl_mape_x2 <- numeric(N)
FEl_mape_x2 <- numeric(N)
FEt_mape_x2 <- numeric(N)
FEsT_mape_x2 <- numeric(N)
FElT_mape_x2 <- numeric(N)
ME_mape_x2 <- numeric(N)
MEl_mape_x2 <- numeric(N)

REmse_x3 <- numeric(N)
FEmse_x3 <- numeric(N)
FE2mse_x3 <- numeric(N)
REl_mse_x3 <- numeric(N)
FEl_mse_x3 <- numeric(N)
FEt_mse_x3 <- numeric(N)
FEsT_mse_x3 <- numeric(N)
FElT_mse_x3 <- numeric(N)
ME_mse_x3 <- numeric(N)
MEl_mse_x3 <- numeric(N)

REmae_x3 <- numeric(N)
FEmae_x3 <- numeric(N)
FE2mae_x3 <- numeric(N)
REl_mae_x3 <- numeric(N)
FEl_mae_x3 <- numeric(N)
FEt_mae_x3 <- numeric(N)
FEsT_mae_x3 <- numeric(N)
FElT_mae_x3 <- numeric(N)
ME_mae_x3 <- numeric(N)
MEl_mae_x3 <- numeric(N)

REmpe_x3 <- numeric(N)
FEmpe_x3 <- numeric(N)
FE2mpe_x3 <- numeric(N)
REl_mpe_x3 <- numeric(N)
FEl_mpe_x3 <- numeric(N)
FEt_mpe_x3 <- numeric(N)
FEsT_mpe_x3 <- numeric(N)
FElT_mpe_x3 <- numeric(N)
ME_mpe_x3 <- numeric(N)
MEl_mpe_x3 <- numeric(N)

REmape_x3 <- numeric(N)
FEmape_x3 <- numeric(N)
FE2mape_x3 <- numeric(N)
REl_mape_x3 <- numeric(N)
FEl_mape_x3 <- numeric(N)
FEt_mape_x3 <- numeric(N)
FEsT_mape_x3 <- numeric(N)
FElT_mape_x3 <- numeric(N)
ME_mape_x3 <- numeric(N)
MEl_mape_x3 <- numeric(N)

FEsT_mse_trend <- numeric(N)
FElT_mse_trend <- numeric(N)
FEsT_mae_trend <- numeric(N)
FElT_mae_trend <- numeric(N)
FEsT_mpe_trend <- numeric(N)
FElT_mpe_trend <- numeric(N)
FEsT_mape_trend <- numeric(N)
FElT_mape_trend <- numeric(N)

CIlo_REb0 <- numeric(N)
CIlo_FEb0 <- numeric(N)
CIlo_FE2b0 <- numeric(N)
CIlo_REl_b0 <- numeric(N)
CIlo_FEl_b0 <- numeric(N)
CIlo_FEt_b0 <- numeric(N)
CIlo_FEsT_b0 <- numeric(N)
CIlo_FElT_b0 <- numeric(N)
CIlo_MEb0 <- numeric(N)
CIlo_MEl_b0 <- numeric(N)
CIhi_REb0 <- numeric(N)
CIhi_FEb0 <- numeric(N)
CIhi_FE2b0 <- numeric(N)
CIhi_REl_b0 <- numeric(N)
CIhi_FEl_b0 <- numeric(N)
CIhi_FEt_b0 <- numeric(N)
CIhi_FEsT_b0 <- numeric(N)
CIhi_FElT_b0 <- numeric(N)
CIhi_MEb0 <- numeric(N)
CIhi_MEl_b0 <- numeric(N)

CIlo_REb1 <- numeric(N)
CIlo_FEb1 <- numeric(N)
CIlo_FE2b1 <- numeric(N)
CIlo_REl_b1 <- numeric(N)
CIlo_FEl_b1 <- numeric(N)
CIlo_FEt_b1 <- numeric(N)
CIlo_FEsT_b1 <- numeric(N)
CIlo_FElT_b1 <- numeric(N)
CIlo_MEb1 <- numeric(N)
CIlo_MEl_b1 <- numeric(N)
CIhi_REb1 <- numeric(N)
CIhi_FEb1 <- numeric(N)
CIhi_FE2b1 <- numeric(N)
CIhi_REl_b1 <- numeric(N)
CIhi_FEl_b1 <- numeric(N)
CIhi_FEt_b1 <- numeric(N)
CIhi_FEsT_b1 <- numeric(N)
CIhi_FElT_b1 <- numeric(N)
CIhi_MEb1 <- numeric(N)
CIhi_MEl_b1 <- numeric(N)

CIlo_REb2 <- numeric(N)
CIlo_FEb2 <- numeric(N)
CIlo_FE2b2 <- numeric(N)
CIlo_REl_b2 <- numeric(N)
CIlo_FEl_b2 <- numeric(N)
CIlo_FEt_b2 <- numeric(N)
CIlo_FEsT_b2 <- numeric(N)
CIlo_FElT_b2 <- numeric(N)
CIlo_MEb2 <- numeric(N)
CIlo_MEl_b2 <- numeric(N)
CIhi_REb2 <- numeric(N)
CIhi_FEb2 <- numeric(N)
CIhi_FE2b2 <- numeric(N)
CIhi_REl_b2 <- numeric(N)
CIhi_FEl_b2 <- numeric(N)
CIhi_FEt_b2 <- numeric(N)
CIhi_FEsT_b2 <- numeric(N)
CIhi_FElT_b2 <- numeric(N)
CIhi_MEb2 <- numeric(N)
CIhi_MEl_b2 <- numeric(N)

CIlo_REb3 <- numeric(N)
CIlo_FEb3 <- numeric(N)
CIlo_FE2b3 <- numeric(N)
CIlo_REl_b3 <- numeric(N)
CIlo_FEl_b3 <- numeric(N)
CIlo_FEt_b3 <- numeric(N)
CIlo_FEsT_b3 <- numeric(N)
CIlo_FElT_b3 <- numeric(N)
CIlo_MEb3 <- numeric(N)
CIlo_MEl_b3 <- numeric(N)
CIhi_REb3 <- numeric(N)
CIhi_FEb3 <- numeric(N)
CIhi_FE2b3 <- numeric(N)
CIhi_REl_b3 <- numeric(N)
CIhi_FEl_b3 <- numeric(N)
CIhi_FEt_b3 <- numeric(N)
CIhi_FEsT_b3 <- numeric(N)
CIhi_FElT_b3 <- numeric(N)
CIhi_MEb3 <- numeric(N)
CIhi_MEl_b3 <- numeric(N)

CIlo_FEsT_bTrend <- numeric(N)
CIhi_FEsT_bTrend <- numeric(N)
CIlo_FElT_bTrend <- numeric(N)
CIhi_FElT_bTrend <- numeric(N)

start<-Sys.time()
for (k in 1:N) {
  # To show which iteration you're on
  cat(k,"/",N,"\n")
  
  #################
  ### SETUP     ###
  #################

  ###### COUNTRY 1 ######
  
  ### 2020 ###
  mu_y <- spreadMuY
  mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # The mean for y changes, x mean remains
  
  z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA) # z acts as a temporary placeholder
  y <- z_data[, 1]
  x <- z_data[, 2:4]

  X <- data.frame(constant = 1, x = x, year = 2020, country = 1, paper = 1, trend = -1)
    
  ### 2021-2024 ###
    
    # 2020 is not in the loop since it's the starting point
  
  for (i in 1:(numYears-1)) {
    mu_y <- mu_y+(timeHet*i) # The mean by _ per year
    mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x mean remains
    
    z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
    y2 <- z_data[, 1]
    x <- z_data[, 2:4]

    X2 <- data.frame(constant = 1, x = x, year = i + 2020, country = 1, 
                     paper = i + 1, trend = -1+i/2)
    
    X <- bind_rows(X, X2)
    y <- c(y, y2)
  }

  
    
  ###### COUNTRY 2 - 5 (or 9) ######
    # Country 1 is not in the loop since it's the starting point
  
  ### 2020-2024 ###
  paper_iterator <- numYears+1 # since Country 1 includes 1 paper per year
  countryIndex=2
  
  if(numCountries==5 & spreadMuY==-2){
    for (i in 2:numCountries) {
      for (j in 0:(numYears-1)) {
        mu_y <- i-3+(timeHet*j) # The mean for y increases by 1 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x means remain
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
    
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)
            
        X <- bind_rows(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
    }
  }
  else if(numCountries==5 & spreadMuY==-10){
    for (i in seq(from=-5, to=10, by=5)) { # (-10), -5, ..., 10
      for (j in 0:(numYears-1)) {
        mu_y <- i+(timeHet*j) # The mean for y increases by 5 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
        
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)
        
        X <- rbind(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
      countryIndex=countryIndex+1
    }
  }
  else if(numCountries==9 & spreadMuY==-2){
    for (i in seq(from=-1.5, to=2, by=0.5)) { # (-2), -1.5, -1, ..., 2
      for (j in 0:(numYears-1)) {
        mu_y <- i+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
        
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)

        X <- rbind(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
      countryIndex=countryIndex+1
    }
  }
  else if(numCountries==9 & spreadMuY==-10){
    for (i in seq(from=-7.5, to=10, by=2.5)) { # (-10), -7.5, ..., 10
      for (j in 0:(numYears-1)) {
        mu_y <- i+(timeHet*j) # The mean for y increases by 2.5 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
        
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)
        
        X <- rbind(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
      countryIndex=countryIndex+1
    }
  }
  else if(numCountries==15 & spreadMuY==-2){
    temp=c(-1.66,-1.33,-1,-0.66,-0.33,-0.15,0,0.15,0.33,0.66,1,1.33,1.66,2) # for calculating spreadMuY
    for (i in 1:(length(temp))){
      for (j in 0:(numYears-1)) {
        mu_y <- temp[i]+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
        
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)
        
        X <- rbind(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
      countryIndex=countryIndex+1
    }
  }
  else if(numCountries==15 & spreadMuY==-10){
    temp=c(-8.5,-7,-5.5,-4,-2.5,-1,0,1,2.5,4,5.5,7,8.5,10) # for calculating spreadMuY
    for (i in 1:(length(temp))){
      for (j in 0:(numYears-1)) {
        mu_y <- temp[i]+(timeHet*j) # The mean for y increases by 0.5 per country, also increases by _ per year
        mu_0 <- c(mu_y, mu_x1, mu_x2, mu_x3) # x remains mean 0
        
        z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
        x <- z_data[, 2:4]
        y2 <- z_data[, 1]
        
        X2 <- data.frame(constant = 1, x = x, year = j + 2020, country = i,
                         paper = paper_iterator, trend = -1+j/2)
      
        X <- rbind(X, X2)
        y <- c(y, y2)
        paper_iterator <- paper_iterator + 1
      }
      countryIndex=countryIndex+1
    }
  }
  
  
  
  ##################
  ### RUN MODELS ###
  ##################
  
  # To ensure NAs don't keep the regressions from running
  y <- na.omit(y)
  X <- na.omit(X)
  
  ### Random Effects (Study Level) ###
  if(model1==TRUE){
    RE <- plm(y ~ x.1 + x.2 + x.3, data=X, index=c("paper"), model="random")  #random model
    #summary(RE)
    
    REr2[k] <- summary(RE)$r.squared[1]
    REa_r2[k] <- summary(RE)$r.squared[2]
    REmse[k] <- mean(RE$residuals^2)                # MSE (df adjusted)
    REmae[k] <- mean(abs(RE$residuals))             # MAE (df adjusted)
    REmpe[k] <- mean(RE$residuals/RE$model$y)       # MPE (df adjusted)
    REmape[k] <- abs(REmpe[k])
    
    REb0[k] <- coef(RE)[1]
    REb1[k] <- coef(RE)[2]
    REb2[k] <- coef(RE)[3]
    REb3[k] <- coef(RE)[4]
    
    REmse_x1[k] <- mean((REb1[k]-b1True)^2)         # MSE (df adjusted)
    REmae_x1[k] <- mean(abs(REb1[k]-b1True))        # MAE (df adjusted)
    REmpe_x1[k] <- mean((b1True-REb1[k])/b1True)   # MPE (df adjusted)
    REmape_x1[k] <- abs(REmpe_x1[k])
    REmse_x2[k] <- mean((REb2[k]-b2True)^2)         # MSE (df adjusted)
    REmae_x2[k] <- mean(abs(REb2[k]-b2True))        # MAE (df adjusted)
    REmpe_x2[k] <- mean((b2True-REb2[k])/b2True)   # MPE (df adjusted)
    REmape_x2[k] <- abs(REmpe_x2[k])
    REmse_x3[k] <- mean((REb3[k]-b3True)^3)         # MSE (df adjusted)
    REmae_x3[k] <- mean(abs(REb3[k]-b3True))        # MAE (df adjusted)
    REmpe_x3[k] <- mean((b3True-REb3[k])/b3True)   # MPE (df adjusted)
    REmape_x3[k] <- abs(REmpe_x3[k])
    
    REsig2[k] <- sigma(RE)^2
    
    REb0SE[k] <- sqrt(diag(vcov(RE)))[1]
    REb1SE[k] <- sqrt(diag(vcov(RE)))[2]
    REb2SE[k] <- sqrt(diag(vcov(RE)))[3]
    REb3SE[k] <- sqrt(diag(vcov(RE)))[4]
    
    CIlo_REb0[k] <- REb0[k] - (1.96*REb0SE[k])
    CIhi_REb0[k] <- REb0[k] + (1.96*REb0SE[k])
    CIlo_REb1[k] <- REb1[k] - (1.96*REb1SE[k])
    CIhi_REb1[k] <- REb1[k] + (1.96*REb1SE[k])
    CIlo_REb2[k] <- REb2[k] - (1.96*REb2SE[k])
    CIhi_REb2[k] <- REb2[k] + (1.96*REb2SE[k])
    CIlo_REb3[k] <- REb3[k] - (1.96*REb3SE[k])
    CIhi_REb3[k] <- REb3[k] + (1.96*REb3SE[k])
    
    p<-length(RE$coefficients)
    obs<-length(RE$residuals)
    REaic[k]<-2*p+obs*(log(2*pi)+sigma(RE)^2)
    
    varRE <- RE$ercomp$sigma2[2] # (between group variation)
    varFE <- RE$ercomp$sigma2[1] # (within group variation)
    RE_i2[k] <- varRE / (varFE+varRE) # Var RE / Var TOTAL
  }
  
  ### Study-level Fixed Effects ###
  if(model3==TRUE){
    FE <- lm(y ~ x.1 + x.2 + x.3 + as.factor(paper), data=X)  #fixed model
    FEresult <- summary(FE)
  
    FEr2[k] <- summary(FE)$r.squared
    FEa_r2[k] <- summary(FE)$adj.r.squared
    FEmse[k] <- mean(FE$residuals^2)                # MSE (df adjusted)
    FEmae[k] <- mean(abs(FE$residuals))             # MAE (df adjusted)
    FEmpe[k] <- mean(FE$residuals/FE$model$y)       # MPE (df adjusted)
    FEmape[k] <- abs(FEmpe[k])
    
    FEb0[k] <- coef(FE)[1]
    FEb1[k] <- coef(FE)[2]
    FEb2[k] <- coef(FE)[3]
    FEb3[k] <- coef(FE)[4]
    
    FEmse_x1[k] <- mean((FEb1[k]-b1True)^2)         # MSE (df adjusted)
    FEmae_x1[k] <- mean(abs(FEb1[k]-b1True))        # MAE (df adjusted)
    FEmpe_x1[k] <- mean((b1True-FEb1[k])/b1True)   # MPE (df adjusted)
    FEmape_x1[k] <- abs(FEmpe_x1[k])
    FEmse_x2[k] <- mean((FEb2[k]-b2True)^2)         # MSE (df adjusted)
    FEmae_x2[k] <- mean(abs(FEb2[k]-b2True))        # MAE (df adjusted)
    FEmpe_x2[k] <- mean((b2True-FEb2[k])/b2True)   # MPE (df adjusted)
    FEmape_x2[k] <- abs(FEmpe_x2[k])
    FEmse_x3[k] <- mean((FEb3[k]-b3True)^3)         # MSE (df adjusted)
    FEmae_x3[k] <- mean(abs(FEb3[k]-b3True))        # MAE (df adjusted)
    FEmpe_x3[k] <- mean((b3True-FEb3[k])/b3True)   # MPE (df adjusted)
    FEmape_x3[k] <- abs(FEmpe_x3[k])
    
    FEsig2[k] <- sigma(FE)^2
    
    FEb0SE[k] <- sqrt(diag(vcov(FE)))[1]
    FEb1SE[k] <- sqrt(diag(vcov(FE)))[2]
    FEb2SE[k] <- sqrt(diag(vcov(FE)))[3]
    FEb3SE[k] <- sqrt(diag(vcov(FE)))[4]
    
    CIlo_FEb0[k] <- FEb0[k] - (1.96*FEb0SE[k])
    CIhi_FEb0[k] <- FEb0[k] + (1.96*FEb0SE[k])
    CIlo_FEb1[k] <- FEb1[k] - (1.96*FEb1SE[k])
    CIhi_FEb1[k] <- FEb1[k] + (1.96*FEb1SE[k])
    CIlo_FEb2[k] <- FEb2[k] - (1.96*FEb2SE[k])
    CIhi_FEb2[k] <- FEb2[k] + (1.96*FEb2SE[k])
    CIlo_FEb3[k] <- FEb3[k] - (1.96*FEb3SE[k])
    CIhi_FEb3[k] <- FEb3[k] + (1.96*FEb3SE[k])
    
    p<-length(FE$coefficients)
    obs<-length(FE$residuals)
    FEaic[k]<-2*p+obs*(log(2*pi)+sigma(FE)^2)
    
    # I^2 (via FE)
    tempN <- length(FE$xlevels$`as.factor(paper)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FEresult$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FE$coefficients[numVars+i]
      var2[i] <- diag(vcov(FE))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FE_dfOfGroups <- tempN
    FE_Q[k] <- Q
    FE_i2[k] <- min(100,(Q-tempN/Q)*100)
    FE_h2[k] <- Q/tempN
  }
  
  ### Two Way Fixed Effects ###
  if(model6==TRUE){
    FE2 <- lm(y ~ x.1 + x.2 + x.3 + as.factor(country) + as.factor(year), data = X)
    FE2result <- summary(FE2)
    
    FE2r2[k] <- summary(FE2)$r.squared
    FE2a_r2[k] <- summary(FE2)$adj.r.squared
    FE2mse[k] <- mean(FE2$residuals^2)                # MSE (df adjusted)
    FE2mae[k] <- mean(abs(FE2$residuals))             # MAE (df adjusted)
    FE2mpe[k] <- mean(FE2$residuals/FE2$model$y)       # MPE (df adjusted)
    FE2mape[k] <- abs(FE2mpe[k])
    
    FE2b0[k] <- coef(FE2)[1]
    FE2b1[k] <- coef(FE2)[2]
    FE2b2[k] <- coef(FE2)[3]
    FE2b3[k] <- coef(FE2)[4]
    
    FE2mse_x1[k] <- mean((FE2b1[k]-b1True)^2)         # MSE (df adjusted)
    FE2mae_x1[k] <- mean(abs(FE2b1[k]-b1True))        # MAE (df adjusted)
    FE2mpe_x1[k] <- mean((b1True-FE2b1[k])/b1True)   # MPE (df adjusted)
    FE2mape_x1[k] <- abs(FE2mpe_x1[k])
    FE2mse_x2[k] <- mean((FE2b2[k]-b2True)^2)         # MSE (df adjusted)
    FE2mae_x2[k] <- mean(abs(FE2b2[k]-b2True))        # MAE (df adjusted)
    FE2mpe_x2[k] <- mean((b2True-FE2b2[k])/b2True)   # MPE (df adjusted)
    FE2mape_x2[k] <- abs(FE2mpe_x2[k])
    FE2mse_x3[k] <- mean((FE2b3[k]-b3True)^3)         # MSE (df adjusted)
    FE2mae_x3[k] <- mean(abs(FE2b3[k]-b3True))        # MAE (df adjusted)
    FE2mpe_x3[k] <- mean((b3True-FE2b3[k])/b3True)   # MPE (df adjusted)
    FE2mape_x3[k] <- abs(FE2mpe_x3[k])
    
    FE2sig2[k] <- sigma(FE2)^2
    
    FE2b0SE[k] <- sqrt(diag(vcov(FE2)))[1]
    FE2b1SE[k] <- sqrt(diag(vcov(FE2)))[2]
    FE2b2SE[k] <- sqrt(diag(vcov(FE2)))[3]
    FE2b3SE[k] <- sqrt(diag(vcov(FE2)))[4]
    
    CIlo_FE2b0[k] <- FE2b0[k] - (1.96*FE2b0SE[k])
    CIhi_FE2b0[k] <- FE2b0[k] + (1.96*FE2b0SE[k])
    CIlo_FE2b1[k] <- FE2b1[k] - (1.96*FE2b1SE[k])
    CIhi_FE2b1[k] <- FE2b1[k] + (1.96*FE2b1SE[k])
    CIlo_FE2b2[k] <- FE2b2[k] - (1.96*FE2b2SE[k])
    CIhi_FE2b2[k] <- FE2b2[k] + (1.96*FE2b2SE[k])
    CIlo_FE2b3[k] <- FE2b3[k] - (1.96*FE2b3SE[k])
    CIhi_FE2b3[k] <- FE2b3[k] + (1.96*FE2b3SE[k])
    
    p<-length(FE2$coefficients)
    obs<-length(FE2$residuals)
    FE2aic[k]<-2*p+obs*(log(2*pi)+sigma(FE2)^2)
    
    # I^2 (via FE2)
    tempN <- length(FE2$xlevels$`as.factor(country)`)-1+length(FE2$xlevels$`as.factor(year)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FE2result$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FE2$coefficients[numVars+i]
      var2[i] <- diag(vcov(FE2))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FE2_dfOfGroups <- tempN
    FE2_Q[k] <- Q
    FE2_i2[k] <- min(100,(Q-tempN/Q)*100)
    FE2_h2[k] <- Q/tempN
  }
  
  ### [LOCATION-level] Random Effects ###
  if(model2==TRUE){
    REl <- plm(y ~ x.1 + x.2 + x.3, data=X, index=c("country"), model="random")  #random model
    #summary(RE)
    
    REl_r2[k] <- summary(REl)$r.squared[1]
    REl_a_r2[k] <- summary(REl)$r.squared[2]
    REl_mse[k] <- mean(REl$residuals^2)                # MSE (df adjusted)
    REl_mae[k] <- mean(abs(REl$residuals))             # MAE (df adjusted)
    REl_mpe[k] <- mean(REl$residuals/REl$model$y)       # MPE (df adjusted)
    REl_mape[k] <- abs(REl_mpe[k])
    
    REl_b0[k] <- coef(REl)[1]
    REl_b1[k] <- coef(REl)[2]
    REl_b2[k] <- coef(REl)[3]
    REl_b3[k] <- coef(REl)[4]
    
    REl_mse_x1[k] <- mean((REl_b1[k]-b1True)^2)         # MSE (df adjusted)
    REl_mae_x1[k] <- mean(abs(REl_b1[k]-b1True))        # MAE (df adjusted)
    REl_mpe_x1[k] <- mean((b1True-REl_b1[k])/b1True)   # MPE (df adjusted)
    REl_mape_x1[k] <- abs(REl_mpe_x1[k])
    REl_mse_x2[k] <- mean((REl_b2[k]-b2True)^2)         # MSE (df adjusted)
    REl_mae_x2[k] <- mean(abs(REl_b2[k]-b2True))        # MAE (df adjusted)
    REl_mpe_x2[k] <- mean((b2True-REl_b2[k])/b2True)   # MPE (df adjusted)
    REl_mape_x2[k] <- abs(REl_mpe_x2[k])
    REl_mse_x3[k] <- mean((REl_b3[k]-b3True)^3)         # MSE (df adjusted)
    REl_mae_x3[k] <- mean(abs(REl_b3[k]-b3True))        # MAE (df adjusted)
    REl_mpe_x3[k] <- mean((b3True-REl_b3[k])/b3True)   # MPE (df adjusted)
    REl_mape_x3[k] <- abs(REl_mpe_x3[k])
    
    REl_sig2[k] <- sigma(REl)^2
    
    REl_b0SE[k] <- sqrt(diag(vcov(REl)))[1]
    REl_b1SE[k] <- sqrt(diag(vcov(REl)))[2]
    REl_b2SE[k] <- sqrt(diag(vcov(REl)))[3]
    REl_b3SE[k] <- sqrt(diag(vcov(REl)))[4]
    
    CIlo_REl_b0[k] <- REl_b0[k] - (1.96*REl_b0SE[k])
    CIhi_REl_b0[k] <- REl_b0[k] + (1.96*REl_b0SE[k])
    CIlo_REl_b1[k] <- REl_b1[k] - (1.96*REl_b1SE[k])
    CIhi_REl_b1[k] <- REl_b1[k] + (1.96*REl_b1SE[k])
    CIlo_REl_b2[k] <- REl_b2[k] - (1.96*REl_b2SE[k])
    CIhi_REl_b2[k] <- REl_b2[k] + (1.96*REl_b2SE[k])
    CIlo_REl_b3[k] <- REl_b3[k] - (1.96*REl_b3SE[k])
    CIhi_REl_b3[k] <- REl_b3[k] + (1.96*REl_b3SE[k])
    
    p<-length(REl$coefficients)
    obs<-length(REl$residuals)
    REl_aic[k]<-2*p+obs*(log(2*pi)+sigma(REl)^2)
    
    varRE <- REl$ercomp$sigma2[2] # (between group variation)
    varFE <- REl$ercomp$sigma2[1] # (within group variation)
    REl_i2[k] <- varRE / (varFE+varRE) # Var RE / Var TOTAL
  }
  
  ### [LOCATION-level] Fixed Effects ###
  if(model4==TRUE){
    FEl <- lm(y ~ x.1 + x.2 + x.3 + as.factor(country), data=X)  #fixed model
    FElresult <- summary(FEl)
    
    FEl_r2[k] <- summary(FEl)$r.squared
    FEl_a_r2[k] <- summary(FEl)$adj.r.squared
    FEl_mse[k] <- mean(FEl$residuals^2)                # MSE (df adjusted)
    FEl_mae[k] <- mean(abs(FEl$residuals))             # MAE (df adjusted)
    FEl_mpe[k] <- mean(FEl$residuals/FEl$model$y)       # MPE (df adjusted)
    FEl_mape[k] <- abs(FEl_mpe[k])
    
    FEl_b0[k] <- coef(FEl)[1]
    FEl_b1[k] <- coef(FEl)[2]
    FEl_b2[k] <- coef(FEl)[3]
    FEl_b3[k] <- coef(FEl)[4]
    
    FEl_mse_x1[k] <- mean((FEl_b1[k]-b1True)^2)         # MSE (df adjusted)
    FEl_mae_x1[k] <- mean(abs(FEl_b1[k]-b1True))        # MAE (df adjusted)
    FEl_mpe_x1[k] <- mean((b1True-FEl_b1[k])/b1True)   # MPE (df adjusted)
    FEl_mape_x1[k] <- abs(FEl_mpe_x1[k])
    FEl_mse_x2[k] <- mean((FEl_b2[k]-b2True)^2)         # MSE (df adjusted)
    FEl_mae_x2[k] <- mean(abs(FEl_b2[k]-b2True))        # MAE (df adjusted)
    FEl_mpe_x2[k] <- mean((b2True-FEl_b2[k])/b2True)   # MPE (df adjusted)
    FEl_mape_x2[k] <- abs(FEl_mpe_x2[k])
    FEl_mse_x3[k] <- mean((FEl_b3[k]-b3True)^3)         # MSE (df adjusted)
    FEl_mae_x3[k] <- mean(abs(FEl_b3[k]-b3True))        # MAE (df adjusted)
    FEl_mpe_x3[k] <- mean((b3True-FEl_b3[k])/b3True)   # MPE (df adjusted)
    FEl_mape_x3[k] <- abs(FEl_mpe_x3[k])
    
    FEl_sig2[k] <- sigma(FEl)^2
    
    FEl_b0SE[k] <- sqrt(diag(vcov(FEl)))[1]
    FEl_b1SE[k] <- sqrt(diag(vcov(FEl)))[2]
    FEl_b2SE[k] <- sqrt(diag(vcov(FEl)))[3]
    FEl_b3SE[k] <- sqrt(diag(vcov(FEl)))[4]
    
    CIlo_FEl_b0[k] <- FEl_b0[k] - (1.96*FEl_b0SE[k])
    CIhi_FEl_b0[k] <- FEl_b0[k] + (1.96*FEl_b0SE[k])
    CIlo_FEl_b1[k] <- FEl_b1[k] - (1.96*FEl_b1SE[k])
    CIhi_FEl_b1[k] <- FEl_b1[k] + (1.96*FEl_b1SE[k])
    CIlo_FEl_b2[k] <- FEl_b2[k] - (1.96*FEl_b2SE[k])
    CIhi_FEl_b2[k] <- FEl_b2[k] + (1.96*FEl_b2SE[k])
    CIlo_FEl_b3[k] <- FEl_b3[k] - (1.96*FEl_b3SE[k])
    CIhi_FEl_b3[k] <- FEl_b3[k] + (1.96*FEl_b3SE[k])
    
    p<-length(FEl$coefficients)
    obs<-length(FEl$residuals)
    FEl_aic[k]<-2*p+obs*(log(2*pi)+sigma(FEl)^2)
    
    # I^2 (via FEl)
    tempN <- length(FEl$xlevels$`as.factor(country)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FElresult$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FEl$coefficients[numVars+i]
      var2[i] <- diag(vcov(FEl))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FEl_dfOfGroups <- tempN
    FEl_Q[k] <- Q
    FEl_i2[k] <- min(100,(Q-tempN/Q)*100)
    FEl_h2[k] <- Q/tempN
  }
  
  ### [TIME-level] Fixed Effects ###
  if(model5==TRUE){
    FEt <- lm(y ~ x.1 + x.2 + x.3 + as.factor(year), data=X)  #fixed model
    FEtresult <- summary(FEt)
    
    FEt_r2[k] <- summary(FEt)$r.squared
    FEt_a_r2[k] <- summary(FEt)$adj.r.squared
    FEt_mse[k] <- mean(FEt$residuals^2)                # MSE (df adjusted)
    FEt_mae[k] <- mean(abs(FEt$residuals))             # MAE (df adjusted)
    FEt_mpe[k] <- mean(FEt$residuals/FEt$model$y)       # MPE (df adjusted)
    FEt_mape[k] <- abs(FEt_mpe[k])
    
    FEt_b0[k] <- coef(FEt)[1]
    FEt_b1[k] <- coef(FEt)[2]
    FEt_b2[k] <- coef(FEt)[3]
    FEt_b3[k] <- coef(FEt)[4]
    
    FEt_mse_x1[k] <- mean((FEt_b1[k]-b1True)^2)         # MSE (df adjusted)
    FEt_mae_x1[k] <- mean(abs(FEt_b1[k]-b1True))        # MAE (df adjusted)
    FEt_mpe_x1[k] <- mean((b1True-FEt_b1[k])/b1True)   # MPE (df adjusted)
    FEt_mape_x1[k] <- abs(FEt_mpe_x1[k])
    FEt_mse_x2[k] <- mean((FEt_b2[k]-b2True)^2)         # MSE (df adjusted)
    FEt_mae_x2[k] <- mean(abs(FEt_b2[k]-b2True))        # MAE (df adjusted)
    FEt_mpe_x2[k] <- mean((b2True-FEt_b2[k])/b2True)   # MPE (df adjusted)
    FEt_mape_x2[k] <- abs(FEt_mpe_x2[k])
    FEt_mse_x3[k] <- mean((FEt_b3[k]-b3True)^3)         # MSE (df adjusted)
    FEt_mae_x3[k] <- mean(abs(FEt_b3[k]-b3True))        # MAE (df adjusted)
    FEt_mpe_x3[k] <- mean((b3True-FEt_b3[k])/b3True)   # MPE (df adjusted)
    FEt_mape_x3[k] <- abs(FEt_mpe_x3[k])
    
    FEt_sig2[k] <- sigma(FEt)^2
    
    FEt_b0SE[k] <- sqrt(diag(vcov(FEt)))[1]
    FEt_b1SE[k] <- sqrt(diag(vcov(FEt)))[2]
    FEt_b2SE[k] <- sqrt(diag(vcov(FEt)))[3]
    FEt_b3SE[k] <- sqrt(diag(vcov(FEt)))[4]
    
    CIlo_FEt_b0[k] <- FEt_b0[k] - (1.96*FEt_b0SE[k])
    CIhi_FEt_b0[k] <- FEt_b0[k] + (1.96*FEt_b0SE[k])
    CIlo_FEt_b1[k] <- FEt_b1[k] - (1.96*FEt_b1SE[k])
    CIhi_FEt_b1[k] <- FEt_b1[k] + (1.96*FEt_b1SE[k])
    CIlo_FEt_b2[k] <- FEt_b2[k] - (1.96*FEt_b2SE[k])
    CIhi_FEt_b2[k] <- FEt_b2[k] + (1.96*FEt_b2SE[k])
    CIlo_FEt_b3[k] <- FEt_b3[k] - (1.96*FEt_b3SE[k])
    CIhi_FEt_b3[k] <- FEt_b3[k] + (1.96*FEt_b3SE[k])
    
    p<-length(FEt$coefficients)
    obs<-length(FEt$residuals)
    FEt_aic[k]<-2*p+obs*(log(2*pi)+sigma(FEt)^2)
    
    # I^2 (via FEt)
    tempN <- length(FEt$xlevels$`as.factor(year)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FEtresult$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FEt$coefficients[numVars+i]
      var2[i] <- diag(vcov(FEt))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FEt_dfOfGroups <- tempN
    FEt_Q[k] <- Q
    FEt_i2[k] <- min(100,(Q-tempN/Q)*100)
    FEt_h2[k] <- Q/tempN
  }
  ### [Study-level] Fixed Effects (with Trend) ###
  if(model7==TRUE){
    FEsT <- lm(y ~ x.1 + x.2 + x.3 + trend + as.factor(paper), data=X)  #fixed model
    FEsTresult <- summary(FEsT)
    
    FEsT_r2[k] <- summary(FEsT)$r.squared
    FEsT_a_r2[k] <- summary(FEsT)$adj.r.squared
    FEsT_mse[k] <- mean(FEsT$residuals^2)                # MSE (df adjusted)
    FEsT_mae[k] <- mean(abs(FEsT$residuals))             # MAE (df adjusted)
    FEsT_mpe[k] <- mean(FEsT$residuals/FEsT$model$y)       # MPE (df adjusted)
    FEsT_mape[k] <- abs(FEsT_mpe[k])
    
    FEsT_b0[k] <- coef(FEsT)[1]
    FEsT_b1[k] <- coef(FEsT)[2]
    FEsT_b2[k] <- coef(FEsT)[3]
    FEsT_b3[k] <- coef(FEsT)[4]
    FEsT_bTrend[k] <- coef(FEsT)[5]
    
    FEsT_mse_x1[k] <- mean((FEsT_b1[k]-b1True)^2)         # MSE (df adjusted)
    FEsT_mae_x1[k] <- mean(abs(FEsT_b1[k]-b1True))        # MAE (df adjusted)
    FEsT_mpe_x1[k] <- mean((b1True-FEsT_b1[k])/b1True)   # MPE (df adjusted)
    FEsT_mape_x1[k] <- abs(FEsT_mpe_x1[k])
    FEsT_mse_trend[k] <- mean((FEsT_bTrend[k]-timeHet)^2)         # MSE (df adjusted)
    FEsT_mae_trend[k] <- mean(abs(FEsT_bTrend[k]-timeHet))        # MAE (df adjusted)
    FEsT_mpe_trend[k] <- mean((timeHet-FEsT_bTrend[k])/timeHet)   # MPE (df adjusted)
    FEsT_mape_trend[k] <- abs(FEsT_mpe_trend[k])
    FEsT_mse_x2[k] <- mean((FEsT_b2[k]-b2True)^2)         # MSE (df adjusted)
    FEsT_mae_x2[k] <- mean(abs(FEsT_b2[k]-b2True))        # MAE (df adjusted)
    FEsT_mpe_x2[k] <- mean((b2True-FEsT_b2[k])/b2True)   # MPE (df adjusted)
    FEsT_mape_x2[k] <- abs(FEsT_mpe_x2[k])
    FEsT_mse_x3[k] <- mean((FEsT_b3[k]-b3True)^3)         # MSE (df adjusted)
    FEsT_mae_x3[k] <- mean(abs(FEsT_b3[k]-b3True))        # MAE (df adjusted)
    FEsT_mpe_x3[k] <- mean((b3True-FEsT_b3[k])/b3True)   # MPE (df adjusted)
    FEsT_mape_x3[k] <- abs(FEsT_mpe_x3[k])
    
    FEsT_sig2[k] <- sigma(FEsT)^2
    
    FEsT_b0SE[k] <- sqrt(diag(vcov(FEsT)))[1]
    FEsT_b1SE[k] <- sqrt(diag(vcov(FEsT)))[2]
    FEsT_b2SE[k] <- sqrt(diag(vcov(FEsT)))[3]
    FEsT_b3SE[k] <- sqrt(diag(vcov(FEsT)))[4]
    FEsT_bTrendSE[k] <- sqrt(diag(vcov(FEsT)))[5]

    CIlo_FEsT_b0[k] <- FEsT_b0[k] - (1.96*FEsT_b0SE[k])
    CIhi_FEsT_b0[k] <- FEsT_b0[k] + (1.96*FEsT_b0SE[k])
    CIlo_FEsT_b1[k] <- FEsT_b1[k] - (1.96*FEsT_b1SE[k])
    CIhi_FEsT_b1[k] <- FEsT_b1[k] + (1.96*FEsT_b1SE[k])
    CIlo_FEsT_b2[k] <- FEsT_b2[k] - (1.96*FEsT_b2SE[k])
    CIhi_FEsT_b2[k] <- FEsT_b2[k] + (1.96*FEsT_b2SE[k])
    CIlo_FEsT_b3[k] <- FEsT_b3[k] - (1.96*FEsT_b3SE[k])
    CIhi_FEsT_b3[k] <- FEsT_b3[k] + (1.96*FEsT_b3SE[k])
    CIlo_FEsT_bTrend[k] <- FEsT_bTrend[k] - (1.96*FEsT_bTrendSE[k])
    CIhi_FEsT_bTrend[k] <- FEsT_bTrend[k] + (1.96*FEsT_bTrendSE[k])
        
    p<-length(FEsT$coefficients)
    obs<-length(FEsT$residuals)
    FEsT_aic[k]<-2*p+obs*(log(2*pi)+sigma(FEsT)^2)
    
    # I^2 (via FEsT)
    tempN <- length(FEsT$xlevels$`as.factor(paper)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FEsTresult$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FEsT$coefficients[numVars+i]
      var2[i] <- diag(vcov(FEsT))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FEsT_dfOfGroups <- tempN
    FEsT_Q[k] <- Q
    FEsT_i2[k] <- min(100,(Q-tempN/Q)*100)
    FEsT_h2[k] <- Q/tempN
  }
  
  ### [Location-level] Fixed Effects (with Trend) ###
  if(model8==TRUE){
    FElT <- lm(y ~ x.1 + x.2 + x.3 + trend + as.factor(country), data=X)  #fixed model
    FElTresult <- summary(FElT)
    
    FElT_r2[k] <- summary(FElT)$r.squared
    FElT_a_r2[k] <- summary(FElT)$adj.r.squared
    FElT_mse[k] <- mean(FElT$residuals^2)                # MSE (df adjusted)
    FElT_mae[k] <- mean(abs(FElT$residuals))             # MAE (df adjusted)
    FElT_mpe[k] <- mean(FElT$residuals/FElT$model$y)       # MPE (df adjusted)
    FElT_mape[k] <- abs(FElT_mpe[k])
    
    FElT_b0[k] <- coef(FElT)[1]
    FElT_b1[k] <- coef(FElT)[2]
    FElT_b2[k] <- coef(FElT)[3]
    FElT_b3[k] <- coef(FElT)[4]
    FElT_bTrend[k] <- coef(FElT)[5]

    FElT_mse_x1[k] <- mean((FElT_b1[k]-b1True)^2)         # MSE (df adjusted)
    FElT_mae_x1[k] <- mean(abs(FElT_b1[k]-b1True))        # MAE (df adjusted)
    FElT_mpe_x1[k] <- mean((b1True-FElT_b1[k])/b1True)   # MPE (df adjusted)
    FElT_mape_x1[k] <- abs(FElT_mpe_x1[k])
    FElT_mse_trend[k] <- mean((FElT_bTrend[k]-timeHet)^2)         # MSE (df adjusted)
    FElT_mae_trend[k] <- mean(abs(FElT_bTrend[k]-timeHet))        # MAE (df adjusted)
    FElT_mpe_trend[k] <- mean((timeHet-FElT_bTrend[k])/timeHet)   # MPE (df adjusted)
    FElT_mape_trend[k] <- abs(FElT_mpe_trend[k])
    FElT_mse_x2[k] <- mean((FElT_b2[k]-b2True)^2)         # MSE (df adjusted)
    FElT_mae_x2[k] <- mean(abs(FElT_b2[k]-b2True))        # MAE (df adjusted)
    FElT_mpe_x2[k] <- mean((b2True-FElT_b2[k])/b2True)   # MPE (df adjusted)
    FElT_mape_x2[k] <- abs(FElT_mpe_x2[k])
    FElT_mse_x3[k] <- mean((FElT_b3[k]-b3True)^3)         # MSE (df adjusted)
    FElT_mae_x3[k] <- mean(abs(FElT_b3[k]-b3True))        # MAE (df adjusted)
    FElT_mpe_x3[k] <- mean((b3True-FElT_b3[k])/b3True)   # MPE (df adjusted)
    FElT_mape_x3[k] <- abs(FElT_mpe_x3[k])
    
    FElT_sig2[k] <- sigma(FElT)^2
    
    FElT_b0SE[k] <- sqrt(diag(vcov(FElT)))[1]
    FElT_b1SE[k] <- sqrt(diag(vcov(FElT)))[2]
    FElT_b2SE[k] <- sqrt(diag(vcov(FElT)))[3]
    FElT_b3SE[k] <- sqrt(diag(vcov(FElT)))[4]
    FElT_bTrendSE[k] <- sqrt(diag(vcov(FElT)))[5]
    
    CIlo_FElT_b0[k] <- FElT_b0[k] - (1.96*FElT_b0SE[k])
    CIhi_FElT_b0[k] <- FElT_b0[k] + (1.96*FElT_b0SE[k])
    CIlo_FElT_b1[k] <- FElT_b1[k] - (1.96*FElT_b1SE[k])
    CIhi_FElT_b1[k] <- FElT_b1[k] + (1.96*FElT_b1SE[k])
    CIlo_FElT_b2[k] <- FElT_b2[k] - (1.96*FElT_b2SE[k])
    CIhi_FElT_b2[k] <- FElT_b2[k] + (1.96*FElT_b2SE[k])
    CIlo_FElT_b3[k] <- FElT_b3[k] - (1.96*FElT_b3SE[k])
    CIhi_FElT_b3[k] <- FElT_b3[k] + (1.96*FElT_b3SE[k])
    CIlo_FElT_bTrend[k] <- FElT_bTrend[k] - (1.96*FElT_bTrendSE[k])
    CIhi_FElT_bTrend[k] <- FElT_bTrend[k] + (1.96*FElT_bTrendSE[k])
    
    p<-length(FElT$coefficients)
    obs<-length(FElT$residuals)
    FElT_aic[k]<-2*p+obs*(log(2*pi)+sigma(FElT)^2)
    
    # I^2 (via FElT)
    tempN <- length(FElT$xlevels$`as.factor(country)`)-1 # !!! EDIT (if grouping != country)
    numVars <- nrow(FElTresult$coefficients)-tempN
    numSum <- 0
    denomSum <- 0
    yi <- numeric(tempN)
    var2 <- numeric(tempN)
    w <- numeric(tempN)
    
    for (i in 1:tempN){
      yi[i] <- FElT$coefficients[numVars+i]
      var2[i] <- diag(vcov(FElT))[numVars+i]
      w[i] <- 1/var2[i]
      
      numSum <- numSum + (w[i]*yi[i])
      denomSum <- denomSum + w[i]
    }
    poolFX <- numSum/denomSum
    #poolFX <- mean(y)
    
    Q <- 0
    for (i in 1:tempN){
      Q <- Q + (yi[i]-poolFX)^2 / var2[i]
    }
    # Values to Export
    FElT_dfOfGroups <- tempN
    FElT_Q[k] <- Q
    FElT_i2[k] <- min(100,(Q-tempN/Q)*100)
    FElT_h2[k] <- Q/tempN
  }
  
  ### Mixed Effects (RE @ Study-level) ###
  if(model9==TRUE){
    ME <- lmer(y ~ x.1 + x.2 + x.3 + (1 | paper), data=X)  #fixed model
    MEresult <- summary(ME)
    
    #ME_r2[k] <- summary(ME)$r.squared
    #ME_a_r2[k] <- summary(ME)$adj.r.squared
    MEmse[k] <- mean(residuals(ME)^2)                # MSE (df adjusted)
    MEmae[k] <- mean(abs(residuals(ME)))             # MAE (df adjusted)
    MEmpe[k] <- mean(residuals(ME)/y)            # MPE (df adjusted)
    MEmape[k] <- abs(MEmpe[k])
    
    MEb0[k] <- summary(ME)$coefficients[1]
    MEb1[k] <- summary(ME)$coefficients[2]
    MEb2[k] <- summary(ME)$coefficients[3]
    MEb3[k] <- summary(ME)$coefficients[4]
    
    ME_mse_x1[k] <- mean((MEb1[k]-b1True)^2)         # MSE (df adjusted)
    ME_mae_x1[k] <- mean(abs(MEb1[k]-b1True))        # MAE (df adjusted)
    ME_mpe_x1[k] <- mean((b1True-MEb1[k])/b1True)   # MPE (df adjusted)
    ME_mape_x1[k] <- abs(ME_mpe_x1[k])
    ME_mse_x2[k] <- mean((MEb2[k]-b2True)^2)         # MSE (df adjusted)
    ME_mae_x2[k] <- mean(abs(MEb2[k]-b2True))        # MAE (df adjusted)
    ME_mpe_x2[k] <- mean((b2True-MEb2[k])/b2True)   # MPE (df adjusted)
    ME_mape_x2[k] <- abs(ME_mpe_x2[k])
    ME_mse_x3[k] <- mean((MEb3[k]-b3True)^3)         # MSE (df adjusted)
    ME_mae_x3[k] <- mean(abs(MEb3[k]-b3True))        # MAE (df adjusted)
    ME_mpe_x3[k] <- mean((b3True-MEb3[k])/b3True)   # MPE (df adjusted)
    ME_mape_x3[k] <- abs(ME_mpe_x3[k])
    
    MEsig2[k] <- sigma(ME)^2
    
    MEb0SE[k] <- summary(ME)$coefficients[5]
    MEb1SE[k] <- summary(ME)$coefficients[6]
    MEb2SE[k] <- summary(ME)$coefficients[7]
    MEb3SE[k] <- summary(ME)$coefficients[8]
    
    CIlo_MEb0[k] <- MEb0[k] - (1.96*MEb0SE[k])
    CIhi_MEb0[k] <- MEb0[k] + (1.96*MEb0SE[k])
    CIlo_MEb1[k] <- MEb1[k] - (1.96*MEb1SE[k])
    CIhi_MEb1[k] <- MEb1[k] + (1.96*MEb1SE[k])
    CIlo_MEb2[k] <- MEb2[k] - (1.96*MEb2SE[k])
    CIhi_MEb2[k] <- MEb2[k] + (1.96*MEb2SE[k])
    CIlo_MEb3[k] <- MEb3[k] - (1.96*MEb3SE[k])
    CIhi_MEb3[k] <- MEb3[k] + (1.96*MEb3SE[k])
    
    p<-1+1
    obs<-length(y)
    ME_aic[k]<-2*p+obs*(log(2*pi)+sigma(ME)^2)
    
    varRE <- VarCorr(ME)$paper[1] # !!! EDIT (if grouping != paper)
    varFE <- 0
    for (i in 1:(ncol(MEresult$coefficients)-1)){
      varFE <- varFE + diag(as.matrix(vcov(ME)))[i]
    }
    ME_i2[k] <- varFE / (varFE+varRE) # Var FE / Var TOTAL
  }
  
  ### Mixed Effects (RE @ Country-level) ###
  if(model10==TRUE){
    MEl <- lmer(y ~ x.1 + x.2 + x.3 + (1 | country), data=X)  #fixed model
    MElresult <- summary(MEl)
    
    #MEl_r2[k] <- summary(MEl)$r.squared
    #MEl_a_r2[k] <- summary(MEl)$adj.r.squared
    MEl_mse[k] <- mean(residuals(MEl)^2)                        # MSE (df adjusted)
    MEl_mae[k] <- mean(abs(residuals(MEl)))             # MAE (df adjusted)
    MEl_mpe[k] <- mean(residuals(MEl)/y)            # MPE (df adjusted)
    MEl_mape[k] <- abs(MEl_mpe[k])
    
    MEl_b0[k] <- summary(MEl)$coefficients[1]
    MEl_b1[k] <- summary(MEl)$coefficients[2]
    MEl_b2[k] <- summary(MEl)$coefficients[3]
    MEl_b3[k] <- summary(MEl)$coefficients[4]
    
    MEl_mse_x1[k] <- mean((MEl_b1[k]-b1True)^2)         # MSE (df adjusted)
    MEl_mae_x1[k] <- mean(abs(MEl_b1[k]-b1True))        # MAE (df adjusted)
    MEl_mpe_x1[k] <- mean((b1True-MEl_b1[k])/b1True)   # MPE (df adjusted)
    MEl_mape_x1[k] <- abs(MEl_mpe_x1[k])
    MEl_mse_x2[k] <- mean((MEl_b2[k]-b2True)^2)         # MSE (df adjusted)
    MEl_mae_x2[k] <- mean(abs(MEl_b2[k]-b2True))        # MAE (df adjusted)
    MEl_mpe_x2[k] <- mean((b2True-MEl_b2[k])/b2True)   # MPE (df adjusted)
    MEl_mape_x2[k] <- abs(MEl_mpe_x2[k])
    MEl_mse_x3[k] <- mean((MEl_b3[k]-b3True)^3)         # MSE (df adjusted)
    MEl_mae_x3[k] <- mean(abs(MEl_b3[k]-b3True))        # MAE (df adjusted)
    MEl_mpe_x3[k] <- mean((b3True-MEl_b3[k])/b3True)   # MPE (df adjusted)
    MEl_mape_x3[k] <- abs(MEl_mpe_x3[k])
    
    MEl_sig2[k] <- sigma(MEl)^2
    
    MEl_b0SE[k] <- summary(MEl)$coefficients[5]
    MEl_b1SE[k] <- summary(MEl)$coefficients[6]
    MEl_b2SE[k] <- summary(MEl)$coefficients[7]
    MEl_b3SE[k] <- summary(MEl)$coefficients[8]
    
    CIlo_MEl_b0[k] <- MEl_b0[k] - (1.96*MEl_b0SE[k])
    CIhi_MEl_b0[k] <- MEl_b0[k] + (1.96*MEl_b0SE[k])
    CIlo_MEl_b1[k] <- MEl_b1[k] - (1.96*MEl_b1SE[k])
    CIhi_MEl_b1[k] <- MEl_b1[k] + (1.96*MEl_b1SE[k])
    CIlo_MEl_b2[k] <- MEl_b2[k] - (1.96*MEl_b2SE[k])
    CIhi_MEl_b2[k] <- MEl_b2[k] + (1.96*MEl_b2SE[k])
    CIlo_MEl_b3[k] <- MEl_b3[k] - (1.96*MEl_b3SE[k])
    CIhi_MEl_b3[k] <- MEl_b3[k] + (1.96*MEl_b3SE[k])
    
    p<-1+1
    obs<-length(y)
    MEl_aic[k]<-2*p+obs*(log(2*pi)+sigma(MEl)^2)
    
    varRE <- VarCorr(MEl)$country[1] # !!! EDIT (if grouping != paper)
    varFE <- 0
    for (i in 1:(ncol(MElresult$coefficients)-1)){
      varFE <- varFE + diag(as.matrix(vcov(MEl)))[i]
    }
    MEl_i2[k] <- varFE / (varFE+varRE) # Var FE / Var TOTAL
  }
}
# Audio to signal Monte Carlo completed
beep(sound = 5)

end<-Sys.time()
print(end-start)

#%%###########################################################################
###                             RESULTS - MONTE CARLO                      ###
##############################################################################

###### BOX & WHISKER PLOTS ######

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Box Whisker - MSE.png", width = 1200, height = 1600)
# 2. Create the plot
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
boxplot(REmse, FEmse, FE2mse, REl_mse, FEl_mse, FEsT_mse, FElT_mse, MEmse, MEl_mse,
        main = paste("n =",n),
        names = c("REs","FEs","FElt","REl","FEl","FEsT","FElT","MEs","MEl"),
        horizontal = TRUE,
        ylab="Method Specifications",
        xlab="MSE Distributions",
        # limits use double min/max functions to set lower upper bounds on the x axis
        ylim=c(min(min(REmse),min(FEmse),min(FE2mse),min(FEl_mse),min(FEsT_mse),min(FElT_mse),min(MEmse),min(MEl_mse)),0.3+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse))),
        bty="l" #to remove top and right border
        
        ### Extra specifications if you wish
        #at = c(1,2,4,5),
        #las = 2,
        #col = c("orange","red"),
        #border = "brown",
        #notch = TRUE
)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 1, paste(round(mean(REmse),2)," [",round(min(REmse),2),",",round(max(REmse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 2, paste(round(mean(FEmse),2)," [",round(min(FEmse),2),",",round(max(FEmse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 3, paste(round(mean(FE2mse),2)," [",round(min(FE2mse),2),",",round(max(FE2mse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 4.2, paste(round(mean(REl_mse),2)," [",round(min(REl_mse),2),",",round(max(REl_mse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 5, paste(round(mean(FEl_mse),2)," [",round(min(FEl_mse),2),",",round(max(FEl_mse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 6, paste(round(mean(FEsT_mse),2)," [",round(min(FEsT_mse),2),",",round(max(FEsT_mse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 7, paste(round(mean(FElT_mse),2)," [",round(min(FElT_mse),2),",",round(max(FElT_mse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 8, paste(round(mean(MEmse),2)," [",round(min(MEmse),2),",",round(max(MEmse),2),"]",sep=""),cex=2)
text(x = 0.2+max(max(REmse),max(FEmse),max(FE2mse),max(FEl_mse),max(FEsT_mse),max(FElT_mse),max(MEmse),max(MEl_mse)), y = 9, paste(round(mean(MEl_mse),2)," [",round(min(MEl_mse),2),",",round(max(MEl_mse),2),"]",sep=""),cex=2)
#abline(v=0,lty=5)
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Box Whisker - b1.png", width = 1200, height = 1600)
# 2. Create the plot
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
boxplot(REb1, FEb1, FE2b1, REl_b1, FEl_b1, FEsT_b1, FElT_b1, MEb1, MEl_b1,
        main = paste("n =",n),
        names = c("REs","FEs","FElt","REl","FEl","FEsT","FElT","MEs","MEl"),
        horizontal = TRUE,
        ylab="Method Specifications",
        xlab="x1 Coefficient Distributions",
        # limits use double min/max functions to set lower upper bounds on the x axis
        ylim=c(b1True-0.1,0.1+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1))),
        bty="l" #to remove top and right border
        
        ### Extra specifications if you wish
        #at = c(1,2,4,5),
        #las = 2,
        #col = c("orange","red"),
        #border = "brown",
        #notch = TRUE
)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 1, paste(round(mean(REb1),2)," [",round(min(REb1),2),",",round(max(REb1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 2, paste(round(mean(FEb1),2)," [",round(min(FEb1),2),",",round(max(FEb1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 3, paste(round(mean(FE2b1),2)," [",round(min(FE2b1),2),",",round(max(FE2b1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 4.1, paste(round(mean(REl_b1),2)," [",round(min(REl_b1),2),",",round(max(REl_b1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 5, paste(round(mean(FEl_b1),2)," [",round(min(FEl_b1),2),",",round(max(FEl_b1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 6, paste(round(mean(FEsT_b1),2)," [",round(min(FEsT_b1),2),",",round(max(FEsT_b1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 7, paste(round(mean(FElT_b1),2)," [",round(min(FElT_b1),2),",",round(max(FElT_b1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 8, paste(round(mean(MEb1),2)," [",round(min(MEb1),2),",",round(max(MEb1),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb1),max(FEb1),max(FE2b1),max(FEl_b1),max(FEsT_b1),max(FElT_b1),max(MEb1),max(MEl_b1)), y = 9, paste(round(mean(MEl_b1),2)," [",round(min(MEl_b1),2),",",round(max(MEl_b1),2),"]",sep=""),cex=2)
abline(v=b1True,lty=5)
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Box Whisker - b2.png", width = 1200, height = 1600)
# 2. Create the plot
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
boxplot(REb2, FEb2, FE2b2, REl_b2, FEl_b2, FEsT_b2, FElT_b2, MEb2, MEl_b2,
        main = paste("n =",n),
        names = c("REs","FEs","FElt","REl","FEl","FEsT","FElT","MEs","MEl"),
        horizontal = TRUE,
        ylab="Method Specifications",
        xlab="x2 Coefficient Distributions",
        # limits use double min/max functions to set lower upper bounds on the x axis
        ylim=c(b2True-0.1,0.1+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2))),
        bty="l" #to remove top and right border
        
        ### Extra specifications if you wish
        #at = c(1,2,4,5),
        #las = 2,
        #col = c("orange","red"),
        #border = "brown",
        #notch = TRUE
)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 1, paste(round(mean(REb2),2)," [",round(min(REb2),2),",",round(max(REb2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 2, paste(round(mean(FEb2),2)," [",round(min(FEb2),2),",",round(max(FEb2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 3, paste(round(mean(FE2b2),2)," [",round(min(FE2b2),2),",",round(max(FE2b2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 4.1, paste(round(mean(REl_b2),2)," [",round(min(REl_b2),2),",",round(max(REl_b2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 5, paste(round(mean(FEl_b2),2)," [",round(min(FEl_b2),2),",",round(max(FEl_b2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 6, paste(round(mean(FEsT_b2),2)," [",round(min(FEsT_b2),2),",",round(max(FEsT_b2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 7, paste(round(mean(FElT_b2),2)," [",round(min(FElT_b2),2),",",round(max(FElT_b2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 8, paste(round(mean(MEb2),2)," [",round(min(MEb2),2),",",round(max(MEb2),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb2),max(FEb2),max(FE2b2),max(FEl_b2),max(FEsT_b2),max(FElT_b2),max(MEb2),max(MEl_b2)), y = 9, paste(round(mean(MEl_b2),2)," [",round(min(MEl_b2),2),",",round(max(MEl_b2),2),"]",sep=""),cex=2)
abline(v=b2True,lty=5)
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Box Whisker - b3.png", width = 1200, height = 1600)
# 2. Create the plot
par(cex.lab=2,cex.axis=2,cex.main=2,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
boxplot(REb3, FEb3, FE2b3, REl_b3, FEl_b3, FEsT_b3, FElT_b3, MEb3, MEl_b3,
        main = paste("n =",n),
        names = c("REs","FEs","FElt","REl","FEl","FEsT","FElT","MEs","MEl"),
        horizontal = TRUE,
        ylab="Method Specifications",
        xlab="x3 Coefficient Distributions",
        # limits use double min/max functions to set lower upper bounds on the x axis
        ylim=c(b3True-0.1,0.1+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3))),
        bty="l" #to remove top and right border
        
        ### Extra specifications if you wish
        #at = c(1,2,4,5),
        #las = 2,
        #col = c("orange","red"),
        #border = "brown",
        #notch = TRUE
)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 1, paste(round(mean(REb3),2)," [",round(min(REb3),2),",",round(max(REb3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 2, paste(round(mean(FEb3),2)," [",round(min(FEb3),2),",",round(max(FEb3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 3, paste(round(mean(FE2b3),2)," [",round(min(FE2b3),2),",",round(max(FE2b3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 4.1, paste(round(mean(REl_b3),2)," [",round(min(REl_b3),2),",",round(max(REl_b3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 5, paste(round(mean(FEl_b3),2)," [",round(min(FEl_b3),2),",",round(max(FEl_b3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 6, paste(round(mean(FEsT_b3),2)," [",round(min(FEsT_b3),2),",",round(max(FEsT_b3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 7, paste(round(mean(FElT_b3),2)," [",round(min(FElT_b3),2),",",round(max(FElT_b3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 8, paste(round(mean(MEb3),2)," [",round(min(MEb3),2),",",round(max(MEb3),2),"]",sep=""),cex=2)
text(x = 0.05+max(max(REb3),max(FEb3),max(FE2b3),max(FEl_b3),max(FEsT_b3),max(FElT_b3),max(MEb3),max(MEl_b3)), y = 9, paste(round(mean(MEl_b3),2)," [",round(min(MEl_b3),2),",",round(max(MEl_b3),2),"]",sep=""),cex=2)
abline(v=b3True,lty=5)
dev.off()



###### HISTOGRAMS ######

# TRUE b1 specified in vertical lines
# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_x1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb1, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEb1, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEb1, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FE2b1, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(REl_b1, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEl_b1, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEl_b1, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEt_b1, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(FEsT_b1, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FElT_b1, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b1True, col = 'black', lwd = 3) #lty = 'dashed' 
mtext("x1 Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# TRUE b2 specified in vertical lines
# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEb2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEb2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FE2b2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(REl_b2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEl_b2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEl_b2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEt_b2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(FEsT_b2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FElT_b2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b2True, col = 'black', lwd = 3) #lty = 'dashed' 
mtext("x2 Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# TRUE b3 specified in vertical lines
# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb3, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEb3, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEb3, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FE2b3, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(REl_b3, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEl_b3, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEl_b3, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEt_b3, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(FEsT_b3, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FElT_b3, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = b3True, col = 'black', lwd = 3) #lty = 'dashed' 
mtext("x3 Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MSE.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmse, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmse, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEmse, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mse, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mse, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mse, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mse, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mse, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mse, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mse, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MSE Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAE.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmae, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmae, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEmae, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mae, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mae, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mae, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mae, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mae, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mae, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mae, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAE Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MPE.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmpe, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmpe, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEmpe, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mpe, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mpe, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mpe, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mpe, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mpe, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mpe, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mpe, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MPE Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAPE.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmape, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmape, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEmape, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mape, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mape, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mape, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mape, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mape, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mape, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mape, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAPE Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MSE_X1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmse_x1, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmse_x1, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mse_x1, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mse_x1, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mse_x1, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mse_x1, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mse_x1, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mse_x1, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mse_x1, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mse_x1, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MSE_X1 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAE_X1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmae_x1, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmae_x1, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mae_x1, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mae_x1, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mae_x1, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mae_x1, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mae_x1, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mae_x1, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mae_x1, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mae_x1, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAE_X1 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MPE_X1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmpe_x1, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmpe_x1, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mpe_x1, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mpe_x1, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mpe_x1, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mpe_x1, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mpe_x1, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mpe_x1, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mpe_x1, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mpe_x1, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MPE_X1 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAPE_X1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmape_x1, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmape_x1, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mape_x1, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mape_x1, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mape_x1, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mape_x1, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mape_x1, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mape_x1, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mape_x1, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mape_x1, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAPE_X1 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MSE_Trend.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(1, 2), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(FEsT_mse_trend, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mse_trend, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MSE_Trend Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAE_Trend.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(1, 2), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(FEsT_mae_trend, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mae_trend, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAE_Trend Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MPE_Trend.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(1, 2), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(FEsT_mpe_trend, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mpe_trend, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MPE_Trend Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAPE_Trend.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(1, 2), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(FEsT_mape_trend, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mape_trend, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAPE_Trend Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MSE_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmse_x2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmse_x2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mse_x2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mse_x2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mse_x2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mse_x2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mse_x2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mse_x2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mse_x2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mse_x2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MSE_x2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAE_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmae_x2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmae_x2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mae_x2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mae_x2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mae_x2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mae_x2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mae_x2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mae_x2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mae_x2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mae_x2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAE_x2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MPE_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmpe_x2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmpe_x2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mpe_x2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mpe_x2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mpe_x2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mpe_x2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mpe_x2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mpe_x2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mpe_x2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mpe_x2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MPE_x2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAPE_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmape_x2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmape_x2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mape_x2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mape_x2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mape_x2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mape_x2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mape_x2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mape_x2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mape_x2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mape_x2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAPE_x2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MSE_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmse_x3, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmse_x3, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mse_x3, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mse_x3, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mse_x3, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mse_x3, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mse_x3, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mse_x3, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mse_x3, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mse_x3, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MSE_x3 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAE_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmae_x3, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmae_x3, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mae_x3, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mae_x3, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mae_x3, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mae_x3, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mae_x3, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mae_x3, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mae_x3, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mae_x3, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAE_x3 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MPE_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmpe_x3, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmpe_x3, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mpe_x3, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mpe_x3, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mpe_x3, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mpe_x3, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mpe_x3, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mpe_x3, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mpe_x3, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mpe_x3, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MPE_x3 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/MAPE_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REmape_x3, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEmape_x3, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(ME_mape_x3, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2mape_x3, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_mape_x3, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_mape_x3, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_mape_x3, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_mape_x3, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_mape_x3, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_mape_x3, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("MAPE_x3 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_Trend.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(1, 2), mgp = c(2, 0.3, 0),cex=2)
hist(FEsT_bTrend, col = "green", main = "", xlab = "Study-level FE (with Trend)", 
     ylab = "Frequency", border = "black")
abline(v = timeHet*numYears/2, col = 'black', lwd = 3) #lty = 'dashed'
hist(FElT_bTrend, col = "red", main = "", xlab = "Location-level FE (with Trend)", 
     ylab = "Frequency", border = "black")
abline(v = timeHet*numYears/2, col = 'black', lwd = 3) #lty = 'dashed'
mtext("Trend Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_Intercept.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb0, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEb0, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb0, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2b0, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_b0, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_b0, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_b0, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_b0, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_b0, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_b0, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("Intercept Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# TRUE sig2 specified in vertical lines
# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/Coeff_sig2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REsig2, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEsig2, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEsig2, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FE2sig2, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(REl_sig2, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEl_sig2, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(MEl_sig2, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FEt_sig2, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
###
hist(FEsT_sig2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
hist(FElT_sig2, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
abline(v = sig2True, col = 'black', lwd = 3) #lty = 'dashed' 
mtext("SigmaHat Squared Coefficient Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/SE_Intercept.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb0SE, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEb0SE, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb0SE, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2b0SE, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_b0SE, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_b0SE, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_b0SE, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_b0SE, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_b0SE, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_b0SE, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("Intercept Standard Error Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/SE_x1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb1SE, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEb1SE, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb1SE, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2b1SE, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_b1SE, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_b1SE, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_b1SE, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_b1SE, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_b1SE, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_b1SE, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("x1 Standard Error Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/SE_x2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb2SE, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEb2SE, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb2SE, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2b2SE, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_b2SE, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_b2SE, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_b2SE, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_b2SE, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_b2SE, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_b2SE, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("x2 Standard Error Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/SE_x3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REb3SE, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEb3SE, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb3SE, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2b3SE, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_b3SE, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_b3SE, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_b3SE, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_b3SE, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_b3SE, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_b3SE, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("x3 Standard Error Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/R2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 3), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REr2, col = "yellow",main = "", xlab = "Random Effects",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FEr2, col = "green", main = "", xlab = "Study Fixed Effects",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FE2r2, col = "cyan", main = "", xlab = "Two-way Fixed Effects",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(REl_r2, col = "yellow",main = "", xlab = "Location Random Effects",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FEl_r2, col = "green", main = "", xlab = "Location Fixed Effects",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FEt_r2, col = "cyan", main = "", xlab = "Time Fixed Effects",
     ylab = "Frequency", border = "black")
###
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FEsT_r2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
hist(FElT_r2, col = "green", main = "", xlab = "Location FE (w/ Trend)",
     ylab = "Frequency", border = "black")
abline(v = R2True, col = 'black', lwd = 3) #lty = 'dashed'
mtext("R^2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/adjR2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 3), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REa_r2, col = "yellow",main = "", xlab = "Random Effects",
     ylab = "Frequency", border = "black")
hist(FEa_r2, col = "green", main = "", xlab = "Study Fixed Effects",
     ylab = "Frequency", border = "black")
hist(FE2a_r2, col = "cyan", main = "", xlab = "Two-way Fixed Effects",
     ylab = "Frequency", border = "black")
hist(REl_a_r2, col = "yellow",main = "", xlab = "Location Random Effects",
     ylab = "Frequency", border = "black")
hist(FEl_a_r2, col = "green", main = "", xlab = "Location Fixed Effects",
     ylab = "Frequency", border = "black")
hist(FEt_a_r2, col = "cyan", main = "", xlab = "Time Fixed Effects",
     ylab = "Frequency", border = "black")
###
hist(FEsT_a_r2, col = "yellow", main = "", xlab = "Study FE (w/ Trend)",
     ylab = "Frequency", border = "black")
hist(FElT_a_r2, col = "green", main = "", xlab = "Location FE (w/ Trend)",
     ylab = "Frequency", border = "black")
mtext("Adjusted R^2 Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()

###

# 1. Open jpeg file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/aic.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4), mgp = c(2, 0.3, 0),cex=2)
#par(mfrow = c(2, 3), mgp = c(2, 0.3, 0),cex=2)
hist(REaic, col = "yellow",main = "", xlab = "Study Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEaic, col = "green", main = "", xlab = "Study Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEb0, col = "cyan", main = "", xlab = "Study Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FE2aic, col = "red", main = "", xlab = "Two-way Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(REl_aic, col = "yellow",main = "", xlab = "Location Random Effects", 
     ylab = "Frequency", border = "black")
hist(FEl_aic, col = "green", main = "", xlab = "Location Fixed Effects", 
     ylab = "Frequency", border = "black")
hist(MEl_aic, col = "cyan", main = "", xlab = "Location Mixed Effects", 
     ylab = "Frequency", border = "black")
hist(FEt_aic, col = "red", main = "", xlab = "Time Fixed Effects", 
     ylab = "Frequency", border = "black")
###
hist(FEsT_aic, col = "yellow", main = "", xlab = "Study FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
hist(FElT_aic, col = "green", main = "", xlab = "Location FE (w/ Trend)", 
     ylab = "Frequency", border = "black")
mtext("AIC Distributions", side = 3, line = -2, outer = TRUE, cex=2)
# 3. Close the file
dev.off()



### Power Calculation (depending on magnitude of discrepancy from H0)) ###
library(stats)
library(exceedProb) # for pnct()

power_curve <- function(trueCoefs, SE, discrep_range, discrep_interval, df, alphaLevel, title, yTitle) {
  ### For Debugging ####
  # trueCoefs<-trueBetas[1]
  # SE<-REb1SE
  # discrep_range <- 100
  # discrep_interval <- 0.1
  # df <- REresult$df.residual
  
  power <- numeric(discrep_range)
  
  for (i in 1:(discrep_range)) {
    beta_alt <- trueCoefs + discrep_interval * i
      # B_alt = B - 0.1*i
    
    # NOTE: Nested loops are introduced since all SEs needs to be included, for each true parameter
    delta_beta <- numeric(length(SE))
    typeII <- numeric(length(SE))
    pow <- numeric(length(SE))
    
    # Find the area under the non-central t distribution to the left of alpha (type 2 error probability)
    crit_value <- qt(1-alphaLevel/2, df) #times 0.5 for the two-sided test
    
    for (j in 1:length(SE)){
      # Calculate non-centrality parameter- note that numerator could be replaced with discrepency!
      delta_beta[j] <- abs(trueCoefs - beta_alt) / SE[j]
        # delta = |B-B_alt|/SE
      
      typeII[j] <- pnct(crit_value, df, delta_beta[j])
      # Power = 1 Type II error probability
      pow[j] <- 1 - typeII[j]
    }
    
    power[i] <- mean(pow)
  }
  
  # Now graph the power
  x <- seq(0, discrep_interval * (discrep_range - 1), by = discrep_interval)
  
  # Create an empty plot (to be filled)
    # If you want to zoom in, change xlim to something smaller (like -1,1)
  plot(x, power, type = "n", lty = 2, col = "red", xlim=c(-10,10) , xlab = "Discrepency", ylab = paste("Power of",yTitle), main = title)
  # Add the 2 lines now
  lines(x, power, type = "l", lty = 2, col = "red")
  lines(-1*x, power, type = "l", lty = 2, col = "red")
  legend("bottomright", legend = c("Power"), col = c("red"), lty = 2)
}

### Calculate power for beta1
# 1. Open image file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/powerX1.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4),cex.axis=2,cex.lab=1.5,cex.title=1.25,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
#par(mfrow = c(3, 3),cex.axis=2,cex.lab=1.5,cex.title=1.25)
power_curve(trueBetas[1], REb1SE, 100, 0.1, RE$df.residual, alphaLevel, "Study Random Effects", "X1")
power_curve(trueBetas[1], FEb1SE, 100, 0.1, FE$df.residual, alphaLevel, "Study Fixed Effects", "X1")
power_curve(trueBetas[1], MEb1SE, 100, 0.1, length(y)-4, alphaLevel, "Study Mixed Effects", "X1")
power_curve(trueBetas[1], FE2b1SE, 100, 0.1, FE2$df.residual, alphaLevel, "Two-way Fixed Effects", "X1")
###
power_curve(trueBetas[1], REl_b1SE, 100, 0.1, REl$df.residual, alphaLevel, "Location Random Effects", "X1")
power_curve(trueBetas[1], FEl_b1SE, 100, 0.1, FEl$df.residual, alphaLevel, "Location Fixed Effects", "X1")
power_curve(trueBetas[1], MEl_b1SE, 100, 0.1, length(y)-4, alphaLevel, "Location Mixed Effects", "X1")
power_curve(trueBetas[1], FEt_b1SE, 100, 0.1, FEt$df.residual, alphaLevel, "Time Fixed Effects", "X1")
###
power_curve(trueBetas[1], FEsT_b1SE, 100, 0.1, FEsT$df.residual, alphaLevel, "Study FE (w/ Trend)", "X1")
power_curve(trueBetas[1], FElT_b1SE, 100, 0.1, FElT$df.residual, alphaLevel, "Location FE (w/ Trend)", "X1")
#mtext("Power of X1", side = 3, line = -2, outer = TRUE, cex=1.5)
# 3. Close the file
dev.off()

### Calculate power for beta2
# 1. Open image file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/powerX2.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4),cex.axis=2,cex.lab=1.5,cex.title=1.25,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
#par(mfrow = c(3, 3),cex.axis=2,cex.lab=1.5,cex.title=1.25)
power_curve(trueBetas[2], REb2SE, 100, 0.1, RE$df.residual, alphaLevel, "Study Random Effects", "X2")
power_curve(trueBetas[2], FEb2SE, 100, 0.1, FE$df.residual, alphaLevel, "Study Fixed Effects", "X2")
power_curve(trueBetas[2], MEb2SE, 100, 0.1, length(y)-4, alphaLevel, "Study Mixed Effects", "X2")
power_curve(trueBetas[2], FE2b2SE, 100, 0.1, FE2$df.residual, alphaLevel, "Two-way Fixed Effects", "X2")
###
power_curve(trueBetas[2], REl_b2SE, 100, 0.1, REl$df.residual, alphaLevel, "Location Random Effects", "X2")
power_curve(trueBetas[2], FEl_b2SE, 100, 0.1, FEl$df.residual, alphaLevel, "Location Fixed Effects", "X2")
power_curve(trueBetas[2], MEl_b2SE, 100, 0.1, length(y)-4, alphaLevel, "Location Mixed Effects", "X2")
power_curve(trueBetas[2], FEt_b2SE, 100, 0.1, FEt$df.residual, alphaLevel, "Time Fixed Effects", "X2")
###
power_curve(trueBetas[2], FEsT_b2SE, 100, 0.1, FEsT$df.residual, alphaLevel, "Study FE (w/ Trend)", "X2")
power_curve(trueBetas[2], FElT_b2SE, 100, 0.1, FElT$df.residual, alphaLevel, "Location FE (w/ Trend)", "X2")
#mtext("Power of X1", side = 3, line = -2, outer = TRUE, cex=1.5)
# 3. Close the file
dev.off()

### Calculate power for beta3
# 1. Open image file
png(filename="G:/SYNC/School/VT/RESEARCH/Dissertation/CH1 - 3YP/Summer 3rd Year Project - Shared Files/Code/Simulation/Simulation Results - R Files/x3/powerX3.png", width = 1200, height = 1600)
# 2. Create the plot
par(mfrow = c(3, 4),cex.axis=2,cex.lab=1.5,cex.title=1.25,mar=c(5,4+2,4,2)+0.1) # mar gives more padding in the left margin
#par(mfrow = c(3, 3),cex.axis=2,cex.lab=1.5,cex.title=1.25)
power_curve(trueBetas[3], REb3SE, 100, 0.1, RE$df.residual, alphaLevel, "Study Random Effects", "X3")
power_curve(trueBetas[3], FEb3SE, 100, 0.1, FE$df.residual, alphaLevel, "Study Fixed Effects", "X3")
power_curve(trueBetas[3], MEb3SE, 100, 0.1, length(y)-4, alphaLevel, "Study Mixed Effects", "X3")
power_curve(trueBetas[3], FE2b3SE, 100, 0.1, FE2$df.residual, alphaLevel, "Two-way Fixed Effects", "X3")
###
power_curve(trueBetas[3], REl_b3SE, 100, 0.1, REl$df.residual, alphaLevel, "Location Random Effects", "X3")
power_curve(trueBetas[3], FEl_b3SE, 100, 0.1, FEl$df.residual, alphaLevel, "Location Fixed Effects", "X3")
power_curve(trueBetas[3], MEl_b3SE, 100, 0.1, length(y)-4, alphaLevel, "Location Mixed Effects", "X3")
power_curve(trueBetas[3], FEt_b3SE, 100, 0.1, FEt$df.residual, alphaLevel, "Time Fixed Effects", "X3")
###
power_curve(trueBetas[3], FEsT_b3SE, 100, 0.1, FEsT$df.residual, alphaLevel, "Study FE (w/ Trend)", "X3")
power_curve(trueBetas[3], FElT_b3SE, 100, 0.1, FElT$df.residual, alphaLevel, "Location FE (w/ Trend)", "X3")
#mtext("Power of X1", side = 3, line = -2, outer = TRUE, cex=1.5)
# 3. Close the file
dev.off()



#%%###########################################################################
###                           EXPORT - Save ALL Data                       ###
##############################################################################

# RE
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(RE$model$y,RE$model$y-RE$residuals,RE$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(REmse,REmae,REmpe,REmape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REb0,REb1,REb2,REb3,REb0SE,REb1SE,REb2SE,REb3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(REaic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_REb0,CIhi_REb0,CIlo_REb1,CIhi_REb1,CIlo_REb2,CIhi_REb2,CIlo_REb3,CIhi_REb3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(RE_i2,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REmse_x1,REmae_x1,REmpe_x1,REmape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REmse_x2,REmae_x2,REmpe_x2,REmape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REmse_x3,REmae_x3,REmpe_x3,REmape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FE
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE$model$y,FE$fitted.values,FE$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FEmse,FEmae,FEmpe,FEmape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEb0,FEb1,FEb2,FEb3,FEb0SE,FEb1SE,FEb2SE,FEb3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FEaic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FEb0,CIhi_FEb0,CIlo_FEb1,CIhi_FEb1,CIlo_FEb2,CIhi_FEb2,CIlo_FEb3,CIhi_FEb3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE_i2,FE_Q,FE_dfOfGroups,FE_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEmse_x1,FEmae_x1,FEmpe_x1,FEmape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEmse_x2,FEmae_x2,FEmpe_x2,FEmape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEmse_x3,FEmae_x3,FEmpe_x3,FEmape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FE2
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2$model$y,FE2$fitted.values,FE2$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FE2mse,FE2mae,FE2mpe,FE2mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2b0,FE2b1,FE2b2,FE2b3,FE2b0SE,FE2b1SE,FE2b2SE,FE2b3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FE2aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FE2b0,CIhi_FE2b0,CIlo_FE2b1,CIhi_FE2b1,CIlo_FE2b2,CIhi_FE2b2,CIlo_FE2b3,CIhi_FE2b3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2_i2,FE2_Q,FE2_dfOfGroups,FE2_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2mse_x1,FE2mae_x1,FE2mpe_x1,FE2mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2mse_x2,FE2mae_x2,FE2mpe_x2,FE2mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FE2mse_x3,FE2mae_x3,FE2mpe_x3,FE2mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# REl
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REl$model$y,REl$model$y-REl$residuals,REl$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(REl_mse,REl_mae,REl_mpe,REl_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REl_b0,REl_b1,REl_b2,REl_b3,REl_b0SE,REl_b1SE,REl_b2SE,REl_b3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(REl_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_REl_b0,CIhi_REl_b0,CIlo_REl_b1,CIhi_REl_b1,CIlo_REl_b2,CIhi_REl_b2,CIlo_REl_b3,CIhi_REl_b3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(REl_i2,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REl_mse_x1,REl_mae_x1,REl_mpe_x1,REl_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REl_mse_x2,REl_mae_x2,REl_mpe_x2,REl_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(REl_mse_x3,REl_mae_x3,REl_mpe_x3,REl_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FEl
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl$model$y,FEl$fitted.values,FEl$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FEl_mse,FEl_mae,FEl_mpe,FEl_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl_b0,FEl_b1,FEl_b2,FEl_b3,FEl_b0SE,FEl_b1SE,FEl_b2SE,FEl_b3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FEl_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FEl_b0,CIhi_FEl_b0,CIlo_FEl_b1,CIhi_FEl_b1,CIlo_FEl_b2,CIhi_FEl_b2,CIlo_FEl_b3,CIhi_FEl_b3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl_i2,FEl_Q,FEl_dfOfGroups,FEl_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl_mse_x1,FEl_mae_x1,FEl_mpe_x1,FEl_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl_mse_x2,FEl_mae_x2,FEl_mpe_x2,FEl_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEl_mse_x3,FEl_mae_x3,FEl_mpe_x3,FEl_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FEt
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt$model$y,FEt$fitted.values,FEt$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FEt_mse,FEt_mae,FEt_mpe,FEt_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt_b0,FEt_b1,FEt_b2,FEt_b3,FEt_b0SE,FEt_b1SE,FEt_b2SE,FEt_b3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FEt_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FEt_b0,CIhi_FEt_b0,CIlo_FEt_b1,CIhi_FEt_b1,CIlo_FEt_b2,CIhi_FEt_b2,CIlo_FEt_b3,CIhi_FEt_b3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt_i2,FEt_Q,FEt_dfOfGroups,FEt_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt_mse_x1,FEt_mae_x1,FEt_mpe_x1,FEt_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt_mse_x2,FEt_mae_x2,FEt_mpe_x2,FEt_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEt_mse_x3,FEt_mae_x3,FEt_mpe_x3,FEt_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FEst
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT$model$y,FEsT$fitted.values,FEsT$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FEsT_mse,FEsT_mae,FEsT_mpe,FEsT_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_b0,FEsT_b1,FEsT_b2,FEsT_b3,FEsT_bTrend,FEsT_b0SE,FEsT_b1SE,FEsT_b2SE,FEsT_b3SE,FEsT_bTrendSE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FEsT_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FEsT_b0,CIhi_FEsT_b0,CIlo_FEsT_b1,CIhi_FEsT_b1,CIlo_FEsT_b2,CIhi_FEsT_b2,CIlo_FEsT_b3,CIhi_FEsT_b3,CIlo_FEsT_bTrend,CIhi_FEsT_bTrend),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_i2,FEsT_Q,FEsT_dfOfGroups,FEsT_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_mse_x1,FEsT_mae_x1,FEsT_mpe_x1,FEsT_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_mse_trend,FEsT_mae_trend,FEsT_mpe_trend,FEsT_mape_trend),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_trend.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_mse_x2,FEsT_mae_x2,FEsT_mpe_x2,FEsT_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FEsT_mse_x3,FEsT_mae_x3,FEsT_mpe_x3,FEsT_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# FElt
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT$model$y,FElT$fitted.values,FElT$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(FElT_mse,FElT_mae,FElT_mpe,FElT_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_b0,FElT_b1,FElT_b2,FElT_b3,FElT_bTrend,FElT_b0SE,FElT_b1SE,FElT_b2SE,FElT_b3SE,FElT_bTrendSE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(FElT_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_FElT_b0,CIhi_FElT_b0,CIlo_FElT_b1,CIhi_FElT_b1,CIlo_FElT_b2,CIhi_FElT_b2,CIlo_FElT_b3,CIhi_FElT_b3,CIlo_FElT_bTrend,CIhi_FElT_bTrend),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_i2,FElT_Q,FElT_dfOfGroups,FElT_h2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_mse_x1,FElT_mae_x1,FElT_mpe_x1,FElT_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_mse_trend,FElT_mae_trend,FElT_mpe_trend,FElT_mape_trend),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_trend.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_mse_x2,FElT_mae_x2,FElT_mpe_x2,FElT_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(FElT_mse_x3,FElT_mae_x3,FElT_mpe_x3,FElT_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# ME
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(y,y-summary(ME)$residuals,summary(ME)$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(MEmse,MEmae,MEmpe,MEmape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(MEb0,MEb1,MEb2,MEb3,MEb0SE,MEb1SE,MEb2SE,MEb3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(ME_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_MEb0,CIhi_MEb0,CIlo_MEb1,CIhi_MEb1,CIlo_MEb2,CIhi_MEb2,CIlo_MEb3,CIhi_MEb3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(ME_i2,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(ME_mse_x1,ME_mae_x1,ME_mpe_x1,ME_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(ME_mse_x2,ME_mae_x2,ME_mpe_x2,ME_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(ME_mse_x3,ME_mae_x3,ME_mpe_x3,ME_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")

# MEl
write.table(X,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_X.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(y,y-summary(MEl)$residuals,summary(MEl)$residuals),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_y(Hat)&Resid.txt", col.names = TRUE,row.names=FALSE, sep=",")
write.table(cbind(MEl_mse,MEl_mae,MEl_mpe,MEl_mape),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(MEl_b0,MEl_b1,MEl_b2,MEl_b3,MEl_b0SE,MEl_b1SE,MEl_b2SE,MEl_b3SE),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_betas&SEs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(MEl_aic,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_aic.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(CIlo_MEl_b0,CIhi_MEl_b0,CIlo_MEl_b1,CIhi_MEl_b1,CIlo_MEl_b2,CIhi_MEl_b2,CIlo_MEl_b3,CIhi_MEl_b3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_CIs.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(MEl_i2,file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_i2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(MEl_mse_x1,MEl_mae_x1,MEl_mpe_x1,MEl_mape_x1),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x1.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(MEl_mse_x2,MEl_mae_x2,MEl_mpe_x2,MEl_mape_x2),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x2.txt", col.names = TRUE, row.names=FALSE, sep=",")
write.table(cbind(MEl_mse_x3,MEl_mae_x3,MEl_mpe_x3,MEl_mape_x3),file="G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x3.txt", col.names = TRUE, row.names=FALSE, sep=",")



### GIVE COLUMN HEADERS NAMES (to exported files) ###

### RE ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat","resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--I^2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_i2.txt")
# Add a header to the data frame
header <- c("I^2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\RE_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)


### FE ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### FE2 ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FE2_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### REl ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat","resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--I^2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_i2.txt")
# Add a header to the data frame
header <- c("I^2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\REl_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### FEl ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEl_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### FEt ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEt_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### FEst ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","bTrend","b0SE","b1SE","b2SE","b3SE","bTrendSE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB","bTrend_LB","bTrend_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEst_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (trend)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_trend.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_trend.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (trend)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_trend.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_trend.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FEsT_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### FElt ###
# y/resid & MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
# MSE
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","bTrend","b0SE","b1SE","b2SE","b3SE","bTrendSE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB","bTrend_LB","bTrend_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElt_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--i2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_i2.txt")
# Add a header to the data frame
header <- c("i2","Q","dfOfGroups","h2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (trend)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_trend.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_trend.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (trend)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_trend.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_trend.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\FElT_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### ME ###
#--y & resid--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
#--MSE--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--I^2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_i2.txt")
# Add a header to the data frame
header <- c("I^2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\ME_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

### MEl ###
#--y & resid--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_y(Hat)&Resid.txt")
# Add a header to the data frame
header <- c("y", "yHat", "resid")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_y(Hat)&Resid.txt", sep = ",", row.names = FALSE, col.names = TRUE, append = FALSE)
#--MSE--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--betas & SEs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_betas&SEs.txt")
# Add a header to the data frame
header <- c("b0","b1","b2","b3","b0SE","b1SE","b2SE","b3SE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_betas&SEs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--AIC--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_aic.txt")
# Add a header to the data frame
header <- c("AIC")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_aic.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--CIs--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_CIs.txt")
# Add a header to the data frame
header <- c("b0_LB","b0_UB","b1_LB","b1_UB","b2_LB","b2_UB","b3_LB","b3_UB")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_CIs.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#--I^2--
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_i2.txt")
# Add a header to the data frame
header <- c("I^2")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_i2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x1)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x1.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x1.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x2)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x2.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x2.txt", sep = ",", row.names = FALSE, col.names = TRUE)
#---MSE (x3)---
# Read the CSV file into a data frame
df <- read.csv("G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x3.txt")
# Add a header to the data frame
header <- c("MSE","MAE","MPE","MAPE")
names(df) <- header
# Write the data frame to the CSV file with header
write.csv(df, "G:\\SYNC\\School\\VT\\RESEARCH\\Dissertation\\CH1 - 3YP\\Summer 3rd Year Project - Shared Files\\Code\\Simulation\\Simulation Results - Data\\x3\\MEl_mse_x3.txt", sep = ",", row.names = FALSE, col.names = TRUE)

# Audio to signal power & exports completed
beep(sound = 5)
