
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
#n <- 50
#
n <- 100
#n <- 150

case <- 3 # Options: 1-12

model1 <- TRUE # RE s
model2 <- TRUE # RE l
model3 <- TRUE # FE s
model4 <- TRUE # FE l
model5 <- TRUE # FE t
model6 <- TRUE # FE lt
model7 <- FALSE # FE s (with Trend)
model8 <- TRUE # FE l (with Trend)

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
  else if(n==100){      # N = 2,500~
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
  else if(n==100){      # N = 7,500~
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# !!! EVERYTHING below is vastly updated code !!!
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(glmmTMB)
set.seed(42)

# All of these empty vectors will store various results for each model
  # (i.e. REs, FEs, FE_lt, ...)

# COEFFICIENT ESTIMATES

# beta_0
RE_beta0   <- numeric(N)
FE_beta0   <- numeric(N)
FE2_beta0  <- numeric(N)
REl_beta0  <- numeric(N)
FEl_beta0  <- numeric(N)
FEt_beta0  <- numeric(N)
FEsT_beta0 <- numeric(N)
FElT_beta0 <- numeric(N)

# beta_1 (x1)
RE_beta1   <- numeric(N)
FE_beta1   <- numeric(N)
FE2_beta1  <- numeric(N)
REl_beta1  <- numeric(N)
FEl_beta1  <- numeric(N)
FEt_beta1  <- numeric(N)
FEsT_beta1 <- numeric(N)
FElT_beta1 <- numeric(N)

# beta_2 (x2)
RE_beta2   <- numeric(N)
FE_beta2   <- numeric(N)
FE2_beta2  <- numeric(N)
REl_beta2  <- numeric(N)
FEl_beta2  <- numeric(N)
FEt_beta2  <- numeric(N)
FEsT_beta2 <- numeric(N)
FElT_beta2 <- numeric(N)

# beta_3 (x3)
RE_beta3   <- numeric(N)
FE_beta3   <- numeric(N)
FE2_beta3  <- numeric(N)
REl_beta3  <- numeric(N)
FEl_beta3  <- numeric(N)
FEt_beta3  <- numeric(N)
FEsT_beta3 <- numeric(N)
FElT_beta3 <- numeric(N)

# TREND COEFFICIENTS
FEsT_betaTrend <- numeric(N)
FElT_betaTrend <- numeric(N)

# RESIDUAL VARIANCE
RE_rep_resid_var   <- numeric(N)
FE_rep_resid_var   <- numeric(N)
FE2_rep_resid_var  <- numeric(N)
REl_rep_resid_var  <- numeric(N)
FEl_rep_resid_var  <- numeric(N)
FEt_rep_resid_var  <- numeric(N)
FEsT_rep_resid_var <- numeric(N)
FElT_rep_resid_var <- numeric(N)

# MODEL-BASED STANDARD ERRORS

# beta_0
RE_SE_beta0   <- numeric(N)
FE_SE_beta0   <- numeric(N)
FE2_SE_beta0  <- numeric(N)
REl_SE_beta0  <- numeric(N)
FEl_SE_beta0  <- numeric(N)
FEt_SE_beta0  <- numeric(N)
FEsT_SE_beta0 <- numeric(N)
FElT_SE_beta0 <- numeric(N)

# beta_1
RE_SE_beta1   <- numeric(N)
FE_SE_beta1   <- numeric(N)
FE2_SE_beta1  <- numeric(N)
REl_SE_beta1  <- numeric(N)
FEl_SE_beta1  <- numeric(N)
FEt_SE_beta1  <- numeric(N)
FEsT_SE_beta1 <- numeric(N)
FElT_SE_beta1 <- numeric(N)

# beta_2
RE_SE_beta2   <- numeric(N)
FE_SE_beta2   <- numeric(N)
FE2_SE_beta2  <- numeric(N)
REl_SE_beta2  <- numeric(N)
FEl_SE_beta2  <- numeric(N)
FEt_SE_beta2  <- numeric(N)
FEsT_SE_beta2 <- numeric(N)
FElT_SE_beta2 <- numeric(N)

# beta_3
RE_SE_beta3   <- numeric(N)
FE_SE_beta3   <- numeric(N)
FE2_SE_beta3  <- numeric(N)
REl_SE_beta3  <- numeric(N)
FEl_SE_beta3  <- numeric(N)
FEt_SE_beta3  <- numeric(N)
FEsT_SE_beta3 <- numeric(N)
FElT_SE_beta3 <- numeric(N)

# TREND STANDARD ERRORS
FEsT_SE_betaTrend <- numeric(N)
FElT_SE_betaTrend <- numeric(N)

# RESIDUAL PERFORMANCE WITHIN EACH REPLICATION

# MSE
RE_rep_resid_mse   <- numeric(N)
FE_rep_resid_mse   <- numeric(N)
FE2_rep_resid_mse  <- numeric(N)
REl_rep_resid_mse  <- numeric(N)
FEl_rep_resid_mse  <- numeric(N)
FEt_rep_resid_mse  <- numeric(N)
FEsT_rep_resid_mse <- numeric(N)
FElT_rep_resid_mse <- numeric(N)

# RMSE
RE_rep_resid_rmse   <- numeric(N)
FE_rep_resid_rmse   <- numeric(N)
FE2_rep_resid_rmse  <- numeric(N)
REl_rep_resid_rmse  <- numeric(N)
FEl_rep_resid_rmse  <- numeric(N)
FEt_rep_resid_rmse  <- numeric(N)
FEsT_rep_resid_rmse <- numeric(N)
FElT_rep_resid_rmse <- numeric(N)

# MAE
RE_rep_resid_mae   <- numeric(N)
FE_rep_resid_mae   <- numeric(N)
FE2_rep_resid_mae  <- numeric(N)
REl_rep_resid_mae  <- numeric(N)
FEl_rep_resid_mae  <- numeric(N)
FEt_rep_resid_mae  <- numeric(N)
FEsT_rep_resid_mae <- numeric(N)
FElT_rep_resid_mae <- numeric(N)

# MPE
RE_rep_resid_mpe   <- numeric(N)
FE_rep_resid_mpe   <- numeric(N)
FE2_rep_resid_mpe  <- numeric(N)
REl_rep_resid_mpe  <- numeric(N)
FEl_rep_resid_mpe  <- numeric(N)
FEt_rep_resid_mpe  <- numeric(N)
FEsT_rep_resid_mpe <- numeric(N)
FElT_rep_resid_mpe <- numeric(N)

# MAPE
RE_rep_resid_mape   <- numeric(N)
FE_rep_resid_mape   <- numeric(N)
FE2_rep_resid_mape  <- numeric(N)
REl_rep_resid_mape  <- numeric(N)
FEl_rep_resid_mape  <- numeric(N)
FEt_rep_resid_mape  <- numeric(N)
FEsT_rep_resid_mape <- numeric(N)
FElT_rep_resid_mape <- numeric(N)


# COEFFICIENT ERRORS WITHIN EACH REPLICATION

# ---------------- x1 ----------------
# Squared error
RE_rep_sqerr_x1   <- numeric(N)
FE_rep_sqerr_x1   <- numeric(N)
FE2_rep_sqerr_x1  <- numeric(N)
REl_rep_sqerr_x1  <- numeric(N)
FEl_rep_sqerr_x1  <- numeric(N)
FEt_rep_sqerr_x1  <- numeric(N)
FEsT_rep_sqerr_x1 <- numeric(N)
FElT_rep_sqerr_x1 <- numeric(N)

# Absolute error
RE_rep_abserr_x1   <- numeric(N)
FE_rep_abserr_x1   <- numeric(N)
FE2_rep_abserr_x1  <- numeric(N)
REl_rep_abserr_x1  <- numeric(N)
FEl_rep_abserr_x1  <- numeric(N)
FEt_rep_abserr_x1  <- numeric(N)
FEsT_rep_abserr_x1 <- numeric(N)
FElT_rep_abserr_x1 <- numeric(N)

# Percentage error
RE_rep_pcterr_x1   <- numeric(N)
FE_rep_pcterr_x1   <- numeric(N)
FE2_rep_pcterr_x1  <- numeric(N)
REl_rep_pcterr_x1  <- numeric(N)
FEl_rep_pcterr_x1  <- numeric(N)
FEt_rep_pcterr_x1  <- numeric(N)
FEsT_rep_pcterr_x1 <- numeric(N)
FElT_rep_pcterr_x1 <- numeric(N)

# Absolute percentage error
RE_rep_abspcterr_x1   <- numeric(N)
FE_rep_abspcterr_x1   <- numeric(N)
FE2_rep_abspcterr_x1  <- numeric(N)
REl_rep_abspcterr_x1  <- numeric(N)
FEl_rep_abspcterr_x1  <- numeric(N)
FEt_rep_abspcterr_x1  <- numeric(N)
FEsT_rep_abspcterr_x1 <- numeric(N)
FElT_rep_abspcterr_x1 <- numeric(N)

# ---------------- x2 ----------------
# Squared error
RE_rep_sqerr_x2   <- numeric(N)
FE_rep_sqerr_x2   <- numeric(N)
FE2_rep_sqerr_x2  <- numeric(N)
REl_rep_sqerr_x2  <- numeric(N)
FEl_rep_sqerr_x2  <- numeric(N)
FEt_rep_sqerr_x2  <- numeric(N)
FEsT_rep_sqerr_x2 <- numeric(N)
FElT_rep_sqerr_x2 <- numeric(N)

# Absolute error
RE_rep_abserr_x2   <- numeric(N)
FE_rep_abserr_x2   <- numeric(N)
FE2_rep_abserr_x2  <- numeric(N)
REl_rep_abserr_x2  <- numeric(N)
FEl_rep_abserr_x2  <- numeric(N)
FEt_rep_abserr_x2  <- numeric(N)
FEsT_rep_abserr_x2 <- numeric(N)
FElT_rep_abserr_x2 <- numeric(N)

# Percentage error
RE_rep_pcterr_x2   <- numeric(N)
FE_rep_pcterr_x2   <- numeric(N)
FE2_rep_pcterr_x2  <- numeric(N)
REl_rep_pcterr_x2  <- numeric(N)
FEl_rep_pcterr_x2  <- numeric(N)
FEt_rep_pcterr_x2  <- numeric(N)
FEsT_rep_pcterr_x2 <- numeric(N)
FElT_rep_pcterr_x2 <- numeric(N)

# Absolute percentage error
RE_rep_abspcterr_x2   <- numeric(N)
FE_rep_abspcterr_x2   <- numeric(N)
FE2_rep_abspcterr_x2  <- numeric(N)
REl_rep_abspcterr_x2  <- numeric(N)
FEl_rep_abspcterr_x2  <- numeric(N)
FEt_rep_abspcterr_x2  <- numeric(N)
FEsT_rep_abspcterr_x2 <- numeric(N)
FElT_rep_abspcterr_x2 <- numeric(N)

# ---------------- x3 ----------------
# Squared error
RE_rep_sqerr_x3   <- numeric(N)
FE_rep_sqerr_x3   <- numeric(N)
FE2_rep_sqerr_x3  <- numeric(N)
REl_rep_sqerr_x3  <- numeric(N)
FEl_rep_sqerr_x3  <- numeric(N)
FEt_rep_sqerr_x3  <- numeric(N)
FEsT_rep_sqerr_x3 <- numeric(N)
FElT_rep_sqerr_x3 <- numeric(N)

# Absolute error
RE_rep_abserr_x3   <- numeric(N)
FE_rep_abserr_x3   <- numeric(N)
FE2_rep_abserr_x3  <- numeric(N)
REl_rep_abserr_x3  <- numeric(N)
FEl_rep_abserr_x3  <- numeric(N)
FEt_rep_abserr_x3  <- numeric(N)
FEsT_rep_abserr_x3 <- numeric(N)
FElT_rep_abserr_x3 <- numeric(N)

# Percentage error
RE_rep_pcterr_x3   <- numeric(N)
FE_rep_pcterr_x3   <- numeric(N)
FE2_rep_pcterr_x3  <- numeric(N)
REl_rep_pcterr_x3  <- numeric(N)
FEl_rep_pcterr_x3  <- numeric(N)
FEt_rep_pcterr_x3  <- numeric(N)
FEsT_rep_pcterr_x3 <- numeric(N)
FElT_rep_pcterr_x3 <- numeric(N)

# Absolute percentage error
RE_rep_abspcterr_x3   <- numeric(N)
FE_rep_abspcterr_x3   <- numeric(N)
FE2_rep_abspcterr_x3  <- numeric(N)
REl_rep_abspcterr_x3  <- numeric(N)
FEl_rep_abspcterr_x3  <- numeric(N)
FEt_rep_abspcterr_x3  <- numeric(N)
FEsT_rep_abspcterr_x3 <- numeric(N)
FElT_rep_abspcterr_x3 <- numeric(N)
# --------------------------------

# TREND ERRORS WITHIN EACH REPLICATION
FEsT_rep_sqerr_trend    <- numeric(N)
FElT_rep_sqerr_trend    <- numeric(N)

FEsT_rep_abserr_trend   <- numeric(N)
FElT_rep_abserr_trend   <- numeric(N)

FEsT_rep_pcterr_trend   <- numeric(N)
FElT_rep_pcterr_trend   <- numeric(N)

FEsT_rep_abspcterr_trend <- numeric(N)
FElT_rep_abspcterr_trend <- numeric(N)


# 95% CONFIDENCE INTERVALS

#beta_0
RE_CIlo_beta0   <- numeric(N)
FE_CIlo_beta0   <- numeric(N)
FE2_CIlo_beta0  <- numeric(N)
REl_CIlo_beta0  <- numeric(N)
FEl_CIlo_beta0  <- numeric(N)
FEt_CIlo_beta0  <- numeric(N)
FEsT_CIlo_beta0 <- numeric(N)
FElT_CIlo_beta0 <- numeric(N)

RE_CIhi_beta0   <- numeric(N)
FE_CIhi_beta0   <- numeric(N)
FE2_CIhi_beta0  <- numeric(N)
REl_CIhi_beta0  <- numeric(N)
FEl_CIhi_beta0  <- numeric(N)
FEt_CIhi_beta0  <- numeric(N)
FEsT_CIhi_beta0 <- numeric(N)
FElT_CIhi_beta0 <- numeric(N)

#beta_1
RE_CIlo_beta1   <- numeric(N)
FE_CIlo_beta1   <- numeric(N)
FE2_CIlo_beta1  <- numeric(N)
REl_CIlo_beta1  <- numeric(N)
FEl_CIlo_beta1  <- numeric(N)
FEt_CIlo_beta1  <- numeric(N)
FEsT_CIlo_beta1 <- numeric(N)
FElT_CIlo_beta1 <- numeric(N)

RE_CIhi_beta1   <- numeric(N)
FE_CIhi_beta1   <- numeric(N)
FE2_CIhi_beta1  <- numeric(N)
REl_CIhi_beta1  <- numeric(N)
FEl_CIhi_beta1  <- numeric(N)
FEt_CIhi_beta1  <- numeric(N)
FEsT_CIhi_beta1 <- numeric(N)
FElT_CIhi_beta1 <- numeric(N)

#beta_2
RE_CIlo_beta2   <- numeric(N)
FE_CIlo_beta2   <- numeric(N)
FE2_CIlo_beta2  <- numeric(N)
REl_CIlo_beta2  <- numeric(N)
FEl_CIlo_beta2  <- numeric(N)
FEt_CIlo_beta2  <- numeric(N)
FEsT_CIlo_beta2 <- numeric(N)
FElT_CIlo_beta2 <- numeric(N)

RE_CIhi_beta2   <- numeric(N)
FE_CIhi_beta2   <- numeric(N)
FE2_CIhi_beta2  <- numeric(N)
REl_CIhi_beta2  <- numeric(N)
FEl_CIhi_beta2  <- numeric(N)
FEt_CIhi_beta2  <- numeric(N)
FEsT_CIhi_beta2 <- numeric(N)
FElT_CIhi_beta2 <- numeric(N)

#beta_3
RE_CIlo_beta3   <- numeric(N)
FE_CIlo_beta3   <- numeric(N)
FE2_CIlo_beta3  <- numeric(N)
REl_CIlo_beta3  <- numeric(N)
FEl_CIlo_beta3  <- numeric(N)
FEt_CIlo_beta3  <- numeric(N)
FEsT_CIlo_beta3 <- numeric(N)
FElT_CIlo_beta3 <- numeric(N)

RE_CIhi_beta3   <- numeric(N)
FE_CIhi_beta3   <- numeric(N)
FE2_CIhi_beta3  <- numeric(N)
REl_CIhi_beta3  <- numeric(N)
FEl_CIhi_beta3  <- numeric(N)
FEt_CIhi_beta3  <- numeric(N)
FEsT_CIhi_beta3 <- numeric(N)
FElT_CIhi_beta3 <- numeric(N)

#trend
FEsT_CIlo_betaTrend <- numeric(N)
FEsT_CIhi_betaTrend <- numeric(N)

FElT_CIlo_betaTrend <- numeric(N)
FElT_CIhi_betaTrend <- numeric(N)



# DEFINE LOCATION-SPECIFIC MEANS (For data generation)
if(numCountries == 5 & spreadMuY == -2){
  
  locationMeans <- c(
    -2, -1, 0, 1, 2
  )
  
} else if(numCountries == 5 & spreadMuY == -10){
  
  locationMeans <- c(
    -10, -5, 0, 5, 10
  )
  
} else if(numCountries == 9 & spreadMuY == -2){
  
  locationMeans <- c(
    -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2
  )
  
} else if(numCountries == 9 & spreadMuY == -10){
  
  locationMeans <- c(
    -10, -7.5, -5, -2.5, 0, 2.5, 5, 7.5, 10
  )
  
} else if(numCountries == 15 & spreadMuY == -2){
  
  locationMeans <- c(
    -2,
    -1.66, -1.33, -1, -0.66, -0.33, -0.15, 0,
    0.15, 0.33, 0.66, 1, 1.33, 1.66, 2
  )
  
} else if(numCountries == 15 & spreadMuY == -10){
  
  locationMeans <- c(
    -10,
    -8.5, -7, -5.5, -4, -2.5, -1, 0,
    1, 2.5, 4, 5.5, 7, 8.5, 10
  )
}



start<-Sys.time()
for (k in 1:N) {
  # To show which iteration you're on
  cat(k,"/",N,"\n")
  
  #################
  ### SETUP     ###
  #################

  # NOW GENERATE THE DATA

  X <- NULL
  y <- numeric(0)
  
  paper_iterator <- 1
  
  for(countryIndex in 1:numCountries){
    
    for(j in 0:(numYears - 1)){
      
      # Location + time heterogeneity
      mu_y <- locationMeans[countryIndex] + (timeHet * j)
      
      mu_0 <- c(
        mu_y,
        mu_x1,
        mu_x2,
        mu_x3
      )
      
      # Generate n observations for this study
      z_data <- rmvnorm(
        n,
        mean = mu_0,
        sigma = SIGMA
      )
      
      y2 <- z_data[, 1]
      x  <- z_data[, 2:4]
      
      # Create study data
      X2 <- data.frame(
        constant = 1,
        x = x,
        year = j + 2020,
        country = countryIndex,
        paper = paper_iterator,
        trend = -1 + j/2
      )
      
      # Add study to full dataset
      X <- bind_rows(X, X2)
      y <- c(y, y2)
      
      # One new paper per country-year combination
      paper_iterator <- paper_iterator + 1
    }
  }
  


  ##################
  ### RUN MODELS ###
  ##################
  
  # To ensure NAs don't keep the regressions from running
  y <- na.omit(y)
  X <- na.omit(X)
  
  ### Random Effects (Study Level) ###
  # paper effects treated as random intercepts
  if(model1==TRUE){
    RE <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + (1 | paper),
      data = X,
      family = gaussian()#,weights = w
    )

    # Extract model objects
    REres  <- residuals(RE, type = "response")
    REy    <- model.response(model.frame(RE))
    REcoef <- fixef(RE)$cond
    REvcov <- vcov(RE)$cond
    REse   <- sqrt(diag(REvcov))
    
    # THEN
    
    # Residual performance (within replication k)
    RE_rep_resid_mse[k]  <- mean(REres^2)                 # Mean Squared Error
    RE_rep_resid_rmse[k] <- sqrt(mean(REres^2))           # Root Mean Squared Error
    RE_rep_resid_mae[k]  <- mean(abs(REres))              # Mean Absolute Error
    RE_rep_resid_mpe[k]  <- mean(REres / REy)             # Mean Percentage Error
    RE_rep_resid_mape[k] <- mean(abs(REres / REy))        # Mean Absolute Percentage Error
    
    
    # Coefficient estimates (within replication k)
    RE_beta0[k] <- REcoef["(Intercept)"]
    RE_beta1[k] <- REcoef["x.1"]
    RE_beta2[k] <- REcoef["x.2"]
    RE_beta3[k] <- REcoef["x.3"]
    
    
    # Coefficient errors within replication k (Used later to calculate Monte Carlo performance)    
    # beta_1
    RE_rep_sqerr_x1[k]     <- (RE_beta1[k] - b1True)^2
    RE_rep_abserr_x1[k]    <- abs(RE_beta1[k] - b1True)
    RE_rep_pcterr_x1[k]    <- (RE_beta1[k] - b1True) / b1True
    RE_rep_abspcterr_x1[k] <- abs(RE_rep_pcterr_x1[k])
    
    # beta_2
    RE_rep_sqerr_x2[k]     <- (RE_beta2[k] - b2True)^2
    RE_rep_abserr_x2[k]    <- abs(RE_beta2[k] - b2True)
    RE_rep_pcterr_x2[k]    <- (RE_beta2[k] - b2True) / b2True
    RE_rep_abspcterr_x2[k] <- abs(RE_rep_pcterr_x2[k])
    
    # beta_3
    RE_rep_sqerr_x3[k]     <- (RE_beta3[k] - b3True)^2
    RE_rep_abserr_x3[k]    <- abs(RE_beta3[k] - b3True)
    RE_rep_pcterr_x3[k]    <- (RE_beta3[k] - b3True) / b3True
    RE_rep_abspcterr_x3[k] <- abs(RE_rep_pcterr_x3[k])
    
    # Residual variance (within replication k)
    RE_rep_resid_var[k] <- sigma(RE)^2
    
    # Model-based standard errors (within replication k)
    RE_SE_beta0[k] <- REse["(Intercept)"]
    RE_SE_beta1[k] <- REse["x.1"]
    RE_SE_beta2[k] <- REse["x.2"]
    RE_SE_beta3[k] <- REse["x.3"]
    
    # 95% confidence intervals (within replication k)
    RE_CIlo_beta0[k] <- RE_beta0[k] - 1.96 * RE_SE_beta0[k]
    RE_CIhi_beta0[k] <- RE_beta0[k] + 1.96 * RE_SE_beta0[k]
    
    RE_CIlo_beta1[k] <- RE_beta1[k] - 1.96 * RE_SE_beta1[k]
    RE_CIhi_beta1[k] <- RE_beta1[k] + 1.96 * RE_SE_beta1[k]
    
    RE_CIlo_beta2[k] <- RE_beta2[k] - 1.96 * RE_SE_beta2[k]
    RE_CIhi_beta2[k] <- RE_beta2[k] + 1.96 * RE_SE_beta2[k]
    
    RE_CIlo_beta3[k] <- RE_beta3[k] - 1.96 * RE_SE_beta3[k]
    RE_CIhi_beta3[k] <- RE_beta3[k] + 1.96 * RE_SE_beta3[k]
  }
  
  ### Study-level Fixed Effects ###
  if(model3==TRUE){
    # FE: paper effects treated as fixed intercepts
    FE <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + factor(paper),
      data = X,
      family = gaussian()#,weights = w
    )
    
    # Extract model objects
    FEres  <- residuals(FE, type = "response")
    FEy    <- model.response(model.frame(FE))
    FEcoef <- fixef(FE)$cond
    FEvcov <- vcov(FE)$cond
    FEse   <- sqrt(diag(FEvcov))

    # THEN
    
    # Residual performance (within replication k)
    FE_rep_resid_mse[k]  <- mean(FEres^2)                 # Mean Squared Error
    FE_rep_resid_rmse[k] <- sqrt(mean(FEres^2))           # Root Mean Squared Error
    FE_rep_resid_mae[k]  <- mean(abs(FEres))               # Mean Absolute Error
    FE_rep_resid_mpe[k]  <- mean(FEres / FEy)              # Mean Percentage Error
    FE_rep_resid_mape[k] <- mean(abs(FEres / FEy))         # Mean Absolute Percentage Error

    # Coefficient estimates (within replication k)
    FE_beta0[k] <- FEcoef["(Intercept)"]
    FE_beta1[k] <- FEcoef["x.1"]
    FE_beta2[k] <- FEcoef["x.2"]
    FE_beta3[k] <- FEcoef["x.3"]

    # Coefficient errors within replication k (Used later to calculate Monte Carlo performance)
    # beta_1
    FE_rep_sqerr_x1[k]  <- (FE_beta1[k] - b1True)^2
    FE_rep_abserr_x1[k] <- abs(FE_beta1[k] - b1True)
    FE_rep_pcterr_x1[k] <- (FE_beta1[k] - b1True) / b1True
    FE_rep_abspcterr_x1[k] <- abs(FE_rep_pcterr_x1[k])
    
    # beta_2
    FE_rep_sqerr_x2[k]  <- (FE_beta2[k] - b2True)^2
    FE_rep_abserr_x2[k] <- abs(FE_beta2[k] - b2True)
    FE_rep_pcterr_x2[k] <- (FE_beta2[k] - b2True) / b2True
    FE_rep_abspcterr_x2[k] <- abs(FE_rep_pcterr_x2[k])
    
    # beta_3
    FE_rep_sqerr_x3[k]  <- (FE_beta3[k] - b3True)^2
    FE_rep_abserr_x3[k] <- abs(FE_beta3[k] - b3True)
    FE_rep_pcterr_x3[k] <- (FE_beta3[k] - b3True) / b3True
    FE_rep_abspcterr_x3[k] <- abs(FE_rep_pcterr_x3[k])

    # Residual variance (within replication k)
    FE_rep_resid_var[k] <- sigma(FE)^2

    # Model-based standard errors (within replication k)    
    FE_SE_beta0[k] <- FEse["(Intercept)"]
    FE_SE_beta1[k] <- FEse["x.1"]
    FE_SE_beta2[k] <- FEse["x.2"]
    FE_SE_beta3[k] <- FEse["x.3"]

    # 95% confidence intervals (within replication k)
    FE_CIlo_beta0[k] <- FE_beta0[k] - 1.96 * FE_SE_beta0[k]
    FE_CIhi_beta0[k] <- FE_beta0[k] + 1.96 * FE_SE_beta0[k]
    
    FE_CIlo_beta1[k] <- FE_beta1[k] - 1.96 * FE_SE_beta1[k]
    FE_CIhi_beta1[k] <- FE_beta1[k] + 1.96 * FE_SE_beta1[k]
    
    FE_CIlo_beta2[k] <- FE_beta2[k] - 1.96 * FE_SE_beta2[k]
    FE_CIhi_beta2[k] <- FE_beta2[k] + 1.96 * FE_SE_beta2[k]
    
    FE_CIlo_beta3[k] <- FE_beta3[k] - 1.96 * FE_SE_beta3[k]
    FE_CIhi_beta3[k] <- FE_beta3[k] + 1.96 * FE_SE_beta3[k]
  }
  
  ### Two Way Fixed Effects ###
  if(model6==TRUE){
    # FE: paper effects treated as fixed intercepts
    FE2 <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + factor(country) + factor(year),
      data = X,
      family = gaussian()#,weights = w
    )
    
    # Extract model objects
    FE2res  <- residuals(FE2, type = "response")
    FE2y    <- model.response(model.frame(FE2))
    FE2coef <- fixef(FE2)$cond
    FE2vcov <- vcov(FE2)$cond
    FE2se   <- sqrt(diag(FE2vcov))
    
    # Residual performance (within replication k)
    FE2_rep_resid_mse[k]  <- mean(FE2res^2)
    FE2_rep_resid_rmse[k] <- sqrt(mean(FE2res^2))
    FE2_rep_resid_mae[k]  <- mean(abs(FE2res))
    FE2_rep_resid_mpe[k]  <- mean(FE2res / FE2y)
    FE2_rep_resid_mape[k] <- mean(abs(FE2res / FE2y))
    
    # Coefficient estimates
    FE2_beta0[k] <- FE2coef["(Intercept)"]
    FE2_beta1[k] <- FE2coef["x.1"]
    FE2_beta2[k] <- FE2coef["x.2"]
    FE2_beta3[k] <- FE2coef["x.3"]
    
    # Coefficient errors: beta_1
    FE2_rep_sqerr_x1[k]     <- (FE2_beta1[k] - b1True)^2
    FE2_rep_abserr_x1[k]    <- abs(FE2_beta1[k] - b1True)
    FE2_rep_pcterr_x1[k]    <- (FE2_beta1[k] - b1True) / b1True
    FE2_rep_abspcterr_x1[k] <- abs(FE2_rep_pcterr_x1[k])
    
    # Coefficient errors: beta_2
    FE2_rep_sqerr_x2[k]     <- (FE2_beta2[k] - b2True)^2
    FE2_rep_abserr_x2[k]    <- abs(FE2_beta2[k] - b2True)
    FE2_rep_pcterr_x2[k]    <- (FE2_beta2[k] - b2True) / b2True
    FE2_rep_abspcterr_x2[k] <- abs(FE2_rep_pcterr_x2[k])
    
    # Coefficient errors: beta_3
    FE2_rep_sqerr_x3[k]     <- (FE2_beta3[k] - b3True)^2
    FE2_rep_abserr_x3[k]    <- abs(FE2_beta3[k] - b3True)
    FE2_rep_pcterr_x3[k]    <- (FE2_beta3[k] - b3True) / b3True
    FE2_rep_abspcterr_x3[k] <- abs(FE2_rep_pcterr_x3[k])
    
    # Residual variance
    FE2_rep_resid_var[k] <- sigma(FE2)^2
    
    # Model-based SEs
    FE2_SE_beta0[k] <- FE2se["(Intercept)"]
    FE2_SE_beta1[k] <- FE2se["x.1"]
    FE2_SE_beta2[k] <- FE2se["x.2"]
    FE2_SE_beta3[k] <- FE2se["x.3"]
    
    # 95% CIs
    FE2_CIlo_beta0[k] <- FE2_beta0[k] - 1.96 * FE2_SE_beta0[k]
    FE2_CIhi_beta0[k] <- FE2_beta0[k] + 1.96 * FE2_SE_beta0[k]
    
    FE2_CIlo_beta1[k] <- FE2_beta1[k] - 1.96 * FE2_SE_beta1[k]
    FE2_CIhi_beta1[k] <- FE2_beta1[k] + 1.96 * FE2_SE_beta1[k]
    
    FE2_CIlo_beta2[k] <- FE2_beta2[k] - 1.96 * FE2_SE_beta2[k]
    FE2_CIhi_beta2[k] <- FE2_beta2[k] + 1.96 * FE2_SE_beta2[k]
    
    FE2_CIlo_beta3[k] <- FE2_beta3[k] - 1.96 * FE2_SE_beta3[k]
    FE2_CIhi_beta3[k] <- FE2_beta3[k] + 1.96 * FE2_SE_beta3[k]
  }
  
  ### [LOCATION-level] Random Effects ###
  if(model2==TRUE){
    # REl: location effects treated as random intercepts
    REl <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + (1 | country),
      data = X,
      family = gaussian()#,weights = w
    )
    
    # Extract model objects
    RElres  <- residuals(REl, type = "response")
    REly    <- model.response(model.frame(REl))
    RElcoef <- fixef(REl)$cond
    RElvcov <- vcov(REl)$cond
    RElse   <- sqrt(diag(RElvcov))
    
    # Residual performance
    REl_rep_resid_mse[k]  <- mean(RElres^2)
    REl_rep_resid_rmse[k] <- sqrt(mean(RElres^2))
    REl_rep_resid_mae[k]  <- mean(abs(RElres))
    REl_rep_resid_mpe[k]  <- mean(RElres / REly)
    REl_rep_resid_mape[k] <- mean(abs(RElres / REly))
    
    # Coefficient estimates
    REl_beta0[k] <- RElcoef["(Intercept)"]
    REl_beta1[k] <- RElcoef["x.1"]
    REl_beta2[k] <- RElcoef["x.2"]
    REl_beta3[k] <- RElcoef["x.3"]
    
    # beta_1
    REl_rep_sqerr_x1[k]     <- (REl_beta1[k] - b1True)^2
    REl_rep_abserr_x1[k]    <- abs(REl_beta1[k] - b1True)
    REl_rep_pcterr_x1[k]    <- (REl_beta1[k] - b1True) / b1True
    REl_rep_abspcterr_x1[k] <- abs(REl_rep_pcterr_x1[k])
    
    # beta_2
    REl_rep_sqerr_x2[k]     <- (REl_beta2[k] - b2True)^2
    REl_rep_abserr_x2[k]    <- abs(REl_beta2[k] - b2True)
    REl_rep_pcterr_x2[k]    <- (REl_beta2[k] - b2True) / b2True
    REl_rep_abspcterr_x2[k] <- abs(REl_rep_pcterr_x2[k])
    
    # beta_3
    REl_rep_sqerr_x3[k]     <- (REl_beta3[k] - b3True)^2
    REl_rep_abserr_x3[k]    <- abs(REl_beta3[k] - b3True)
    REl_rep_pcterr_x3[k]    <- (REl_beta3[k] - b3True) / b3True
    REl_rep_abspcterr_x3[k] <- abs(REl_rep_pcterr_x3[k])
    
    # Residual variance
    REl_rep_resid_var[k] <- sigma(REl)^2
    
    # Model-based SEs
    REl_SE_beta0[k] <- RElse["(Intercept)"]
    REl_SE_beta1[k] <- RElse["x.1"]
    REl_SE_beta2[k] <- RElse["x.2"]
    REl_SE_beta3[k] <- RElse["x.3"]
    
    # 95% CIs
    REl_CIlo_beta0[k] <- REl_beta0[k] - 1.96 * REl_SE_beta0[k]
    REl_CIhi_beta0[k] <- REl_beta0[k] + 1.96 * REl_SE_beta0[k]
    
    REl_CIlo_beta1[k] <- REl_beta1[k] - 1.96 * REl_SE_beta1[k]
    REl_CIhi_beta1[k] <- REl_beta1[k] + 1.96 * REl_SE_beta1[k]
    
    REl_CIlo_beta2[k] <- REl_beta2[k] - 1.96 * REl_SE_beta2[k]
    REl_CIhi_beta2[k] <- REl_beta2[k] + 1.96 * REl_SE_beta2[k]
    
    REl_CIlo_beta3[k] <- REl_beta3[k] - 1.96 * REl_SE_beta3[k]
    REl_CIhi_beta3[k] <- REl_beta3[k] + 1.96 * REl_SE_beta3[k]
  }
  
  ### [LOCATION-level] Fixed Effects ###
  if(model4==TRUE){
    # FE: location effects treated as fixed intercepts
    FEl <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + factor(country),
      data = X,
      family = gaussian()#,weights = w
    )

    # Extract model objects
    FElres  <- residuals(FEl, type = "response")
    FEly    <- model.response(model.frame(FEl))
    FElcoef <- fixef(FEl)$cond
    FElvcov <- vcov(FEl)$cond
    FElse   <- sqrt(diag(FElvcov))
    
    # Residual performance
    FEl_rep_resid_mse[k]  <- mean(FElres^2)
    FEl_rep_resid_rmse[k] <- sqrt(mean(FElres^2))
    FEl_rep_resid_mae[k]  <- mean(abs(FElres))
    FEl_rep_resid_mpe[k]  <- mean(FElres / FEly)
    FEl_rep_resid_mape[k] <- mean(abs(FElres / FEly))
    
    # Coefficient estimates
    FEl_beta0[k] <- FElcoef["(Intercept)"]
    FEl_beta1[k] <- FElcoef["x.1"]
    FEl_beta2[k] <- FElcoef["x.2"]
    FEl_beta3[k] <- FElcoef["x.3"]
    
    # beta_1
    FEl_rep_sqerr_x1[k]     <- (FEl_beta1[k] - b1True)^2
    FEl_rep_abserr_x1[k]    <- abs(FEl_beta1[k] - b1True)
    FEl_rep_pcterr_x1[k]    <- (FEl_beta1[k] - b1True) / b1True
    FEl_rep_abspcterr_x1[k] <- abs(FEl_rep_pcterr_x1[k])
    
    # beta_2
    FEl_rep_sqerr_x2[k]     <- (FEl_beta2[k] - b2True)^2
    FEl_rep_abserr_x2[k]    <- abs(FEl_beta2[k] - b2True)
    FEl_rep_pcterr_x2[k]    <- (FEl_beta2[k] - b2True) / b2True
    FEl_rep_abspcterr_x2[k] <- abs(FEl_rep_pcterr_x2[k])
    
    # beta_3
    FEl_rep_sqerr_x3[k]     <- (FEl_beta3[k] - b3True)^2
    FEl_rep_abserr_x3[k]    <- abs(FEl_beta3[k] - b3True)
    FEl_rep_pcterr_x3[k]    <- (FEl_beta3[k] - b3True) / b3True
    FEl_rep_abspcterr_x3[k] <- abs(FEl_rep_pcterr_x3[k])
    
    # Residual variance
    FEl_rep_resid_var[k] <- sigma(FEl)^2
    
    # Model-based SEs
    FEl_SE_beta0[k] <- FElse["(Intercept)"]
    FEl_SE_beta1[k] <- FElse["x.1"]
    FEl_SE_beta2[k] <- FElse["x.2"]
    FEl_SE_beta3[k] <- FElse["x.3"]
    
    # 95% CIs
    FEl_CIlo_beta0[k] <- FEl_beta0[k] - 1.96 * FEl_SE_beta0[k]
    FEl_CIhi_beta0[k] <- FEl_beta0[k] + 1.96 * FEl_SE_beta0[k]
    
    FEl_CIlo_beta1[k] <- FEl_beta1[k] - 1.96 * FEl_SE_beta1[k]
    FEl_CIhi_beta1[k] <- FEl_beta1[k] + 1.96 * FEl_SE_beta1[k]
    
    FEl_CIlo_beta2[k] <- FEl_beta2[k] - 1.96 * FEl_SE_beta2[k]
    FEl_CIhi_beta2[k] <- FEl_beta2[k] + 1.96 * FEl_SE_beta2[k]
    
    FEl_CIlo_beta3[k] <- FEl_beta3[k] - 1.96 * FEl_SE_beta3[k]
    FEl_CIhi_beta3[k] <- FEl_beta3[k] + 1.96 * FEl_SE_beta3[k]
  }
  
  ### [TIME-level] Fixed Effects ###
  if(model5==TRUE){
    # FEt: year effects treated as fixed intercepts
    FEt <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + factor(year),
      data = X,
      family = gaussian()#,weights = w
    )

    # Extract model objects
    FEtres  <- residuals(FEt, type = "response")
    FEty    <- model.response(model.frame(FEt))
    FEtcoef <- fixef(FEt)$cond
    FEtvcov <- vcov(FEt)$cond
    FEtse   <- sqrt(diag(FEtvcov))
    
    # Residual performance
    FEt_rep_resid_mse[k]  <- mean(FEtres^2)
    FEt_rep_resid_rmse[k] <- sqrt(mean(FEtres^2))
    FEt_rep_resid_mae[k]  <- mean(abs(FEtres))
    FEt_rep_resid_mpe[k]  <- mean(FEtres / FEty)
    FEt_rep_resid_mape[k] <- mean(abs(FEtres / FEty))
    
    # Coefficient estimates
    FEt_beta0[k] <- FEtcoef["(Intercept)"]
    FEt_beta1[k] <- FEtcoef["x.1"]
    FEt_beta2[k] <- FEtcoef["x.2"]
    FEt_beta3[k] <- FEtcoef["x.3"]
    
    # beta_1
    FEt_rep_sqerr_x1[k]     <- (FEt_beta1[k] - b1True)^2
    FEt_rep_abserr_x1[k]    <- abs(FEt_beta1[k] - b1True)
    FEt_rep_pcterr_x1[k]    <- (FEt_beta1[k] - b1True) / b1True
    FEt_rep_abspcterr_x1[k] <- abs(FEt_rep_pcterr_x1[k])
    
    # beta_2
    FEt_rep_sqerr_x2[k]     <- (FEt_beta2[k] - b2True)^2
    FEt_rep_abserr_x2[k]    <- abs(FEt_beta2[k] - b2True)
    FEt_rep_pcterr_x2[k]    <- (FEt_beta2[k] - b2True) / b2True
    FEt_rep_abspcterr_x2[k] <- abs(FEt_rep_pcterr_x2[k])
    
    # beta_3
    FEt_rep_sqerr_x3[k]     <- (FEt_beta3[k] - b3True)^2
    FEt_rep_abserr_x3[k]    <- abs(FEt_beta3[k] - b3True)
    FEt_rep_pcterr_x3[k]    <- (FEt_beta3[k] - b3True) / b3True
    FEt_rep_abspcterr_x3[k] <- abs(FEt_rep_pcterr_x3[k])
    
    # Residual variance
    FEt_rep_resid_var[k] <- sigma(FEt)^2
    
    # Model-based SEs
    FEt_SE_beta0[k] <- FEtse["(Intercept)"]
    FEt_SE_beta1[k] <- FEtse["x.1"]
    FEt_SE_beta2[k] <- FEtse["x.2"]
    FEt_SE_beta3[k] <- FEtse["x.3"]
    
    # 95% CIs
    FEt_CIlo_beta0[k] <- FEt_beta0[k] - 1.96 * FEt_SE_beta0[k]
    FEt_CIhi_beta0[k] <- FEt_beta0[k] + 1.96 * FEt_SE_beta0[k]
    
    FEt_CIlo_beta1[k] <- FEt_beta1[k] - 1.96 * FEt_SE_beta1[k]
    FEt_CIhi_beta1[k] <- FEt_beta1[k] + 1.96 * FEt_SE_beta1[k]
    
    FEt_CIlo_beta2[k] <- FEt_beta2[k] - 1.96 * FEt_SE_beta2[k]
    FEt_CIhi_beta2[k] <- FEt_beta2[k] + 1.96 * FEt_SE_beta2[k]
    
    FEt_CIlo_beta3[k] <- FEt_beta3[k] - 1.96 * FEt_SE_beta3[k]
    FEt_CIhi_beta3[k] <- FEt_beta3[k] + 1.96 * FEt_SE_beta3[k]
  }

  ### [Study-level] Fixed Effects (with Trend) ###
  if(model7==TRUE){
    # FEsT: paper effects treated as fixed intercepts
    FEsT <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + trend + factor(paper),
      data = X,
      family = gaussian()#,weights = w
    )
  }
  
  ### [Location-level] Fixed Effects (with Trend) ###
  if(model8==TRUE){
    # FElT: location effects treated as fixed intercepts
    FElT <- glmmTMB(
      y ~ x.1 + x.2 + x.3 + trend + factor(country),
      data = X,
      family = gaussian()#,weights = w
    )

    # Extract model objects
    FElTres  <- residuals(FElT, type = "response")
    FElTy    <- model.response(model.frame(FElT))
    FElTcoef <- fixef(FElT)$cond
    FElTvcov <- vcov(FElT)$cond
    FElTse   <- sqrt(diag(FElTvcov))
    
    # Residual performance
    FElT_rep_resid_mse[k]  <- mean(FElTres^2)
    FElT_rep_resid_rmse[k] <- sqrt(mean(FElTres^2))
    FElT_rep_resid_mae[k]  <- mean(abs(FElTres))
    FElT_rep_resid_mpe[k]  <- mean(FElTres / FElTy)
    FElT_rep_resid_mape[k] <- mean(abs(FElTres / FElTy))
    
    # Coefficient estimates
    FElT_beta0[k]     <- FElTcoef["(Intercept)"]
    FElT_beta1[k]     <- FElTcoef["x.1"]
    FElT_beta2[k]     <- FElTcoef["x.2"]
    FElT_beta3[k]     <- FElTcoef["x.3"]
    FElT_betaTrend[k] <- FElTcoef["trend"]
    
    # beta_1
    FElT_rep_sqerr_x1[k]     <- (FElT_beta1[k] - b1True)^2
    FElT_rep_abserr_x1[k]    <- abs(FElT_beta1[k] - b1True)
    FElT_rep_pcterr_x1[k]    <- (FElT_beta1[k] - b1True) / b1True
    FElT_rep_abspcterr_x1[k] <- abs(FElT_rep_pcterr_x1[k])
    
    # beta_2
    FElT_rep_sqerr_x2[k]     <- (FElT_beta2[k] - b2True)^2
    FElT_rep_abserr_x2[k]    <- abs(FElT_beta2[k] - b2True)
    FElT_rep_pcterr_x2[k]    <- (FElT_beta2[k] - b2True) / b2True
    FElT_rep_abspcterr_x2[k] <- abs(FElT_rep_pcterr_x2[k])
    
    # beta_3
    FElT_rep_sqerr_x3[k]     <- (FElT_beta3[k] - b3True)^2
    FElT_rep_abserr_x3[k]    <- abs(FElT_beta3[k] - b3True)
    FElT_rep_pcterr_x3[k]    <- (FElT_beta3[k] - b3True) / b3True
    FElT_rep_abspcterr_x3[k] <- abs(FElT_rep_pcterr_x3[k])
    
    # Trend
    trendTrue <- 2 * timeHet
    
    FElT_rep_sqerr_trend[k] <-
      (FElT_betaTrend[k] - trendTrue)^2
    
    FElT_rep_abserr_trend[k] <-
      abs(FElT_betaTrend[k] - trendTrue)
    
    FElT_rep_pcterr_trend[k] <-
      (FElT_betaTrend[k] - trendTrue) / trendTrue
    
    FElT_rep_abspcterr_trend[k] <-
      abs(FElT_rep_pcterr_trend[k])
    
    # Residual variance
    FElT_rep_resid_var[k] <- sigma(FElT)^2
    
    # Model-based SEs
    FElT_SE_beta0[k]     <- FElTse["(Intercept)"]
    FElT_SE_beta1[k]     <- FElTse["x.1"]
    FElT_SE_beta2[k]     <- FElTse["x.2"]
    FElT_SE_beta3[k]     <- FElTse["x.3"]
    FElT_SE_betaTrend[k] <- FElTse["trend"]
    
    # 95% CIs
    FElT_CIlo_beta0[k] <- FElT_beta0[k] - 1.96 * FElT_SE_beta0[k]
    FElT_CIhi_beta0[k] <- FElT_beta0[k] + 1.96 * FElT_SE_beta0[k]
    
    FElT_CIlo_beta1[k] <- FElT_beta1[k] - 1.96 * FElT_SE_beta1[k]
    FElT_CIhi_beta1[k] <- FElT_beta1[k] + 1.96 * FElT_SE_beta1[k]
    
    FElT_CIlo_beta2[k] <- FElT_beta2[k] - 1.96 * FElT_SE_beta2[k]
    FElT_CIhi_beta2[k] <- FElT_beta2[k] + 1.96 * FElT_SE_beta2[k]
    
    FElT_CIlo_beta3[k] <- FElT_beta3[k] - 1.96 * FElT_SE_beta3[k]
    FElT_CIhi_beta3[k] <- FElT_beta3[k] + 1.96 * FElT_SE_beta3[k]
    
    FElT_CIlo_betaTrend[k] <-
      FElT_betaTrend[k] - 1.96 * FElT_SE_betaTrend[k]
    
    FElT_CIhi_betaTrend[k] <-
      FElT_betaTrend[k] + 1.96 * FElT_SE_betaTrend[k]
  }
}
# Audio to signal Monte Carlo completed
beep(sound = 5)

end<-Sys.time()
print(end-start)



#######################################
### MONTE CARLO PERFORMANCE SUMMARY ### (no longer just k)
#######################################

### REs ###
# beta_1
RE_MC_bias_beta1 <- mean(RE_beta1, na.rm = TRUE) - b1True
RE_MC_empSD_beta1 <- sd(RE_beta1, na.rm = TRUE)
RE_MC_meanSE_beta1 <- mean(RE_SE_beta1, na.rm = TRUE)
RE_MC_mse_beta1 <- mean(RE_rep_sqerr_x1, na.rm = TRUE)
RE_MC_rmse_beta1 <- sqrt(RE_MC_mse_beta1)

RE_MC_CIcoverage_beta1 <- mean(
  RE_CIlo_beta1 <= b1True & RE_CIhi_beta1 >= b1True,
  na.rm = TRUE)

RE_MC_CIwidth_beta1 <- mean(
  RE_CIhi_beta1 - RE_CIlo_beta1,
  na.rm = TRUE)

# beta_2
RE_MC_bias_beta2 <- mean(RE_beta2, na.rm = TRUE) - b2True
RE_MC_empSD_beta2 <- sd(RE_beta2, na.rm = TRUE)
RE_MC_meanSE_beta2 <- mean(RE_SE_beta2, na.rm = TRUE)
RE_MC_mse_beta2 <- mean(RE_rep_sqerr_x2, na.rm = TRUE)
RE_MC_rmse_beta2 <- sqrt(RE_MC_mse_beta2)

RE_MC_CIcoverage_beta2 <- mean(
  RE_CIlo_beta2 <= b2True & RE_CIhi_beta2 >= b2True,
  na.rm = TRUE)

RE_MC_CIwidth_beta2 <- mean(
  RE_CIhi_beta2 - RE_CIlo_beta2,
  na.rm = TRUE)

# beta_3 
RE_MC_bias_beta3 <- mean(RE_beta3, na.rm = TRUE) - b3True
RE_MC_empSD_beta3 <- sd(RE_beta3, na.rm = TRUE)
RE_MC_meanSE_beta3 <- mean(RE_SE_beta3, na.rm = TRUE)
RE_MC_mse_beta3 <- mean(RE_rep_sqerr_x3, na.rm = TRUE)
RE_MC_rmse_beta3 <- sqrt(RE_MC_mse_beta3)

RE_MC_CIcoverage_beta3 <- mean(
  RE_CIlo_beta3 <= b3True & RE_CIhi_beta3 >= b3True,
  na.rm = TRUE)

RE_MC_CIwidth_beta3 <- mean(
  RE_CIhi_beta3 - RE_CIlo_beta3,
  na.rm = TRUE)


### FEs ###
# beta_1

FE_MC_bias_beta1 <- mean(FE_beta1, na.rm = TRUE) - b1True
FE_MC_empSD_beta1 <- sd(FE_beta1, na.rm = TRUE)
FE_MC_meanSE_beta1 <- mean(FE_SE_beta1, na.rm = TRUE)
FE_MC_mse_beta1 <- mean(FE_rep_sqerr_x1, na.rm = TRUE)
FE_MC_rmse_beta1 <- sqrt(FE_MC_mse_beta1)

FE_MC_CIcoverage_beta1 <- mean(
  FE_CIlo_beta1 <= b1True & FE_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

FE_MC_CIwidth_beta1 <- mean(
  FE_CIhi_beta1 - FE_CIlo_beta1,
  na.rm = TRUE
)

# beta_2 ----------------
FE_MC_bias_beta2 <- mean(FE_beta2, na.rm = TRUE) - b2True
FE_MC_empSD_beta2 <- sd(FE_beta2, na.rm = TRUE)
FE_MC_meanSE_beta2 <- mean(FE_SE_beta2, na.rm = TRUE)
FE_MC_mse_beta2 <- mean(FE_rep_sqerr_x2, na.rm = TRUE)
FE_MC_rmse_beta2 <- sqrt(FE_MC_mse_beta2)

FE_MC_CIcoverage_beta2 <- mean(
  FE_CIlo_beta2 <= b2True & FE_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

FE_MC_CIwidth_beta2 <- mean(
  FE_CIhi_beta2 - FE_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
FE_MC_bias_beta3 <- mean(FE_beta3, na.rm = TRUE) - b3True
FE_MC_empSD_beta3 <- sd(FE_beta3, na.rm = TRUE)
FE_MC_meanSE_beta3 <- mean(FE_SE_beta3, na.rm = TRUE)
FE_MC_mse_beta3 <- mean(FE_rep_sqerr_x3, na.rm = TRUE)
FE_MC_rmse_beta3 <- sqrt(FE_MC_mse_beta3)

FE_MC_CIcoverage_beta3 <- mean(
  FE_CIlo_beta3 <= b3True & FE_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

FE_MC_CIwidth_beta3 <- mean(
  FE_CIhi_beta3 - FE_CIlo_beta3,
  na.rm = TRUE
)


### FE2 ###
# beta_1
FE2_MC_bias_beta1 <- mean(FE2_beta1, na.rm = TRUE) - b1True
FE2_MC_empSD_beta1 <- sd(FE2_beta1, na.rm = TRUE)
FE2_MC_meanSE_beta1 <- mean(FE2_SE_beta1, na.rm = TRUE)
FE2_MC_mse_beta1 <- mean(FE2_rep_sqerr_x1, na.rm = TRUE)
FE2_MC_rmse_beta1 <- sqrt(FE2_MC_mse_beta1)

FE2_MC_CIcoverage_beta1 <- mean(
  FE2_CIlo_beta1 <= b1True & FE2_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

FE2_MC_CIwidth_beta1 <- mean(
  FE2_CIhi_beta1 - FE2_CIlo_beta1,
  na.rm = TRUE
)

# beta_2
FE2_MC_bias_beta2 <- mean(FE2_beta2, na.rm = TRUE) - b2True
FE2_MC_empSD_beta2 <- sd(FE2_beta2, na.rm = TRUE)
FE2_MC_meanSE_beta2 <- mean(FE2_SE_beta2, na.rm = TRUE)
FE2_MC_mse_beta2 <- mean(FE2_rep_sqerr_x2, na.rm = TRUE)
FE2_MC_rmse_beta2 <- sqrt(FE2_MC_mse_beta2)

FE2_MC_CIcoverage_beta2 <- mean(
  FE2_CIlo_beta2 <= b2True & FE2_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

FE2_MC_CIwidth_beta2 <- mean(
  FE2_CIhi_beta2 - FE2_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
FE2_MC_bias_beta3 <- mean(FE2_beta3, na.rm = TRUE) - b3True
FE2_MC_empSD_beta3 <- sd(FE2_beta3, na.rm = TRUE)
FE2_MC_meanSE_beta3 <- mean(FE2_SE_beta3, na.rm = TRUE)
FE2_MC_mse_beta3 <- mean(FE2_rep_sqerr_x3, na.rm = TRUE)
FE2_MC_rmse_beta3 <- sqrt(FE2_MC_mse_beta3)

FE2_MC_CIcoverage_beta3 <- mean(
  FE2_CIlo_beta3 <= b3True & FE2_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

FE2_MC_CIwidth_beta3 <- mean(
  FE2_CIhi_beta3 - FE2_CIlo_beta3,
  na.rm = TRUE
)

### REl ###
# beta_1
REl_MC_bias_beta1 <- mean(REl_beta1, na.rm = TRUE) - b1True
REl_MC_empSD_beta1 <- sd(REl_beta1, na.rm = TRUE)
REl_MC_meanSE_beta1 <- mean(REl_SE_beta1, na.rm = TRUE)
REl_MC_mse_beta1 <- mean(REl_rep_sqerr_x1, na.rm = TRUE)
REl_MC_rmse_beta1 <- sqrt(REl_MC_mse_beta1)

REl_MC_CIcoverage_beta1 <- mean(
  REl_CIlo_beta1 <= b1True & REl_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

REl_MC_CIwidth_beta1 <- mean(
  REl_CIhi_beta1 - REl_CIlo_beta1,
  na.rm = TRUE
)

# beta_2
REl_MC_bias_beta2 <- mean(REl_beta2, na.rm = TRUE) - b2True
REl_MC_empSD_beta2 <- sd(REl_beta2, na.rm = TRUE)
REl_MC_meanSE_beta2 <- mean(REl_SE_beta2, na.rm = TRUE)
REl_MC_mse_beta2 <- mean(REl_rep_sqerr_x2, na.rm = TRUE)
REl_MC_rmse_beta2 <- sqrt(REl_MC_mse_beta2)

REl_MC_CIcoverage_beta2 <- mean(
  REl_CIlo_beta2 <= b2True & REl_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

REl_MC_CIwidth_beta2 <- mean(
  REl_CIhi_beta2 - REl_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
REl_MC_bias_beta3 <- mean(REl_beta3, na.rm = TRUE) - b3True
REl_MC_empSD_beta3 <- sd(REl_beta3, na.rm = TRUE)
REl_MC_meanSE_beta3 <- mean(REl_SE_beta3, na.rm = TRUE)
REl_MC_mse_beta3 <- mean(REl_rep_sqerr_x3, na.rm = TRUE)
REl_MC_rmse_beta3 <- sqrt(REl_MC_mse_beta3)

REl_MC_CIcoverage_beta3 <- mean(
  REl_CIlo_beta3 <= b3True & REl_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

REl_MC_CIwidth_beta3 <- mean(
  REl_CIhi_beta3 - REl_CIlo_beta3,
  na.rm = TRUE
)

### FEl ###
# beta_1
FEl_MC_bias_beta1 <- mean(FEl_beta1, na.rm = TRUE) - b1True
FEl_MC_empSD_beta1 <- sd(FEl_beta1, na.rm = TRUE)
FEl_MC_meanSE_beta1 <- mean(FEl_SE_beta1, na.rm = TRUE)
FEl_MC_mse_beta1 <- mean(FEl_rep_sqerr_x1, na.rm = TRUE)
FEl_MC_rmse_beta1 <- sqrt(FEl_MC_mse_beta1)

FEl_MC_CIcoverage_beta1 <- mean(
  FEl_CIlo_beta1 <= b1True & FEl_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

FEl_MC_CIwidth_beta1 <- mean(
  FEl_CIhi_beta1 - FEl_CIlo_beta1,
  na.rm = TRUE
)

# beta_2
FEl_MC_bias_beta2 <- mean(FEl_beta2, na.rm = TRUE) - b2True
FEl_MC_empSD_beta2 <- sd(FEl_beta2, na.rm = TRUE)
FEl_MC_meanSE_beta2 <- mean(FEl_SE_beta2, na.rm = TRUE)
FEl_MC_mse_beta2 <- mean(FEl_rep_sqerr_x2, na.rm = TRUE)
FEl_MC_rmse_beta2 <- sqrt(FEl_MC_mse_beta2)

FEl_MC_CIcoverage_beta2 <- mean(
  FEl_CIlo_beta2 <= b2True & FEl_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

FEl_MC_CIwidth_beta2 <- mean(
  FEl_CIhi_beta2 - FEl_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
FEl_MC_bias_beta3 <- mean(FEl_beta3, na.rm = TRUE) - b3True
FEl_MC_empSD_beta3 <- sd(FEl_beta3, na.rm = TRUE)
FEl_MC_meanSE_beta3 <- mean(FEl_SE_beta3, na.rm = TRUE)
FEl_MC_mse_beta3 <- mean(FEl_rep_sqerr_x3, na.rm = TRUE)
FEl_MC_rmse_beta3 <- sqrt(FEl_MC_mse_beta3)

FEl_MC_CIcoverage_beta3 <- mean(
  FEl_CIlo_beta3 <= b3True & FEl_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

FEl_MC_CIwidth_beta3 <- mean(
  FEl_CIhi_beta3 - FEl_CIlo_beta3,
  na.rm = TRUE
)

### FEt ###
# beta_1
FEt_MC_bias_beta1 <- mean(FEt_beta1, na.rm = TRUE) - b1True
FEt_MC_empSD_beta1 <- sd(FEt_beta1, na.rm = TRUE)
FEt_MC_meanSE_beta1 <- mean(FEt_SE_beta1, na.rm = TRUE)
FEt_MC_mse_beta1 <- mean(FEt_rep_sqerr_x1, na.rm = TRUE)
FEt_MC_rmse_beta1 <- sqrt(FEt_MC_mse_beta1)

FEt_MC_CIcoverage_beta1 <- mean(
  FEt_CIlo_beta1 <= b1True & FEt_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

FEt_MC_CIwidth_beta1 <- mean(
  FEt_CIhi_beta1 - FEt_CIlo_beta1,
  na.rm = TRUE
)

# beta_2
FEt_MC_bias_beta2 <- mean(FEt_beta2, na.rm = TRUE) - b2True
FEt_MC_empSD_beta2 <- sd(FEt_beta2, na.rm = TRUE)
FEt_MC_meanSE_beta2 <- mean(FEt_SE_beta2, na.rm = TRUE)
FEt_MC_mse_beta2 <- mean(FEt_rep_sqerr_x2, na.rm = TRUE)
FEt_MC_rmse_beta2 <- sqrt(FEt_MC_mse_beta2)

FEt_MC_CIcoverage_beta2 <- mean(
  FEt_CIlo_beta2 <= b2True & FEt_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

FEt_MC_CIwidth_beta2 <- mean(
  FEt_CIhi_beta2 - FEt_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
FEt_MC_bias_beta3 <- mean(FEt_beta3, na.rm = TRUE) - b3True
FEt_MC_empSD_beta3 <- sd(FEt_beta3, na.rm = TRUE)
FEt_MC_meanSE_beta3 <- mean(FEt_SE_beta3, na.rm = TRUE)
FEt_MC_mse_beta3 <- mean(FEt_rep_sqerr_x3, na.rm = TRUE)
FEt_MC_rmse_beta3 <- sqrt(FEt_MC_mse_beta3)

FEt_MC_CIcoverage_beta3 <- mean(
  FEt_CIlo_beta3 <= b3True & FEt_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

FEt_MC_CIwidth_beta3 <- mean(
  FEt_CIhi_beta3 - FEt_CIlo_beta3,
  na.rm = TRUE
)

### FElT ###
# beta_1
FElT_MC_bias_beta1 <- mean(FElT_beta1, na.rm = TRUE) - b1True
FElT_MC_empSD_beta1 <- sd(FElT_beta1, na.rm = TRUE)
FElT_MC_meanSE_beta1 <- mean(FElT_SE_beta1, na.rm = TRUE)
FElT_MC_mse_beta1 <- mean(FElT_rep_sqerr_x1, na.rm = TRUE)
FElT_MC_rmse_beta1 <- sqrt(FElT_MC_mse_beta1)

FElT_MC_CIcoverage_beta1 <- mean(
  FElT_CIlo_beta1 <= b1True & FElT_CIhi_beta1 >= b1True,
  na.rm = TRUE
)

FElT_MC_CIwidth_beta1 <- mean(
  FElT_CIhi_beta1 - FElT_CIlo_beta1,
  na.rm = TRUE
)

# beta_2
FElT_MC_bias_beta2 <- mean(FElT_beta2, na.rm = TRUE) - b2True
FElT_MC_empSD_beta2 <- sd(FElT_beta2, na.rm = TRUE)
FElT_MC_meanSE_beta2 <- mean(FElT_SE_beta2, na.rm = TRUE)
FElT_MC_mse_beta2 <- mean(FElT_rep_sqerr_x2, na.rm = TRUE)
FElT_MC_rmse_beta2 <- sqrt(FElT_MC_mse_beta2)

FElT_MC_CIcoverage_beta2 <- mean(
  FElT_CIlo_beta2 <= b2True & FElT_CIhi_beta2 >= b2True,
  na.rm = TRUE
)

FElT_MC_CIwidth_beta2 <- mean(
  FElT_CIhi_beta2 - FElT_CIlo_beta2,
  na.rm = TRUE
)

# beta_3
FElT_MC_bias_beta3 <- mean(FElT_beta3, na.rm = TRUE) - b3True
FElT_MC_empSD_beta3 <- sd(FElT_beta3, na.rm = TRUE)
FElT_MC_meanSE_beta3 <- mean(FElT_SE_beta3, na.rm = TRUE)
FElT_MC_mse_beta3 <- mean(FElT_rep_sqerr_x3, na.rm = TRUE)
FElT_MC_rmse_beta3 <- sqrt(FElT_MC_mse_beta3)

FElT_MC_CIcoverage_beta3 <- mean(
  FElT_CIlo_beta3 <= b3True & FElT_CIhi_beta3 >= b3True,
  na.rm = TRUE
)

FElT_MC_CIwidth_beta3 <- mean(
  FElT_CIhi_beta3 - FElT_CIlo_beta3,
  na.rm = TRUE
)

#%%###########################################################################
###                             RESULTS - MONTE CARLO                      ###
##############################################################################

directory <- "C:/Users/gendronj1/OneDrive - Xavier University/SYNC/School/VT/RESEARCH/Dissertation/CH2 - 3YP/Code/Simulation/Simulation Results - R Files/x3/"

### To generate the last dataset (X and y)
# write.csv(
#   X,
#   file = file.path(
#     directory,
#     paste0("Generated_X_(Case", case, ").csv")
#   ),
#   row.names = FALSE
# )
# 
# write.csv(
#   data.frame(y = y),
#   file = file.path(
#     directory,
#     paste0("Generated_y_(Case", case, ").csv")
#   ),
#   row.names = FALSE
# )

###############
### TABLE 1 ### (Overall, all performance criteria for all models)
###############
models <- c(
  "RE",
  "FE",
  "FE2",
  "REl",
  "FEl",
  "FEt",
  "FElT"
)

true_betas <- c(
  b1True,
  b2True,
  b3True
)

MC_results <- purrr::map_dfr(
  models,
  function(model_name) {
    
    purrr::map_dfr(
      1:3,
      function(beta_num) {
        
        # Pull corresponding Monte Carlo vectors
        beta_est <- get(
          paste0(model_name, "_beta", beta_num)
        )
        
        beta_se <- get(
          paste0(model_name, "_SE_beta", beta_num)
        )
        
        sqerr <- get(
          paste0(model_name, "_rep_sqerr_x", beta_num)
        )
        
        CIlo <- get(
          paste0(model_name, "_CIlo_beta", beta_num)
        )
        
        CIhi <- get(
          paste0(model_name, "_CIhi_beta", beta_num)
        )
        
        beta_true <- true_betas[beta_num]
        
        # Create one row per model × beta
        data.frame(
          Case = case,
          #Countries = numCountries,
          #Years = numYears,
          #Study_n = n,
          #TimeHet = timeHet,
          #LocationSpread = spreadMuY,
          
          Model = model_name,
          Beta = paste0("Beta ", beta_num),
          
          TrueValue = beta_true,
          
          Bias =
            mean(beta_est, na.rm = TRUE) - beta_true,
          
          EmpSD =
            sd(beta_est, na.rm = TRUE),
          
          MeanSE =
            mean(beta_se, na.rm = TRUE),
          
          MSE =
            mean(sqerr, na.rm = TRUE),
          
          RMSE =
            sqrt(mean(sqerr, na.rm = TRUE)),
          
          CI_Coverage =
            mean(
              CIlo <= beta_true &
                CIhi >= beta_true,
              na.rm = TRUE
            ),
          
          CI_Width =
            mean(
              CIhi - CIlo,
              na.rm = TRUE
            )
        )
      }
    )
  }
)

# Set display order
MC_results <- MC_results %>%
  mutate(
    Model = factor(
      Model,
      levels = c(
        "RE",
        "FE",
        "FE2",
        "REl",
        "FEl",
        "FEt",
        "FElT"
      )
    ),
    
    Beta = factor(
      Beta,
      levels = c(
        "Beta 1",
        "Beta 2",
        "Beta 3"
      )
    )
  )

# For testing
#MC_results

file_name <- paste0(
  "Table_MC_results_",
  "(Case", case,
  #"_countries", numCountries,
  #"_n", n,
  ").csv"
)

write.csv(
  MC_results,
  file = file.path(directory, file_name),
  row.names = FALSE
)

###############
### FIGURES ### (per performance criteria)
###############

### SET UP ###

library(ggplot2)
library(dplyr)

### Figure 1 - Bias ###
plot_bias <- ggplot(
  MC_results,
  aes(
    x = Bias,
    y = Model
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Beta,
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    title = "Monte Carlo Bias by Model",
    subtitle = paste0(
      "Case ", case#,
      #"; Countries = ", numCountries,
      #"; n = ", n
    ),
    x = "Bias",
    y = "Model"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# For testing
#plot_bias

ggsave(
  filename = file.path(
    directory,
    paste0(
      "Figure_Bias_(Case",
      case,
      #"_n",
      #n,
      ").pdf"
    )
  ),
  plot = plot_bias,
  width = 10,
  height = 5
)

# ggsave(
#   filename = file.path(
#     directory,
#     paste0(
#       "Figure_Bias_case",
#       case,
#       #"_n",
#       #n,
#       ".png"
#     )
#   ),
#   plot = plot_bias,
#   width = 10,
#   height = 5,
#   dpi = 300
# )



### Figure 2 - RMSE ###
plot_rmse <- ggplot(
  MC_results,
  aes(
    x = RMSE,
    y = Model
  )
) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Beta,
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    title = "Monte Carlo RMSE by Model",
    subtitle = paste0(
      "Case ", case#,
      #"; Countries = ", numCountries,
      #"; n = ", n
    ),
    x = "Root Mean Squared Error",
    y = "Model"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# For testing
#plot_rmse


# Save as PDF
ggsave(
  filename = file.path(
    directory,
    paste0(
      "Figure_RMSE_(Case",
      case,
      #"_n",
      #n,
      ").pdf"
    )
  ),
  plot = plot_rmse,
  width = 10,
  height = 5
)

# Save as PNG
# ggsave(
#   filename = file.path(
#     directory,
#     paste0(
#       "Figure_RMSE_case",
#       case,
#       #"_n",
#       #n,
#       ".png"
#     )
#   ),
#   plot = plot_rmse,
#   width = 10,
#   height = 5,
#   dpi = 300
# )


### Figure 3 - CI coverage ###
plot_coverage <- ggplot(
  MC_results,
  aes(
    x = CI_Coverage,
    y = Model
  )
) +
  geom_vline(
    xintercept = 0.95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Beta,
    nrow = 1
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.1)
  ) +
  labs(
    title = "95% Confidence Interval Coverage by Model",
    subtitle = paste0(
      "Case ", case#,
      #"; Countries = ", numCountries,
      #"; n = ", n
    ),
    x = "Empirical Coverage Probability",
    y = "Model"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# For testing
#plot_coverage


# Save as PDF
ggsave(
  filename = file.path(
    directory,
    paste0(
      "Figure_CI_Coverage_(Case",
      case,
      #"_n",
      #n,
      ").pdf"
    )
  ),
  plot = plot_coverage,
  width = 10,
  height = 5
)

# Save as PNG
# ggsave(
#   filename = file.path(
#     directory,
#     paste0(
#       "Figure_CI_Coverage_case",
#       case,
#       #"_n",
#       #n,
#       ".png"
#     )
#   ),
#   plot = plot_coverage,
#   width = 10,
#   height = 5,
#   dpi = 300
# )


### Figure 4 - Empirical SD vs Mean SE ###
plot_SEcalibration <- ggplot(
  MC_results,
  aes(
    x = EmpSD,
    y = MeanSE,
    shape = Model
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  geom_point(
    size = 3
  ) +
  scale_shape_manual(
    values = c(
      "RE"   = 16,
      "FE"   = 17,
      "FE2"  = 15,
      "REl"  = 18,
      "FEl"  = 3,
      "FEt"  = 4,
      "FElT" = 8
    )
  ) +
  facet_wrap(
    ~ Beta,
    nrow = 1,
    scales = "free"
  ) +
  labs(
    title = "Model-Based Standard Errors versus Empirical Monte Carlo SD",
    subtitle = paste0(
      "Case ", case#,
      #"; Countries = ", numCountries,
      #"; n = ", n
    ),
    x = "Empirical Standard Deviation",
    y = "Mean Model-Based Standard Error",
    shape = "Model"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# For testing
#plot_SEcalibration

# Save as PDF
ggsave(
  filename = file.path(
    directory,
    paste0(
      "Figure_SE_Calibration_(Case",
      case,
      #"_n",
      #n,
      ").pdf"
    )
  ),
  plot = plot_SEcalibration,
  width = 10,
  height = 5
)

# Save as PNG
# ggsave(
#   filename = file.path(
#     directory,
#     paste0(
#       "Figure_SE_Calibration_case",
#       case,
#       #"_n",
#       #n,
#       ".png"
#     )
#   ),
#   plot = plot_SEcalibration,
#   width = 10,
#   height = 5,
#   dpi = 300
# )


### Figure 5 - CI width ###
plot_CIwidth <- ggplot(
  MC_results,
  aes(
    x = CI_Width,
    y = Model
  )
) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ Beta,
    nrow = 1,
    scales = "free_x"
  ) +
  labs(
    title = "Average 95% Confidence Interval Width by Model",
    subtitle = paste0(
      "Case ", case#,
      #"; Countries = ", numCountries,
      #"; n = ", n
    ),
    x = "Average 95% CI Width",
    y = "Model"
  ) +
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

# For testing
#plot_CIwidth

# Save as PDF
ggsave(
  filename = file.path(
    directory,
    paste0(
      "Figure_CI_Width_(Case",
      case,
      #"_n",
      #n,
      ").pdf"
    )
  ),
  plot = plot_CIwidth,
  width = 10,
  height = 5
)

# Save as PNG
# ggsave(
#   filename = file.path(
#     directory,
#     paste0(
#       "Figure_CI_Width_case",
#       case,
#       #"_n",
#       #n,
#       ".png"
#     )
#   ),
#   plot = plot_CIwidth,
#   width = 10,
#   height = 5,
#   dpi = 300
# )



#%%###########################################################################
###                           EXPORT - Save ALL Data                       ###
##############################################################################
scenario_id <- paste0(
  "(Case", case#,
  #"_countries", numCountries,
  #"_n", n
)

# Identify all Monte Carlo vectors generated for the retained models
MC_vectors <- ls(
  pattern = "^(RE|FE|FE2|REl|FEl|FEt|FElT)_"
)

# Add simulation settings, true values, and summary results
objects_to_save <- unique(
  c(
    MC_vectors,
    
    # Simulation settings
    "case",
    "N",
    "n",
    "numCountries",
    "numYears",
    "timeHet",
    "spreadMuY",
    
    # True parameter values
    "b1True",
    "b2True",
    "b3True",
    "sig2True",
    "R2True",
    
    # Convenient finished summary table
    "MC_results"
  )
)

# Save everything in ONE file
save(
  list = objects_to_save,
  file = file.path(
    directory,
    paste0(
      "DATA-MC_ALL_DATA_",
      scenario_id,
      ").RData"
    )
  )
)



#%%########################################################
###                         POWER                       ###
###########################################################

#################
### STRUCTURE ###
#################

# 1. POWER FUNCTIONS
#   1A. Calculate empirical power
#   1B. Plot empirical power
#   1C. Save power figure/data
# 
# 2. POWER SETUP
#   2A. Define models
#   2B. Define discrepancy range
#   2C. Define alpha levels
#   2D. Define output directory
# 
# 3. RUN Adjusted alpha
#
# 4. RUN Traditional alpha (fixed @ 0.05)

##########################
### SET UP - Functions ###
##########################

### 1A. Calculate Empirical Power ###

empirical_power_curve <- function(
    beta_est,
    se_est,
    true_beta,
    discrepancy_grid,
    alpha_level
) {
  
  power_vals <- sapply(
    discrepancy_grid,
    function(d) {
      
      beta_null <- true_beta + d
      
      z_stat <- (beta_est - beta_null) / se_est
      
      p_val <- 2 * pnorm(-abs(z_stat))
      
      mean(
        p_val < alpha_level,
        na.rm = TRUE
      )
    }
  )
  
  data.frame(
    discrepancy = discrepancy_grid,
    power = power_vals
  )
}


### 1B. Plot Empirical Power ###

plot_empirical_power_panel <- function(
    beta_est,
    se_est,
    true_beta,
    discrepancy_grid,
    alpha_level,
    main_title,
    y_title
) {
  
  curve_df <- empirical_power_curve(
    beta_est,
    se_est,
    true_beta,
    discrepancy_grid,
    alpha_level
  )
  
  plot(
    curve_df$discrepancy,
    curve_df$power,
    type = "l",
    lwd = 2,
    col = "steelblue4",
    xlab = "Discrepancy from True Coefficient",
    ylab = paste("Rejection Probability for", y_title),
    main = main_title,
    ylim = c(0, 1),
    xlim = range(discrepancy_grid)
  )
  
  # Horizontal line = nominal significance level
  abline(
    h = alpha_level,
    lty = 3,
    col = "gray40"
  )
  
  # Vertical line = true coefficient / zero discrepancy
  abline(
    v = 0,
    lty = 3,
    col = "gray40"
  )
  
  # Display alpha below model title
  mtext(
    bquote(alpha == .(format(alpha_level, scientific = TRUE))),
    #paste0(
    #  "alpha = ",
    #  format(alpha_level, scientific = TRUE)
    #),
    side = 3,
    line = -0.15,
    cex = 0.65
  )
}



### 1C. Save Power Figure ###

save_empirical_power_figure <- function(
    model_list,
    beta_num,
    true_beta,
    alpha_level,
    scenario_label,
    out_dir,
    file_stub,
    discrepancy_grid
) {
  
  # Publication-quality PDF
  pdf(
    file = file.path(
      out_dir,
      paste0(file_stub, ".pdf")
    ),
    width = 14,
    height = 8,
    pointsize = 11,
    useDingbats = FALSE
  )
  
  par(
    mfrow = c(2, 4),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0),
    cex.axis = 0.9,
    cex.lab = 0.95,
    cex.main = 0.95
  )
  
  for (m in model_list) {
    
    plot_empirical_power_panel(
      beta_est = m[[paste0("beta", beta_num)]],
      se_est = m[[paste0("se", beta_num)]],
      true_beta = true_beta,
      discrepancy_grid = discrepancy_grid,
      alpha_level = alpha_level,
      main_title = m$title,
      y_title = paste0("Beta ", beta_num)
    )
  }
  
  # Blank eighth panel
  plot.new()
  
  mtext(
    paste0(
      "Empirical Power Curves: Beta ",
      beta_num,
      " (",
      scenario_label,
      ")"
    ),
    outer = TRUE,
    cex = 1.3,
    font = 2
  )
  
  dev.off()
}



##########################
### SET UP - Remainder ###
##########################

### 2A. Define Models ###
model_results <- list(
  
  RE = list(
    title = "Study Random Effects",
    beta1 = RE_beta1,
    se1 = RE_SE_beta1,
    beta2 = RE_beta2,
    se2 = RE_SE_beta2,
    beta3 = RE_beta3,
    se3 = RE_SE_beta3
  ),
  
  FE = list(
    title = "Study Fixed Effects",
    beta1 = FE_beta1,
    se1 = FE_SE_beta1,
    beta2 = FE_beta2,
    se2 = FE_SE_beta2,
    beta3 = FE_beta3,
    se3 = FE_SE_beta3
  ),
  
  FE2 = list(
    title = "Two-Way Fixed Effects",
    beta1 = FE2_beta1,
    se1 = FE2_SE_beta1,
    beta2 = FE2_beta2,
    se2 = FE2_SE_beta2,
    beta3 = FE2_beta3,
    se3 = FE2_SE_beta3
  ),
  
  REl = list(
    title = "Location Random Effects",
    beta1 = REl_beta1,
    se1 = REl_SE_beta1,
    beta2 = REl_beta2,
    se2 = REl_SE_beta2,
    beta3 = REl_beta3,
    se3 = REl_SE_beta3
  ),
  
  FEl = list(
    title = "Location Fixed Effects",
    beta1 = FEl_beta1,
    se1 = FEl_SE_beta1,
    beta2 = FEl_beta2,
    se2 = FEl_SE_beta2,
    beta3 = FEl_beta3,
    se3 = FEl_SE_beta3
  ),
  
  FEt = list(
    title = "Time Fixed Effects",
    beta1 = FEt_beta1,
    se1 = FEt_SE_beta1,
    beta2 = FEt_beta2,
    se2 = FEt_SE_beta2,
    beta3 = FEt_beta3,
    se3 = FEt_SE_beta3
  ),
  
  FElT = list(
    title = "Location FE (w/ Trend)",
    beta1 = FElT_beta1,
    se1 = FElT_SE_beta1,
    beta2 = FElT_beta2,
    se2 = FElT_SE_beta2,
    beta3 = FElT_beta3,
    se3 = FElT_SE_beta3
  )
)


### 2B. Define Discrepancy Range ###

discrepancy_grid <- seq(
  -1,
  1,
  by = 0.02
)


### 2C. Define Alpha Levels ###

# Scenario-specific alpha used to address Large-N problem
adjusted_alpha <- alphaLevel

# Conventional benchmark
fixed_alpha <- 0.05


### 2D. Define Output Directory ###

power_dir <- file.path(
  directory,
  "Power Curves"
)

dir.create(
  power_dir,
  recursive = TRUE,
  showWarnings = FALSE
)



##################################
### SET UP - Traditional Alpha ###
##################################
save_empirical_power_figure(
  model_results,
  1,
  b1True,
  fixed_alpha,
  "Fixed Alpha = 0.05",
  power_dir,
  paste0("Power_Beta1_FixedAlpha005_(Case", case,")"),
  discrepancy_grid
)

save_empirical_power_figure(
  model_results,
  2,
  b2True,
  fixed_alpha,
  "Fixed Alpha = 0.05",
  power_dir,
  paste0("Power_Beta2_FixedAlpha005_(Case", case,")"),
  discrepancy_grid
)

save_empirical_power_figure(
  model_results,
  3,
  b3True,
  fixed_alpha,
  "Fixed Alpha = 0.05",
  power_dir,
  paste0("Power_Beta3_FixedAlpha005_(Case", case,")"),
  discrepancy_grid
)

###############################
### SET UP - Adjusted Alpha ###
###############################
save_empirical_power_figure(
  model_results,
  1,
  b1True,
  adjusted_alpha,
  "Adjusted Alpha",
  power_dir,
  paste0("Power_Beta1_AdjustedAlpha_(Case", case,")"),
  discrepancy_grid
)

save_empirical_power_figure(
  model_results,
  2,
  b2True,
  adjusted_alpha,
  "Adjusted Alpha",
  power_dir,
  paste0("Power_Beta2_AdjustedAlpha_(Case", case,")"),
  discrepancy_grid
)

save_empirical_power_figure(
  model_results,
  3,
  b3True,
  adjusted_alpha,
  "Adjusted Alpha",
  power_dir,
  paste0("Power_Beta3_AdjustedAlpha_(Case", case,")"),
  discrepancy_grid
)

# Audio to signal power & exports completed
beep(sound = 5)



##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

### If you need to load Rdata (excluding X & y since they're generated each iteration)

# directory <- "C:/Users/gendronj1/OneDrive - Xavier University/SYNC/School/VT/RESEARCH/Dissertation/CH2 - 3YP/Code/Simulation/Simulation Results - R Files/x3/"
# 
# load(
#   file.path(
#     directory,
#     "DATA-MC_ALL_DATA_(Case1).RData"
#   )
# )


