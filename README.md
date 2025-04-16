## <em>Each R file does the following:</em>
1. Generates data depending on which case it is (1-12) for either 1, 3 or 5 covariates (depending on the file)
2. Runs the data through all ten meta-regression methodology specifications in a Monte Carlo simulation
3. Generates graphs to illustrate the results for all ten methodology specifications side by side.
4. Calculates power, and then exports all data from the simulation (including the data itself, parameters, standard errors, etc.)


<em>Above is done to test which meta-regression methodology specification performs optimally in the face of joint heterogeneity in the location and time that each study was conducted in.</em>
>NOTE: The paper this code refers to will be available on arXiv very soon.

## Sample Code

This illustrates how the code operates when there is one covariate

#### Set Up - Find True Parameters
```

# Start by setting up variance covariance matrix (presumed same for all studies)
library(mvtnorm)

S00 <- matrix(1, nrow = 1, ncol = 1)                            # top left
S01 <- matrix(0.5, nrow = 1, ncol = 1)                          # top right
S10 <- t(S01)                                                   # bottom left
S11 <- S00                                                      # bottom right

SIGMA <- rbind(cbind(S00, S01), cbind(S10, S11))                # VarCov matrix

# Elements used to calculate true parameters
varY <- S00[1,1]
covY <- SIGMA[2:length(SIGMA[1,]),1:1]
covX <- SIGMA[2:length(SIGMA[1,]),2:length(SIGMA[1,])]   
  # the length() extracts the length of one row
trueBetas <- solve(covX) %*% t(t(covY)) # t(t()) to transpose once, likely because of datatype

### The true parameters
sig2True <- varY - covY %*% trueBetas                           # Omega
b1True <- trueBetas[1,1]

R2True <- 1 - sig2True/varY
```

#### Set Up - Input Simulation Parameters
```
mu_x1 <- 1      # Means
numYears <- 5
case <- 1       # Options: 1-12

# Sample size per study
#
n <- 50
#n <- 100
#n <- 150

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
```

#### Data Generation
```
library(tidyverse)
library(lmtest)
library(plm)
library(lme4)
library(sandwich)

#HOW THE DATA IS GENERATED
  #-Each (country) has 1 study ran per year for 5 years
  #-Heterogeneity in countries is captured by the different mean value for each country
  #-Heterogeneity in time is captured by the trend in y=y+0.5*i

#################
### COUNTRY 1 ###
#################

### 2020 ###
mu_y <- -2

mu_0 <- c(mu_y, mu_x1) # The mean for y changes, x remains mean 1
z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
y <- z_data[, 1]
x <- z_data[, 2]

### Comments here were testing adding a trend variable
X <- data.frame(constant = 1, x = x, year = 2020, country = 1, paper = 1, trend = -1)

### 2021-2024 ###
for (i in 1:(numYears-1)) {
  mu_y <- mu_y+(timeHet*i) # The mean by __ per year
  mu_0 <- c(mu_y, mu_x1) # x remains mean 0
  
  z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
  y2 <- z_data[, 1]
  x <- z_data[, 2]

  ### Comments here were testing adding a trend variable
  X2 <- data.frame(constant = 1, x=x, year = 2020 + i, country = 1, paper = i + 1, trend = -1+i/2)

  X <- bind_rows(X, X2)
  y <- c(y, y2)
}

#############################
### COUNTRY 2 - 5 (or 9) ###
#############################

### 2020-2024 ###
paper_iterator <- numYears+1

if(numCountries==5){
  for (i in 2:numCountries) {
    for (j in 0:(numYears-1)) {
      mu_y <- i-3+(timeHet*j) # The mean for y increases by 1 per country, also increases by __ per year
      mu_0 <- c(mu_y, mu_x1) # x remains mean 1
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2]
      y2 <- z_data[, 1]

      ### Comments here were testing adding a trend variable  
      X2 <- data.frame(constant = 1, x = x, year = 2020 + j, country = i, paper = paper_iterator, trend = -1+j/2)
      
      X <- bind_rows(X, X2)
      y <- c(y, y2)
      
      paper_iterator <- paper_iterator + 1
    }
  }
}else if(numCountries==9){
  countryIndex=2
  for (i in seq(from=-1.5, to=2, by=0.5)) { # (-2), -1.5, -1, ..., 2
    for (j in 0:(numYears-1)) {
      mu_y <- i+(timeHet*j) # The mean for y increases by 1 per country, also increases by __ per year
      mu_0 <- c(mu_y, mu_x1) # x remains mean 0
      
      z_data <- rmvnorm(n, mean = mu_0, sigma = SIGMA)
      x <- z_data[, 2]
      y2 <- z_data[, 1]
      
      ### Comments here were testing adding a trend variable
      X2 <- data.frame(constant = 1, x = x, year = 2020 + j, country = i, paper = paper_iterator, trend = -1+j/2)
      
      paper_iterator <- paper_iterator + 1
      X <- bind_rows(X, X2)
      y <- c(y, y2)
    }
    countryIndex=countryIndex+1
  }
}
```

#### Run Models
```
# To ensure NAs don't keep the regressions from running
y <- na.omit(y)
X <- na.omit(X)

### Random Effects (Study Level) ###
RE <- plm(y ~ x, data=X, index=c("paper"), model="random")  #random model
#summary(RE)

REmse <- mean(RE$residuals^2)                # MSE (df adjusted)
REmae <- mean(abs(RE$residuals))             # MAE (df adjusted)
REmpe <- mean(RE$residuals/RE$model$y)       # MPE (df adjusted)
REmape <- abs(REmpe)

REb0 <- coef(RE)[1]
REb1 <- coef(RE)[2]

REmse_x1 <- mean((REb1-b1True)^2)         # MSE (df adjusted)
REmae_x1 <- mean(abs(REb1-b1True))        # MAE (df adjusted)
REmpe_x1 <- mean((b1True-REb1)/b1True)   # MPE (df adjusted)
REmape_x1 <- abs(REmpe_x1)

REsig2 <- sigma(RE)^2

REb0SE <- sqrt(diag(vcov(RE)))[1]
REb1SE <- sqrt(diag(vcov(RE)))[2]

CIlo_REb0 <- REb0 - (1.96*REb0SE)
CIhi_REb0 <- REb0 + (1.96*REb0SE)
CIlo_REb1 <- REb1 - (1.96*REb1SE)
CIhi_REb1 <- REb1 + (1.96*REb1SE)



### Study-level Fixed Effects ###
FE <- lm(y ~ x + as.factor(paper), data=X)  #fixed model
FEresult <- summary(FE)

FEmse <- mean(FE$residuals^2)                # MSE (df adjusted)
FEmae <- mean(abs(FE$residuals))             # MAE (df adjusted)
FEmpe <- mean(FE$residuals/FE$model$y)       # MPE (df adjusted)
FEmape <- abs(FEmpe)

FEb0 <- coef(FE)[1]
FEb1 <- coef(FE)[2]

FEmse_x1 <- mean((FEb1-b1True)^2)         # MSE (df adjusted)
FEmae_x1 <- mean(abs(FEb1-b1True))        # MAE (df adjusted)
FEmpe_x1 <- mean((b1True-FEb1)/b1True)   # MPE (df adjusted)
FEmape_x1 <- abs(FEmpe_x1)

FEsig2 <- sigma(FE)^2

FEb0SE <- sqrt(diag(vcov(FE)))[1]
FEb1SE <- sqrt(diag(vcov(FE)))[2]

CIlo_FEb0 <- FEb0 - (1.96*FEb0SE)
CIhi_FEb0 <- FEb0 + (1.96*FEb0SE)
CIlo_FEb1 <- FEb1 - (1.96*FEb1SE)
CIhi_FEb1 <- FEb1 + (1.96*FEb1SE)



### Two Way Fixed Effects ###
FE2 <- lm(y ~ x + as.factor(country) + as.factor(year), data = X)
FE2result <- summary(FE2)

FE2mse <- mean(FE2$residuals^2)                # MSE (df adjusted)
FE2mae <- mean(abs(FE2$residuals))             # MAE (df adjusted)
FE2mpe <- mean(FE2$residuals/FE2$model$y)       # MPE (df adjusted)
FE2mape <- abs(FE2mpe)

FE2b0 <- coef(FE2)[1]
FE2b1 <- coef(FE2)[2]

FE2mse_x1 <- mean((FE2b1-b1True)^2)         # MSE (df adjusted)
FE2mae_x1 <- mean(abs(FE2b1-b1True))        # MAE (df adjusted)
FE2mpe_x1 <- mean((b1True-FE2b1)/b1True)   # MPE (df adjusted)
FE2mape_x1 <- abs(FE2mpe_x1)

FE2sig2 <- sigma(FE2)^2

FE2b0SE <- sqrt(diag(vcov(FE2)))[1]
FE2b1SE <- sqrt(diag(vcov(FE2)))[2]

CIlo_FE2b0 <- FE2b0 - (1.96*FE2b0SE)
CIhi_FE2b0 <- FE2b0 + (1.96*FE2b0SE)
CIlo_FE2b1 <- FE2b1 - (1.96*FE2b1SE)
CIhi_FE2b1 <- FE2b1 + (1.96*FE2b1SE)



### [LOCATION-level] Random Effects ###
REl <- plm(y ~ x, data=X, index=c("country"), model="random")  #random model
#summary(RE)

REl_mse <- mean(REl$residuals^2)                # MSE (df adjusted)
REl_mae <- mean(abs(REl$residuals))             # MAE (df adjusted)
REl_mpe <- mean(REl$residuals/REl$model$y)       # MPE (df adjusted)
REl_mape <- abs(REl_mpe)

REl_b0 <- coef(REl)[1]
REl_b1 <- coef(REl)[2]

REl_mse_x1 <- mean((REl_b1-b1True)^2)         # MSE (df adjusted)
REl_mae_x1 <- mean(abs(REl_b1-b1True))        # MAE (df adjusted)
REl_mpe_x1 <- mean((b1True-REl_b1)/b1True)   # MPE (df adjusted)
REl_mape_x1 <- abs(REl_mpe_x1)

REl_sig2 <- sigma(REl)^2

REl_b0SE <- sqrt(diag(vcov(REl)))[1]
REl_b1SE <- sqrt(diag(vcov(REl)))[2]

CIlo_REl_b0 <- REl_b0 - (1.96*REl_b0SE)
CIhi_REl_b0 <- REl_b0 + (1.96*REl_b0SE)
CIlo_REl_b1 <- REl_b1 - (1.96*REl_b1SE)
CIhi_REl_b1 <- REl_b1 + (1.96*REl_b1SE)



### [LOCATION-level] Fixed Effects ###
FEl <- lm(y ~ x + as.factor(country), data=X)  #fixed model
FElresult <- summary(FEl)

FEl_mse <- mean(FEl$residuals^2)                # MSE (df adjusted)
FEl_mae <- mean(abs(FEl$residuals))             # MAE (df adjusted)
FEl_mpe <- mean(FEl$residuals/FEl$model$y)       # MPE (df adjusted)
FEl_mape <- abs(FEl_mpe)

FEl_b0 <- coef(FEl)[1]
FEl_b1 <- coef(FEl)[2]

FEl_mse_x1 <- mean((FEl_b1-b1True)^2)         # MSE (df adjusted)
FEl_mae_x1 <- mean(abs(FEl_b1-b1True))        # MAE (df adjusted)
FEl_mpe_x1 <- mean((b1True-FEl_b1)/b1True)   # MPE (df adjusted)
FEl_mape_x1 <- abs(FEl_mpe_x1)

FEl_sig2 <- sigma(FEl)^2

FEl_b0SE <- sqrt(diag(vcov(FEl)))[1]
FEl_b1SE <- sqrt(diag(vcov(FEl)))[2]

CIlo_FEl_b0 <- FEl_b0 - (1.96*FEl_b0SE)
CIhi_FEl_b0 <- FEl_b0 + (1.96*FEl_b0SE)
CIlo_FEl_b1 <- FEl_b1 - (1.96*FEl_b1SE)
CIhi_FEl_b1 <- FEl_b1 + (1.96*FEl_b1SE)



### [TIME-level] Fixed Effects ###
FEt <- lm(y ~ x + as.factor(year), data=X)  #fixed model
FEtresult <- summary(FEt)

FEt_mse <- mean(FEt$residuals^2)                # MSE (df adjusted)
FEt_mae <- mean(abs(FEt$residuals))             # MAE (df adjusted)
FEt_mpe <- mean(FEt$residuals/FEt$model$y)       # MPE (df adjusted)
FEt_mape <- abs(FEt_mpe)

FEt_b0 <- coef(FEt)[1]
FEt_b1 <- coef(FEt)[2]

FEt_mse_x1 <- mean((FEt_b1-b1True)^2)         # MSE (df adjusted)
FEt_mae_x1 <- mean(abs(FEt_b1-b1True))        # MAE (df adjusted)
FEt_mpe_x1 <- mean((b1True-FEt_b1)/b1True)   # MPE (df adjusted)
FEt_mape_x1 <- abs(FEt_mpe_x1)

FEt_sig2 <- sigma(FEt)^2

FEt_b0SE <- sqrt(diag(vcov(FEt)))[1]
FEt_b1SE <- sqrt(diag(vcov(FEt)))[2]

CIlo_FEt_b0 <- FEt_b0 - (1.96*FEt_b0SE)
CIhi_FEt_b0 <- FEt_b0 + (1.96*FEt_b0SE)
CIlo_FEt_b1 <- FEt_b1 - (1.96*FEt_b1SE)
CIhi_FEt_b1 <- FEt_b1 + (1.96*FEt_b1SE)



### [Study-level] Fixed Effects (with Trend) ###
FEsT <- lm(y ~ x + trend + as.factor(paper), data=X)  #fixed model
FEsTresult <- summary(FEsT)

FEsT_mse <- mean(FEsT$residuals^2)                # MSE (df adjusted)
FEsT_mae <- mean(abs(FEsT$residuals))             # MAE (df adjusted)
FEsT_mpe <- mean(FEsT$residuals/FEsT$model$y)       # MPE (df adjusted)
FEsT_mape <- abs(FEsT_mpe)

FEsT_b0 <- coef(FEsT)[1]
FEsT_b1 <- coef(FEsT)[2]
FEsT_bTrend <- coef(FEsT)[3]

FEsT_mse_x1 <- mean((FEsT_b1-b1True)^2)         # MSE (df adjusted)
FEsT_mae_x1 <- mean(abs(FEsT_b1-b1True))        # MAE (df adjusted)
FEsT_mpe_x1 <- mean((b1True-FEsT_b1)/b1True)   # MPE (df adjusted)
FEsT_mape_x1 <- abs(FEsT_mpe_x1)
FEsT_mse_trend <- mean((FEsT_bTrend-timeHet)^2)         # MSE (df adjusted)
FEsT_mae_trend <- mean(abs(FEsT_bTrend-timeHet))        # MAE (df adjusted)
FEsT_mpe_trend <- mean((timeHet-FEsT_bTrend)/timeHet)   # MPE (df adjusted)
FEsT_mape_trend <- abs(FEsT_mpe_trend)

FEsT_sig2 <- sigma(FEsT)^2

FEsT_b0SE <- sqrt(diag(vcov(FEsT)))[1]
FEsT_b1SE <- sqrt(diag(vcov(FEsT)))[2]
FEsT_bTrendSE <- sqrt(diag(vcov(FEsT)))[3]

CIlo_FEsT_b0 <- FEsT_b0 - (1.96*FEsT_b0SE)
CIhi_FEsT_b0 <- FEsT_b0 + (1.96*FEsT_b0SE)
CIlo_FEsT_b1 <- FEsT_b1 - (1.96*FEsT_b1SE)
CIhi_FEsT_b1 <- FEsT_b1 + (1.96*FEsT_b1SE)
CIlo_FEsT_bTrend <- FEsT_bTrend - (1.96*FEsT_bTrendSE)
CIhi_FEsT_bTrend <- FEsT_bTrend + (1.96*FEsT_bTrendSE)  



### [Location-level] Fixed Effects (with Trend) ###
FElT <- lm(y ~ x + trend + as.factor(country), data=X)  #fixed model
FElTresult <- summary(FElT)

FElT_mse <- mean(FElT$residuals^2)                # MSE (df adjusted)
FElT_mae <- mean(abs(FElT$residuals))             # MAE (df adjusted)
FElT_mpe <- mean(FElT$residuals/FElT$model$y)       # MPE (df adjusted)
FElT_mape <- abs(FElT_mpe)

FElT_b0 <- coef(FElT)[1]
FElT_b1 <- coef(FElT)[2]
FElT_bTrend <- coef(FElT)[3]

FElT_mse_x1 <- mean((FElT_b1-b1True)^2)         # MSE (df adjusted)
FElT_mae_x1 <- mean(abs(FElT_b1-b1True))        # MAE (df adjusted)
FElT_mpe_x1 <- mean((b1True-FElT_b1)/b1True)   # MPE (df adjusted)
FElT_mape_x1 <- abs(FElT_mpe_x1)
FElT_mse_trend <- mean((FElT_bTrend-timeHet)^2)         # MSE (df adjusted)
FElT_mae_trend <- mean(abs(FElT_bTrend-timeHet))        # MAE (df adjusted)
FElT_mpe_trend <- mean((timeHet-FElT_bTrend)/timeHet)   # MPE (df adjusted)
FElT_mape_trend <- abs(FElT_mpe_trend)
    
FElT_sig2 <- sigma(FElT)^2

FElT_b0SE <- sqrt(diag(vcov(FElT)))[1]
FElT_b1SE <- sqrt(diag(vcov(FElT)))[2]
FElT_bTrendSE <- sqrt(diag(vcov(FElT)))[3]

CIlo_FElT_b0 <- FElT_b0 - (1.96*FElT_b0SE)
CIhi_FElT_b0 <- FElT_b0 + (1.96*FElT_b0SE)
CIlo_FElT_b1 <- FElT_b1 - (1.96*FElT_b1SE)
CIhi_FElT_b1 <- FElT_b1 + (1.96*FElT_b1SE)
CIlo_FElT_bTrend <- FElT_bTrend - (1.96*FElT_bTrendSE)
CIhi_FElT_bTrend <- FElT_bTrend + (1.96*FElT_bTrendSE)   



### Mixed Effects (RE @ Study-level) ###
ME <- lmer(y ~ x + (1 | paper), data=X)  #fixed model
MEresult <- summary(ME)

MEmse <- mean(residuals(ME)^2)                # MSE (df adjusted)
MEmae <- mean(abs(residuals(ME)))             # MAE (df adjusted)
MEmpe <- mean(residuals(ME)/y)               # MPE (df adjusted)
MEmape <- abs(MEmpe)

MEb0 <- summary(ME)$coefficients[1]
MEb1 <- summary(ME)$coefficients[2]

ME_mse_x1 <- mean((MEb1-b1True)^2)         # MSE (df adjusted)
ME_mae_x1 <- mean(abs(MEb1-b1True))        # MAE (df adjusted)
ME_mpe_x1 <- mean((b1True-MEb1)/b1True)   # MPE (df adjusted)
ME_mape_x1 <- abs(ME_mpe_x1)

MEsig2 <- sigma(ME)^2

MEb0SE <- summary(ME)$coefficients[3]
MEb1SE <- summary(ME)$coefficients[4]

CIlo_MEb0 <- MEb0 - (1.96*MEb0SE)
CIhi_MEb0 <- MEb0 + (1.96*MEb0SE)
CIlo_MEb1 <- MEb1 - (1.96*MEb1SE)
CIhi_MEb1 <- MEb1 + (1.96*MEb1SE)



### Mixed Effects (RE @ Country-level) ###
MEl <- lmer(y ~ x + (1 | country), data=X)  #fixed model
MElresult <- summary(MEl)

MEl_mse <- mean(residuals(MEl)^2)                        # MSE (df adjusted)
MEl_mae <- mean(abs(residuals(MEl)))             # MAE (df adjusted)
MEl_mpe <- mean(residuals(MEl)/y)                # MPE (df adjusted)
MEl_mape <- abs(MEl_mpe)

MEl_b0 <- summary(MEl)$coefficients[1]
MEl_b1 <- summary(MEl)$coefficients[2]

MEl_mse_x1 <- mean((MEl_b1-b1True)^2)         # MSE (df adjusted)
MEl_mae_x1 <- mean(abs(MEl_b1-b1True))        # MAE (df adjusted)
MEl_mpe_x1 <- mean((b1True-MEl_b1)/b1True)   # MPE (df adjusted)
MEl_mape_x1 <- abs(MEl_mpe_x1)

MEl_sig2 <- sigma(MEl)^2

MEl_b0SE <- summary(MEl)$coefficients[3]
MEl_b1SE <- summary(MEl)$coefficients[4]

CIlo_MEl_b0 <- MEl_b0 - (1.96*MEl_b0SE)
CIhi_MEl_b0 <- MEl_b0 + (1.96*MEl_b0SE)
CIlo_MEl_b1 <- MEl_b1 - (1.96*MEl_b1SE)
CIhi_MEl_b1 <- MEl_b1 + (1.96*MEl_b1SE)
```
