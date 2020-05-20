## load libraries
library(nimble)
library(coda)
library(mcmcplots)

source("Dist_Siler.R")

CH<-readRDS("CP.CJS_array.rds")
age<-readRDS("CP.age_array.rds")
tKD<-readRDS("CP.dead.rds")
CP.CVdata<-readRDS("CP.CVdata.rds")
f<-readRDS("CP.f.rds")

CH[260,53]<-1
CH[99,80]<-1

## read in data
tB <- CP.CVdata$birth_yr
names(tB) <- NULL

## extract max possible death time
tM <- ifelse(is.na(tKD), ncol(CH), tKD)

## extract last alive time
tL <- apply(CH, 1, function(x) max(which(x == 1)))

## normalise to survival times
## (necessary at the moment due to censoring
## constraints)
tM <- tM - tB
tL <- tL - tB
tKD <- tKD - tB

## define censoring matrices
cint <- cbind(tL, tM)
censored <- ifelse(!is.na(tKD), 1, 2)
tD <- rep(NA, length(tKD))
dind <- rep(1, length(tKD))

## extract number of captures
y <- apply(CH, 1, sum)
names(y) <- NULL

## set up nind
nind <- length(y)

## set up initial values
tinit <- apply(cint, 1, function(x) {
    runif(1, x[1], x[2])
})


## code for NIMBLE model with censoring
CJS.code <- nimbleCode({

    ## survival components for dead badgers
    for (i in 1:nind) {
      
        ## likelihood for interval-truncated Siler
        censored[i] ~ dinterval(tD[i], cint[i, ])
        tD[i] ~ dsiler(a1, a2, b1, b2, c)

        ## sampling component
        pd[i] <- exp(y[i] * log(mean.p) + (min(floor(tD[i]), cint[i, 2]) - y[i]) * log(1 - mean.p))
        dind[i] ~ dbern(pd[i])
    }

    ## priors
    a1 ~ dexp(1)
    a2 ~ dexp(1)
    b1 ~ dexp(1)
    b2 ~ dexp(1)
    c ~ dexp(1)
    mean.p ~ dunif(0, 1)
    
})

## set up other components of model
CJS.Consts <- list(nind = nind)
CJS.data <- list(y = y, cint = cint, 
    censored = censored, tD = tD, dind = dind)
CJS.inits <- list(
    tD = tinit,
    a1 = 0.1, 
    a2 = 0.1, 
    b1 = 0.1,
    b2 = 0.1,
    mean.p = runif(1, 0, 1),
    c = 0.05
)

## define the model, data, inits and constants
CJSModel <- nimbleModel(code = CJS.code, constants = CJS.Consts, data = CJS.data, inits = CJS.inits, name = "CJS")

## compile the model
cCJSModel <- compileNimble(CJSModel, showCompilerOutput = TRUE)

## try with adaptive slice sampler
CJSconfig <- configureMCMC(cCJSModel, monitors = c("a1", "a2", "b1", "b2", "c", "mean.p"), thin = 1)
CJSconfig$removeSamplers(c("a1", "a2", "b1", "c","b2"))
CJSconfig$addSampler(target = c("a1", "a2", "b1", "c", "b2"), type = 'AF_slice')

#Check monitors and samplers
CJSconfig$printMonitors()
CJSconfig$printSamplers(c("a1", "a2", "b1", "b2", "c"))

#Build the model
CJSbuilt <- buildMCMC(CJSconfig)
cCJSbuilt <- compileNimble(CJSbuilt)

#Run the model
system.time(runAF <- runMCMC(cCJSbuilt, 
    niter = 10000, 
    nburnin = 2000, 
    nchains = 2, 
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))
runAF$summary

#Plot mcmcm
samples.uninf <- runAF$samples
mcmcplot(samples.uninf)
png("traceAF%d.png")
plot(samples)
dev.off()

png("pairsAF%d.png")
pairs(samples)
dev.off()

## save samples
saveRDS(samples, "newSiler.rds")

#Set age variable (quarter years)
x <- 0:80

## extract samples
samples <- as.matrix(samples)[, 1:5]

#Siler Mortality rate
mort <- apply(samples, 1, function(pars, x) {
    ## extract pars
    a1 <- pars[1]
    a2 <- pars[2]
    b1 <- pars[3]
    b2 <- pars[4]
    c <- pars[5]
    
    ## return predictions
    exp(a1-(b1*x)) + c + exp(a2+(b2*x))
}, x = x)

## extract mean and 95% intervals
mort <- apply(mort, 1, function(x) {
    c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

#Siler survival function
surv <- apply(samples, 1, function(pars, x) {
    ## extract pars
    a1 <- pars[1]
    a2 <- pars[2]
    b1 <- pars[3]
    b2 <- pars[4]
    c <- pars[5]
    
    ## return predictions
    exp(((exp(a1))/b1)*(exp(-b1*x)-1) - c*x + ((exp(a2))/b2)*(1-exp(b2*x)))
}, x = x)

## extract mean and 95% intervals
surv <- apply(surv, 1, function(x) {
    c(mean = mean(x), LCI = quantile(x, probs = 0.025), UCI = quantile(x, probs = 0.975))
})

## produce plots
pdf("survcurves.pdf", width = 10, height = 5)

par(mfrow = c(1, 2))

#Draw mortality curve
plot(x, mort[1, ], type = 'l', main = "Hazard function")
lines(x, mort[2, ], lty = 2)
lines(x, mort[3, ], lty = 2)

#Draw survival curve
plot(x, surv[1, ], type = "l", main = "Survivor function")
lines(x, surv[2, ], lty = 2)
lines(x, surv[3, ], lty = 2)

dev.off()
