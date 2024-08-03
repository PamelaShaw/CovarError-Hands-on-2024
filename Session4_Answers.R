library(MASS)
library(mvtnorm) ## for data generation
library(survey) ### need to use sandwich2stage
library(devtools) ### need to install sandwich2stage
library(boot) # for bootstrap SE
library(simex)
library(simexaft)
library(openxlsx) ## to read excel files

#Loading required package: first time need to install from github using this line
## Note this is currently not working for Mac users. Working on a fix.
#install_github("lboe23/sandwich2stage")

library(sandwich2stage)


#### Regression Calibration: Simple linear regression Example 

### Set up some simulated data and fix seed for reproducibility
	set.seed(353)
	n <- 1000  ### Cohort sample size
	n1 <- 100  ### Calibration subset sample size

### anti-logit utility function
	expit=function(x){exp(x)/(1+exp(x))}
	
	xzcorr <- 0.25
	XZ <- rmvnorm(n, mean=c(0,0), sigma=array(c(1,xzcorr,xzcorr,1),dim=c(2,2)))
	x <- XZ[,1]
	z <- XZ[,2]

	y <- rbinom(n,1,expit(-1+log(1.5)*x))	  ### simple logistic regression model
	mean(y)

	reliability <- 0.5	### Set up measurement error with 0.5 Attenuation coef
	sigma_u_sq <- 1/reliability - 1
	xstar1 <- x+rnorm(n, sd=sigma_u_sq^0.5)
	xstar2 <- x+rnorm(n, sd=sigma_u_sq^0.5)

	mydat <- data.frame(y,z,xstar1,xstar2)

####Answers
#1a
	##### 1b Works if you have repeat meausures: just need a second measure with classical measurement error in a subset. 
	###Regression calibration model
	#perform regression calibration
	rcfit<-lm(xstar2~xstar1,data=mydat)
	mydat$xhat<-predict(rcfit,newdata=mydat)
	final.model=glm(y~xhat,family="binomial",data=mydat)

#### Notice the coefficients are similar to those on slide 60
	summary(final.model)
	
###1c: Standard Error Method 1: bootstrap 

bootstrap.function<-function(dat,inds){
	mydat.boot<-dat[inds,]
	rcfit.boot<-lm(xstar2~xstar1,data=mydat.boot)
	mydat.boot$xhat<-predict(rcfit.boot,newdata=mydat.boot)
	
	final.model=glm(y~xhat,family="binomial",data=mydat.boot)
	return(final.model$coef)
}

#### Notice the SE for X is a bit bigger than slide 60

my.boot<-boot(mydat,bootstrap.function,R=1000)
bsSD<- apply(my.boot$t,2,sd)
bsSD

t.stat<-coef(final.model)/bsSD
t.stat

###1d
varU.xstar<-var(mydat$xstar1-mydat$xstar2)/2
varU.xstarAvg<-varU.xstar/2
sd_me<-sqrt(varU.xstarAvg)

###1e
# Fit models
#naive model: ignores error in exposure
## Suppose have xstar1,xstar2 on everyone, where xstar=X+u
mydat$xstarbar<-(mydat$xstar1+mydat$xstar2)/2
naive.model=glm(y~xstarbar,family="binomial",data=mydat)
summary(naive.model)
	
model.simex <-simex(naive.model, SIMEXvariable = "xstarbar", measurement.error = sd_me,asymptotic = FALSE)
#### looks like only jacknife works. maybe becasue there were no other covariates.
summary(model.simex)
par(mfrow=c(2,2))
plot(model.simex,mfrow=c(2,2)) 


## 1f
### true SE
model.simex2 <-simex(naive.model, SIMEXvariable = "xstarbar", measurement.error = sqrt(.5),asymptotic = FALSE)
summary(model.simex2)

### too high SE
model.simex3 <-simex(naive.model, SIMEXvariable = "xstarbar", measurement.error = 1,asymptotic = FALSE)
summary(model.simex3)

### too low SE
model.simex4 <-simex(naive.model, SIMEXvariable = "xstarbar", measurement.error = .25,asymptotic = FALSE)
summary(model.simex4)

###### Question 2
### Inspired by WHI dietary data 
data<-read.xlsx("~/Library/CloudStorage/OneDrive-KaiserPermanente/Meetings/JSM 2024/Rcode/simulated_whidata.xlsx")
names(data)
### true log energy intake = logEn
### selfreported log energy intake = logEn_sr (systematic measurement error)
### biomarker log energy intake = logEn_bm  (classical measurement error)
### repeat measures of these variables have a 2 after them
### Physical activity = step_1000 (units 1000 steps)
### Body mass index = BMI
### age at baseline = age
### race = race
### time = censored event time to first cancer
### delta = observed event indicator

### Note sandwich2stage and simexFit did not seem to like factor variables, so create the binary dummary variables for each race
data$black<-ifelse(data$race=="black",1,0)
data$other<-ifelse(data$race=="other",1,0)
data$white<-ifelse(data$race=="white",1,0)

### turn race into a factor and make the biggest racial group the reference category
data$myrace<-relevel(factor(data$race),ref="white")

### 2a
### Fit calibration
calfit<-lm(logEn_bm~logEn_sr + age+BMI+step_1000+ myrace,data=data)
data$logEn.hat<- predict(calfit,newdata=data)
summary(calfit)

### 2b
true<-coxph(Surv(time,delta)~logEn+age+BMI+step_1000+myrace, data=data)
summary(true)

### 2c
### Fit ignoring meas error
naive<-coxph(Surv(time,delta)~logEn_sr+age+BMI+step_1000+myrace, data=data)
summary(naive)

### 2d
### Perform regression calibration
rc<-coxph(Surv(time,delta)~logEn.hat+age+BMI+step_1000+myrace, data=data)
summary(rc)

### Get bootstrapped SE and confidence intervals
set.seed(806)
bootstrap.function<-function(dat,inds){
    mydat.boot<-dat[inds,]
    rc.model=lm(logEn_bm~ logEn_sr+ age + BMI+step_1000+myrace,data=mydat.boot)
    xhat=predict(rc.model,newdata=mydat.boot)
    final.model=coxph(Surv(time,delta)~logEn.hat+age+BMI+step_1000+myrace, data=mydat.boot)
    return(final.model$coef)
}
data$subset<-ifelse(data$subset==1,1,0)  ### Strata variable needed for bootsrap
my.boot<-boot(data,bootstrap.function, strata=data$subset,R=250) ### R low to speed up for class, I would choose R=500 or 1000
my.boot
### confidence interval for 1st coef... can modify for others
boot.ci(my.boot,type=c("norm","perc"))

bsSD<- apply(my.boot$t,2,sd)
bsSD

### can get sandwich SE
### Declare a simple random sample design for svyglm()
names(data)
datadesign <- svydesign(id=~1, data=data)

## calibration model is the stage 1 model
stage1.model<-svyglm(logEn_bm~logEn_sr+age+BMI+step_1000+black+other, design=datadesign, family=gaussian(),data=data)
datadesign <- update(datadesign, xhat =predict(stage1.model,newdata=datadesign$variables),ID=1:nrow(data))

### outcome model is the stage 2 model
stage2.model<-  svycoxph(Surv(time,delta)~ xhat+age+BMI+step_1000+black+other, design=datadesign)

sandwich.object<-sandwich2stage(stage1.model,stage2.model,xstar="logEn_sr",xhat="xhat",Stage1ID="ID",Stage2ID="ID")
swvar<-vcov(sandwich.object) 
### swvar is 13 by 13 var-covar matrix (3 calibration params and 3 outcome model params )
#swvar 
### These are the SE for the 6 parameters in the stage2.model
sqrt(diag(swvar)[8:13])


### 2e and 2f
##############################Explore SIMEX#############################
#### Simex assumes classical measruement error and a known error function
#### Have have the systematic error variable on everyone. 
#### Have the classical measurment error variable on subset, but very few events on subset
### Have no choice but to apply SIMEX to the systematic error variable
#### Do not expect it to work becasue SIMEX needs classical measurement error
### Note: SIMEX works for glm, need SIMEXAFT for survival models

###fit a AFT model and use simex for accelerated failure time models-- simexaft package
set.seed(120) ## Set random seed so answer is reproducible.
### Dont have repeats of the self reported diet, so dont know measurement error
### Could guess the measurement error variance is = to total variance in exposure
ind <- c("logEn_sr") 
err.mat <- var(data$logEn_sr)

var(data$logEn_sr)
var(data$logEn)


### saving some typing by creating formula variables
true.formula<- "Surv(time, delta) ~ logEn+ age + BMI + step_1000+black+other"
naive.formula<- "Surv(time, delta) ~ logEn_sr + age + BMI + step_1000+black+other"

## Since data are comptuer generated we can consider the true exposure, for a benchmark
 trueAFT<-survreg(formula = formula(true.formula), data = data, dist = "weibull")

#### Naive model using errorprone self-reported energy with no correction for error
naiveAFT<-survreg(formula = formula(naive.formula), data = data, dist = "weibull")

### Look at SIMEX for weibull model, using quadratic exptrapolation
simexFit <- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
summary(simexFit)
plotsimexaft(simexFit,"logEn_sr",ylimit=c(-3,1),extrapolation="quadratic")

### Try linear extrapolation
simexFitLin <- simexaft(formula = formula(naiveAFT), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "linear", dist = "weibull") 
summary(simexFitLin)
plotsimexaft(simexFitLin,"logEn_sr",ylimit=c(-3,1),extrapolation="linear")

### Hazard ratio (HR) estimated by true Cox model
trueCox<-coxph(as.formula(true.formula), data=data)
exp(trueCox$coef)

###  HR estimated by true AFT model, similar to Cox model
exp(-trueAFT$coef/trueAFT$scale)


###Naive
exp(-naiveAFT$coef/naiveAFT$scale)
### HR estimated by simex AFT model: Does terribly because systematic and not classical error in the self-reported energy variable
### Some inflation here, due to the error being more complicated than the simple classical error SIMEX assumes. Also had to guess at error variance.
exp(-simexFit$coef/simexFit$scale)

##########
#### Look at sensitiivty to the assumed error variance.
##########
err.mat2<-.5*err.mat
simexFit2 <- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat2, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
exp(-simexFit2$coef/simexFit2$scale)

##### Could look at actual error variance, but simex error model still wrong.
err.mat3<- var(data$logEn - data$logEn_sr)
simexFit3<- simexaft(formula = formula(naive.formula), data = data, SIMEXvariable = ind, repeated = FALSE, err.mat = err.mat3, B = 50, lambda = seq(0, 2, 0.1),extrapolation = "quadratic", dist = "weibull") 
exp(-simexFit3$coef/simexFit3$scale)





