library(simex)
library(mvtnorm)  


#### Example 1 Contnuous exposure
## Note typo in Seed in code slides 
### Lecture results are from this seed    set.seed(9494)
    sd_me <- 0.6
    n<-100  ### Cohort size
    
# Generate data
    XZ<- rmvnorm(n,c(1,1),sigma=cbind(c(1,.5),c(.5,1)))     
    X<-XZ[,1]
    Z<-XZ[,2]
    Xstar<-X+rnorm(n,0,sd_me)
    Y <- X + 2*Z +rnorm(n, sd= 1)
    
# Fit models
par(mfrow=c(2,2))
    model_true<- lm(Y~X+Z)
    model_naiv <- lm(Y~Xstar+Z,x = TRUE)
    summary(model_true)
    summary(model_naiv)
    model_simex <-simex(model_naiv, SIMEXvariable = "Xstar", measurement.error = sd_me)
    plot(model_simex,mfrow=c(2,2))
    summary(model_simex)
    
    
#### Example 2: MC SIMEX and discrete exposure
## Again  typo in Seed in code in lecture notes slides for this example
### But this code reproduces results in slides
set.seed(4003)
n<-200
alpha<- -.1
beta<-log(1.5)
smoking <- ifelse(runif(n)<.35,1,0)
smokingStar<-ifelse((smoking==1 & (runif(n)<.15)),0,smoking)
pY<-exp(alpha+beta*smoking)/(1+exp(alpha+beta*smoking))
Y<-ifelse(runif(n)<pY,1,0)
smoking<-factor(smoking,levels=c(0,1))
smokingStar<-factor(smokingStar,levels=c(0,1))

Pi <- matrix(c(1,0,0.15,0.85),nrow=2)dimnames(Pi) <- list(c(0,1), c(0,1))
Pi


naive.model <- glm(Y ~ smokingStar , family = binomial, x = TRUE, y = TRUE)
true.model  <- glm(Y ~ smoking , family = binomial)
simex.model <- mcsimex(naive.model, mc.matrix = Pi, SIMEXvariable = "smokingStar")
summary(true.model)
summary(naive.model)
summary(simex.model)

op <- par(mfrow = c(2, 2))
invisible(lapply(simex.model$theta, boxplot, notch = TRUE, outline = FALSE,
    names = c(0.5, 1, 1.5, 2)))
plot(simex.model)

simex.model2 <- refit(simex.model, "line")
plot(simex.model2)
summary(simex.model2)

#### Misspecified misclassification matrix
### note results dont match slides, only due to random seed issues.
Pi2 <- matrix(c(1,0,0.08,0.92),nrow=2)

dimnames(Pi2) <- list(c(0,1), c(0,1))
summary(mcsimex(naive.model, mc.matrix = Pi2, SIMEXvariable = "smokingStar"))

