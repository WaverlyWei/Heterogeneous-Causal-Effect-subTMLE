library(dplyr)
library(tidyverse)
library(glmnet)
library(caret)

source("/Users/waverlywei/Desktop/streaming_data/streaming data project/data/utils.R")
## load data
dat <- readRDS("processed_CAD_dat_for_Jingshen_TMLE.rds")

## data cleaning ##

# remove col 2:16(lifestyles) except col8 (oily fish intake)
# remove 17:23(family history)
# remove 101:110(unused PRS)
dat <- dat %>% select(-c(2:7,9:16,17:23,101:110))
# remove chest pain
dat <- dat %>% select(-chest_pain_angia)


## data splitting ##
train_dat <- dat %>% filter(Phase==1)
test_dat <- dat %>% filter(Phase==2)

# remove phase variable 
train_dat <- train_dat %>% select(-Phase)
test_dat <- test_dat %>% select(-Phase)

## Check for data balancing ##
sum(train_dat$cad==1)/nrow(train_dat)
#[1] 0.05993409
sum(test_dat$cad==1)/nrow(test_dat)
#[1] 0.05398269

## create interaction variable ##

# all the covariates except PRS, intervention
cov <- train_dat %>% select(-LDpred,-oily_fish_intake_f1329_0_0)
interactions <- model.matrix(cad ~ .^2, data = cov)
# remove interaction term
interactions <- interactions[,-1]

## Lasso 
# fit by cross validation 
fit <- cv.glmnet(interactions, train_dat$cad,family = "binomial")
#save(fit,file = "cvLassoFit.RData")
plot(fit)
# min lambda
fit$lambda.min
# extract coefficient
lassocoef <- coef(fit, s = "lambda.min")
lassocoef <- as.matrix(lassocoef)
varIdx <- which(lassocoef[,1]!=0)
varIdx <- varIdx[-1]
#save(varIdx,file = "varIndex.RData")
# selected 14 variables 

## work on test set
test_cov <- test_dat %>% select(-LDpred,-oily_fish_intake_f1329_0_0)
test_interactions <- model.matrix(cad ~ .^2, data = test_cov)
test_interactions <- test_interactions[,-1]

## predict on test set
# remove intercept 
pred <- predict(fit,test_interactions,s = "lambda.min",type="response")

class_pred <-ifelse(pred > 0.50,1,0)
class_actual <-ifelse(test_dat$cad > 0.50,1,0)


# check prediction performance
cfmat <- confusionMatrix(factor(class_pred),factor(class_actual))
cfmat
#save(cfmat,file = "confusionMat.RData")

## Select covariates ##
X <- test_interactions[,varIdx]
Y <- test_dat$cad
D <- dat %>% filter(Phase==2) %>% 
  select(oily_fish_intake_f1329_0_0) 
D <- as.matrix(D)

# code D to binary by mean
# >=mean: 1, < mean = 0
D <- ifelse(D>=mean(D),1,0)
#[1] 0.5651583

# conditioning variable: polygenic risk score
X1 <- dat %>% filter(Phase==2) %>% 
  select(LDpred) 
X1 <- as.matrix(X1)
# compute quantile threshold
C <- quantile(X1,seq(0.1,1,by = 0.1))


## Single Subgroup TMLE ##
newdat <- as.data.frame(cbind(Y,D,X))
names(newdat)[1:2] <- c("Y","D")

res <- matrix(NA,nrow = 10,ncol = 5)
for(i in 1:length(C)){
  res[i,] <- subTMLE(D,Y,X,X1,C[i])
}

res <- as.data.frame(res)
colnames(res) <- c("OddsRatio","upperCI","lowerCI","Variance","pval")


## Plot
ggplot(data=res, aes(x = seq(0.1,1,0.1), y = OddsRatio))+
  geom_point(shape = 16,color = "blue",size = 2) +
  geom_line(color = "blue",size = 1)+
  geom_ribbon(data=res,aes(ymin=lowerCI,
                              ymax=upperCI),
              alpha = 0.4,fill = "blue")+
  theme_bw()+
  scale_x_continuous(breaks = seq(0.1, 1, by = 0.1))+
  theme(text = element_text(size=28),
        plot.title = element_text(hjust = 0.5,size = 18,face = "bold"),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +xlab("PRS Percentile")+
  geom_hline(yintercept = 1, linetype = "dashed",size = 1.5)+
  ggtitle("Odds Ratio of Fish Oil Intake Above Average vs. Below Average")


## initial e
# psModel <- glm(D~., data=newdat[,-Y], family="binomial")
# e1 <- predict(psModel,data = data.frame(D=1,X),type="response")
# e0 <- predict(psModel,data = data.frame(D=0,X),type="response")
# # initial mu 
# muModel <- glm(Y~., data=newdat, family="binomial")
# mu <- predict(muModel,data = newdat,type="response")
# mu1 <- predict(muModel,data = data.frame(Y,D=1,X),type="response")
# mu0 <- predict(muModel,data = data.frame(Y,D=0,X),type="response")
# initMu <- cbind(mu,mu0,mu1)
# names(initMu) <- c("mu", "mu0", "mu1")
# # estimate CDF of 20% PRS
# f <- ecdf(X1) 
# P <- f(prs20)
# # compute perturbing direction
# S1 <- (X1<=prs20)/P * D/e1
# S0 <- (X1<=prs20)/P * (1-D)/e0
# # compute epsilon
# epsilon <- coef(glm(Y~-1 + offset(mu) + S1 + S0, 
#                     family="binomial"))
# 
# mustar <- initMu+ c((epsilon[1]*S0 + epsilon[2]*S1), 
#                   epsilon[1]/e0, 
#                   epsilon[2]/e1)
# 
# mustar <- plogis(mustar)
# # calc parameters
# mu1 <- mean(mustar[,"mu1"])
# mu0 <- mean(mustar[,"mu0"])
# psiATE <- mu1-mu0
# psiOdds <- mu1/mu0
# # influence function 
# IC.ATE <- (X1<=prs80)/P*(D/e1 - (1-D)/e0)*(Y-mustar[,"mu"]) +
#   mustar[,"mu1"] - mustar[,"mu0"] - (mu1-mu0)
# # var 
# var.psi <- var(IC.ATE)/length(Y)
# # CI
# upperCI <- psi + 1.96*sqrt(var.psi)
# lowerCI <- psi - 1.96*sqrt(var.psi)
# # p value
# pvalue <- 2*pnorm(-abs(psi/sqrt(var.psi)))
# 
