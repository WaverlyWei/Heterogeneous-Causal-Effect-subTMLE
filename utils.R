subTMLE <- function(D,Y,X,X1,C){
  # initial e
  psModel <- glm(D~., data=newdat[,-Y], family="binomial")
  e1 <- predict(psModel,data = data.frame(D=1,X),type="response")
  e0 <- predict(psModel,data = data.frame(D=0,X),type="response")
  # initial mu 
  muModel <- glm(Y~., data=newdat, family="binomial")
  mu <- predict(muModel,data = newdat,type="response")
  mu1 <- predict(muModel,data = data.frame(Y,D=1,X),type="response")
  mu0 <- predict(muModel,data = data.frame(Y,D=0,X),type="response")
  initMu <- cbind(mu,mu0,mu1)
  names(initMu) <- c("mu", "mu0", "mu1")
  # estimate CDF of 20% PRS
  f <- ecdf(X1) 
  P <- f(C)
  # compute perturbing direction
  S1 <- (X1<=C)/P * D/e1
  S0 <- (X1<=C)/P * (1-D)/e0
  # compute epsilon
  epsilon <- coef(glm(Y~-1 + offset(mu) + S1 + S0, 
                      family="binomial"))
  
  mustar <- initMu+ c((epsilon[1]*S0 + epsilon[2]*S1), 
                      epsilon[1]/e0, 
                      epsilon[2]/e1)
  
  mustar <- plogis(mustar)
  # calc parameters
  mu1 <- mean(mustar[,"mu1"])
  mu0 <- mean(mustar[,"mu0"])
  #psiATE <- mu1-mu0
  psi <- mu1/mu0
  # influence function 
  IC.ATE <- (X1<=C)/P*(D/e1 - (1-D)/e0)*(Y-mustar[,"mu"]) +
    mustar[,"mu1"] - mustar[,"mu0"] - (mu1-mu0)
  # var 
  var.psi <- var(IC.ATE)/length(Y)
  # CI
  upperCI <- psi + 1.96*sqrt(var.psi)
  lowerCI <- psi - 1.96*sqrt(var.psi)
  # p value
  pvalue <- 2*pnorm(-abs(psi/sqrt(var.psi)))
  
  return(c(psi,upperCI,
             lowerCI,var.psi,
             pvalue))
}
