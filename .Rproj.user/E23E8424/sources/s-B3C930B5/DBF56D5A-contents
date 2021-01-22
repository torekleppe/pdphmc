library(pdphmc)
library(rstan)

Tmax <- 25000

dta <- list(n=11L,vr=0.5^2)


mdl <- build("smile.cpp",lambda = "lambda_arclength", massMatrix = "diagMassISG")
fit.a <- run(mdl,seed=1,Tmax=Tmax,data=dta,chains=10,control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,2.0)))
fit.al <- run(mdl,seed=1,Tmax=Tmax,data=dta,chains=10,control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,10.0)))

mdl.c <- build("smile.cpp",lambda = "lambda_constant", massMatrix = "diagMassISG")
fit.c <- run(mdl.c,seed=1,Tmax=Tmax,data=dta,chains=10,control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,2.0)))
fit.cl <- run(mdl.c,seed=1,Tmax=Tmax,data=dta,chains=10,control=list(absTol=1.0e-3,relTol=1.0e-3,lambdaAdapt=c(1.0,1.0,10.0)))

s.fit <- stan("smile.stan",seed=1,data=dta,chains=10,iter=10000L)


save.image("Computations")
