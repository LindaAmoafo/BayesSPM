#!/bin/sh
#SBATCH --job-name=BayesSPM-sourceSlurm
#SBATCH --time=01:10:00

# Model parameters
nsims <- 1e3

devtools::install_github("LindaAmoafo/BayesSPM")
library(BayesSPM)

# Function to simulate data
Sim.Func <- function(i, True_beta=NULL){
  grp=rep(c(0,1,2),each=20)
  grp.matrix <- model.matrix(~factor(grp))
  if(is.null(True_beta)){True_beta <- c(0.5,1.5,2.5)}
  time=1:10

  Y <- matrix(NA, nrow=60, ncol=length(time))
  for(i in 1:length(time)){
    if (i==5|i ==7){
      Y[,i] = rnorm(60, mean=grp.matrix %*% True_beta, sd=0.5)
    }else{Y[,i] = rnorm(60, mean=0, sd=0.5)}
  }
  out <- SPMBayes(Y,group=grp, time, threshold=0.95, Q.threshold=0.05, Hypothesis="alt")

  list(PPH1=summary(out)$PPH1, QH1=summary(out)$QH1)
}

# Approximation
set.seed(12322)
ans <- sapply(1:nsims, \(x) out<- Sim.Func)

saveRDS(ans, "05-sapply.rds")
