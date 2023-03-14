#' -SPMBayes-
#'
#' An R-package for conducting Bayesian SPM Analysis
#'
#' @param Y  (J x T) array specifying observations for T time points.
#' @param group (J x 1) factor vector specify the groups to which each observation belongs.
#' @param time (T x 1)  vector specify the times (% stance)
#' @param threshold scalar specifying the posterior probability threshold.
#' @param Q.threshold scalar specifying a less conservative threshold using the Q-value
#' @param Hypothesis for specifying which alternative we are interested in, between "alt" and "null".
#' @param paired logical only used for ttestBF to specify if observations are paired
#' @examples
#' Y = mvtnorm::rmvnorm(40,mean=rep(0,10))
#' grp=rep(c(1,2),each=20)
#' time=1:10
#' out <- SPMBayes(Y,group=grp, time, threshold=0.95, Q.threshold=0.05, Hypothesis="alt")
#' plot(x=out, summary = TRUE) ## Get summary plot
#' plot(x=out, summary = FALSE) ## Get SPM Map plot
#'
#' @export
SPMBayes <- function(Y,group,time, threshold=0.95, Q.threshold=0.05, Hypothesis="alt", paired = FALSE){

  ## Summaries of the data for preliminary plot
  MEAN   = t(apply(Y, 2,  \(i) {agg <- stats::aggregate(i, list(group), FUN=mean, na.rm=T)
  agg.out <- agg[,2]; names(agg.out)=paste0("Mean_",agg[,1]); return(agg.out)}))
  SD     = t(apply(Y, 2,  \(i) {agg <- stats::aggregate(i, list(group), FUN=stats::sd, na.rm=T)
  agg.out <- agg[,2]; names(agg.out)=paste0("SD_",agg[,1]); return(agg.out)}))

  SUMMARY = data.table::data.table(time=time, MEAN, SD)

  ## Provide SPM posterior maps
  gp <- factor(group)

  ## Two groups
  if(length(unique(group))==2){
    bfi = apply(Y, MARGIN = 2, FUN = \(y){data = data.table::data.table(YY=y, gp=group)
    BayesFactor::ttestBF(formula=YY~gp,data = data, nullInterval = c(-0.20, 0.20), rscale = "medium", paired = paired)})
    # division of bf(!interval vs 0)/bf(interval vs 0)
    bft = lapply(bfi, FUN = function(bfi){bfi[2]/bfi[1]})
    }else{

  bfi = apply(Y, MARGIN = 2, FUN = \(y){data =data.table::data.table(YY=y,gp=gp)
  BayesFactor::anovaBF(YY~gp,data = data, progress = FALSE)}) ## Get the Bayes factor for comparing the groups at each time point.

  # division of bf(!interval vs 0)/bf(interval vs 0); but already done in anovaBF
  bft = bfi
    }

  # use a prior odds of 1 (equal probability for H0 and H1)
  prior.odds = lapply(bft, FUN = function(bft){BayesFactor::newPriorOdds(bft, type = "equal")})
  # convert the prior odds to posterior odds and then to posterior probability
  post.odds = mapply(FUN = function(prior.odds, bft){prior.odds*bft}, prior.odds, bft)
  post.prob = lapply(post.odds, FUN = function(post.odds){BayesFactor::as.BFprobability(post.odds)})

  # add the Bayes factors and posterior probabilities to the SUMMARY data frame
  SUMMARY$BF   = sapply(bft, FUN = function(bft){as.numeric(as.vector(bft))})
  SUMMARY$PPH1 = sapply(post.prob, FUN = function(post.prob){as.numeric(as.vector(post.prob))})[1,]
  SUMMARY$PPH0 = sapply(post.prob, FUN = function(post.prob){as.numeric(as.vector(post.prob))})[2,]

  # calculate the Q-values for each time point (Q = cumulative mean of sorted posterior error probabilities)
  PEPH0       = cbind(time, 1 - SUMMARY$PPH0)            # calculate PEP (= 1 - PP) per time point
  PEPH0s      = PEPH0[order(PEPH0[,2]),]                 # sort the PEP in ascending order
  QH0         = cumsum(PEPH0s[,2])/seq_along(PEPH0s[,2]) # Q = cumulative mean of sorted PEP values
  QH0         = cbind(PEPH0s,QH0)    # attach Q-values to the PEP and time variable
  QH0s        = QH0[order(QH0[,1]),] # order the Q-values back according to the time variable
  SUMMARY$QH0 = QH0s[,3]             # add Q-values to the SUMMARY dataframe
  PEPH1       = cbind(time, 1 - SUMMARY$PPH1)
  PEPH1s      = PEPH1[order(PEPH1[,2]),]
  QH1         = cumsum(PEPH1s[,2])/seq_along(PEPH1s[,2])
  QH1         = cbind(PEPH1s,QH1)
  QH1s        = QH1[order(QH1[,1]),]
  SUMMARY$QH1 = QH1s[,3]

  # cluster identification: where is the evidence strong enough to accept H0 or H1?
  # conservative threshold: posterior error probability <= 0.05 (i.e. posterior probability >= 0.95)
  SUMMARY$tH0c = SUMMARY$PPH0 >= threshold
  SUMMARY$tH1c = SUMMARY$PPH1 >= threshold
  # less conservative threshold: Q-value < 0.05
  SUMMARY$tH0  = SUMMARY$QH0 <= Q.threshold
  SUMMARY$tH1  = SUMMARY$QH1 <= Q.threshold

  # PP for the largest Q below 0.05
  if(Hypothesis=="alt"){
    if(sum(SUMMARY$tH1 == 1)>0){
      PPq = min(SUMMARY$PPH1[which(SUMMARY$tH1 == 1)])
    }else{
       # less conservative threshold: Q-value < 0.05
      Q.threshold <- min(SUMMARY$QH1)
      SUMMARY$tH1  = SUMMARY$QH1 <= Q.threshold
      PPq = min(SUMMARY$PPH1[which(SUMMARY$tH1 == 1)])
      cat("No significance at the specified Q-value setting to a Q-value =", min(SUMMARY$QH1),"that has at least 1 significant timepoint")
    }

  }

  if(Hypothesis=="Null"){
    if(sum(SUMMARY$tH1 == 1)>0){
      PPq = min(SUMMARY$PPH1[which(SUMMARY$tH0 == 1)])
    }else{

      # less conservative threshold: Q-value < 0.05
      Q.threshold <- min(SUMMARY$QH0)
      SUMMARY$tH0  = SUMMARY$QH0 <= Q.threshold
      PPq = min(SUMMARY$PPH1[which(SUMMARY$tH0 == 1)])
      cat("No significance at the specified Q-value setting to a Q-value =", min(SUMMARY$QH0),"that has at least 1 significant timepoint")
    }

  }

  structure(list(SUMMARY=SUMMARY, PPq=PPq, threshold=threshold,Q.threshold=Q.threshold), class = "SPMBayes")
}


#' @rdname SPMBayes
#' @param x Summary object of class \code{SPMBayes}.
#' @param summary Logical stating whether to only give summary plot with mean and standad deviation
#' @param ... Further arguments passed to ggplot2
#'
#' @export
plot.SPMBayes <- function(x, summary = TRUE,  ...) {
  Mean <- sd <- time <- name <- SD <- PPH1 <- tH1c <- NULL

  SUMMARY <- x$SUMMARY; PPq <- x$PPq; threshold <- x$threshold; Q.threshold <- x$Q.threshold

  dat0Mean <- data.table::data.table(time=SUMMARY$time,SUMMARY[,startsWith(names(SUMMARY), c("Mean")), with=FALSE])
  dat0Sd <- data.table::data.table(time=SUMMARY$time,SUMMARY[,startsWith(names(SUMMARY), c("SD")), with=FALSE])
  dat <- data.table::data.table(tidyr::pivot_longer(dat0Mean,cols=-c(time), values_to = "Mean"), tidyr::pivot_longer(dat0Sd,cols=-c(time), values_to="SD")[,"SD"])

   if(summary){
     ggplot2::ggplot(data = dat, ggplot2::aes(x = time, y = Mean, color=name), linewidth = 1) +
       ggplot2::geom_line() +
       ggplot2::geom_ribbon(ggplot2::aes(x = time, ymin = Mean - SD, ymax = Mean + SD, fill = name),
                   linetype = 2, alpha = 0.3) +
       ggplot2::xlab("time (%)") +
       ggplot2::ggtitle("mean \u00B1 SD time series") +
       ggplot2::scale_x_continuous(expand = c(0,0)) +
       ggplot2::scale_y_continuous(expand = c(0,0)) +
       ggplot2::theme_linedraw()
   }else{
     ggplot2::ggplot(data = SUMMARY, ggplot2::aes(x = time, y = PPH1)) +
       ggplot2::geom_line(ggplot2::aes(y = PPH1), colour = "black", linewidth = 1) +
       ggplot2::geom_hline(yintercept = PPq, linewidth = 0.5, linetype = "dashed") +
       ggplot2::geom_hline(yintercept = threshold, linewidth = 0.5, linetype = "dashed") +
       ggplot2::geom_hline(yintercept = Q.threshold, linewidth = 0.5, linetype = "dashed") +
       ggplot2::annotate("text", x = 10, y = Q.threshold-0.01, label = paste0("Q <", round(Q.threshold,2))) +
       ggplot2::annotate("text", x = 10, y = threshold+0.02, label =  paste("P(H1 | data) > ", threshold, sep="")) +
       ggplot2::annotate("text", x = 10, y = 0.00, label =  paste("P(H0 | data) > ", threshold, sep="")) +
       ggplot2::geom_ribbon(data = subset(SUMMARY, tH1c), ggplot2::aes(ymin = threshold, ymax = PPH1),
                   fill = "gray", alpha = 0.4) +
       ggplot2::geom_ribbon(data = subset(SUMMARY, tH1c), ggplot2::aes(ymin = PPq, ymax = PPH1),
                   fill = "gray", alpha = 0.4) +
       ggplot2::ylim(0,1) +
       ggplot2::ylab("posterior probability") +
       ggplot2::xlab("time (%)") +
       ggplot2::ggtitle("Bayesian SPM: Posterior Probability Map") +
       ggplot2::scale_x_continuous(expand = c(0,0)) +
       ggplot2::theme_linedraw()
   }

}
