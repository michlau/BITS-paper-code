# Analyze how many terms are discarded by the branch-and-bound interaction
# search optimization in BITS
# Supply simulation setting number (4,...,6) as command-line argument

study <- "complex"
# Set number of cores
use_cores <- 32
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE
source("../common/prep.R")

# Define custom BITS model fitting function to obtain the fraction
# of discarded (i.e., skipped) terms
BITS.fun <- function() {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  # Parallel over all repetitions:
  scores.orig <- foreach(i=1:n_reps,
                         .combine=rbind,
                         .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "eval.performance", "E",
                                     "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt = TRUE)
      split_data <- dat$train; split_data_val <- dat$val
      # Combine training and validation data sets as in the final evaluation
      split_data <- rbind(split_data, split_data_val)
      # Split off environmental variable (if it exists) and obtain modes of
      # inheritance for all SNPs
      Z <- NULL -> Z_val
      if(E) {
        Z <- split_data[,"E",drop=FALSE]
        Z_val <- split_data_val[,"E",drop=FALSE]
      }
      if(snps) {
        mois <- BITS::MOI(split_data[,-c(1, rem_E)], split_data[,1])
        split_data_bin <- BITS::applyMOI(split_data[,-c(1, rem_E)], mois)
        split_data_val_bin <- BITS::applyMOI(split_data_val[,-c(1, rem_E)], mois)
      } else {
        split_data_bin <- split_data[,-c(1, rem_E)]
        split_data_val_bin <- split_data_val[,-c(1, rem_E)]
      }

      # Potentially concatenate environmental data back
      if(E) {
        split_data_bin <- cbind(split_data_bin, Z); split_data_val_bin <- cbind(split_data_val_bin, Z_val)
        Z <- NULL -> Z_val
      }

      s <- NULL
      model <- BITS::BITS.complete(split_data_bin, split_data[,1],
                                   max.vars = 3,
                                   boosting.iter = 50, learning.rate = 0.1,
                                   lambda = NULL, alpha = 1, nfolds = 0,
                                   term.select = "elnet",
                                   relax = TRUE, gamma2 = 0.5,
                                   max.iter = function(p) -1, parallel = FALSE,
                                   gmin = function(gmax) 0.001 * gmax, gsteps = 50)
      tmp <- NULL
      for(j in 1:length(model)) {
        tmp <- rbind(tmp, data.frame(Iter = i, gamma = model[[j]]$gamma,
                                     evaluated.terms = model[[j]]$evaluated.terms,
                                     possible.terms = model[[j]]$possible.terms,
                                     boosting.iter = 1:length(model[[j]]$evaluated.terms)))
        # Save gamma, number of evaluated terms, and number of total possible terms.
        # Also record boosting.iter, as algorithm might stop early.
      }
      model <- NULL; gc()
      tmp
    }
  return(list(Score = scores.orig))
}

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  scenes <- as.numeric(args[1])
}
cat(scenes); cat("\n")

for(s in scenes) {
  data.complete <<- data.list[[s]]
  scen <<- s

  res <- data.frame(s = s, BITS.fun())

  load("discarded_terms.RData")
  res[[s]] <- res2
  save(res, file="discarded_terms.RData")
}

stopCluster(cl)
quit("no")

##################################

# Locally process results and create plot

load("discarded_terms.RData")
res2 <- NULL
for(s in 4:6) {
  tmp <- res[[s]]
  tmp <- tmp[abs(tmp$Score.gamma) > 1e-10,]
  res2 <- rbind(res2, tmp)
}
res <- res2
res$discarded.fraction <- 1-res$Score.evaluated.terms/res$Score.possible.terms
scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
res <- cbind(res, scenarios[res$s,])

# Function to get unique double values
# The purpose is simply to unequality due to numeric issues/round errors
real.unique <- function(x) {
  tmp <- sort(unique(x))
  diffs <- diff(tmp)
  c(tmp[which(diffs > 1e-10)], tmp[length(tmp)])
}

# Function for replacing the original gamma values with its corresponding
# approximate unique values.
# This function is written C++ to reduce computation time.
library(Rcpp)
cppFunction('NumericVector unifyGamma(NumericVector gamma, NumericVector unique_gammas) {
  NumericVector gamma2 = clone(gamma);
  for(int i = 0; i < gamma.length(); i++) {
    for(int j = 0; j < unique_gammas.length(); j++) {
      if(fabs(gamma[i] - unique_gammas[j]) < 1e-10) {
        gamma2[i] = unique_gammas[j];
        break;
      }
    }
  }
  return gamma2;
}')

# Scale all gamma values to [1e-3, 1] to enable comparison
# This scaling is performed by creating a log-equidistant gamma-path
# on [1e-3, 1] and mapping the unscaled values to the scaled values
# by their order
res$gamma.scaled <- NA
max.gamma <- 1
min.gamma <- 1e-3
gamma.steps <- exp(seq(log(max.gamma), log(min.gamma), length.out = 50))
for(s in unique(res$s)) {
  for(i in unique(res$Score.Iter)) {
    tmp.res <- res[res$s == s & res$Score.Iter == i,]
    tmp.gamma <- tmp.res$Score.gamma

    unique.gamma <- sort(real.unique(tmp.gamma), decreasing = TRUE)
    tmp.gamma2 <- unifyGamma(tmp.gamma, unique.gamma)
    tmp <- match(tmp.gamma2, unique.gamma)
    gamma2 <- gamma.steps[tmp]

    res$gamma.scaled[res$s == s & res$Score.Iter == i] <- gamma2
  }
}

#########################################

# Plot

res2 <- res
res2 <- res2[!is.na(res2$discarded.fraction),]
res2$snr.factor <- factor(paste("SNR =", res2$snr), levels = unique(paste("SNR =", res2$snr)))
res2 <- res2[log(res2$gamma.scaled) > -4,]
res2$gamma.factor <- factor(res2$gamma.scaled, levels=unique(res2$gamma.scaled))
library(ggplot2)
p <- ggplot(res2, aes(x=Score.boosting.iter, y=discarded.fraction, color = log(gamma.scaled), group = gamma.factor)) +
  geom_smooth(method = "loess", se=FALSE) +
  theme_bw() +
  xlab("Boosting Iteration") + ylab("Fraction of Discarded Interaction Terms") +
  scale_colour_gradientn(name = expression(log(gamma)), colours = rev(hcl.colors(7))) +
  facet_grid(cols = vars(snr.factor))
p
ggsave(plot = p, filename = "discarded_terms.pdf", width = 8, height = 3)
