# Common data loading/splitting, model fitting and prediction, and evaluation
# code used for both simulation study and real data applications
# This file is called 'prep.R', as it loaded as preparation beforehand in other files.

n_reps <- 100 # Number of repetitions
set.seed(123)
train_inds <- list() # Data indices
test_inds <- list()
val_inds <- list()

prepareSNPData <- function(snpData) {
  x <- snpData["x"]$x
  y <- snpData["y"]$y
  preparedSNPData <- data.frame(y = y, x)
  return(preparedSNPData)
}

eval.performance <- function(preds, y) {
  if (!setequal(unique(y), c(0, 1))) {
    max(1-BITS::calcNMSE(preds, y), 0)
  } else {
    if(length(y) > 90000)
      pROC::auc(pROC::roc(y, preds, direction = "<", levels = 0:1))
    else
      logicDT::calcAUC(preds, y) # Fast AUC computation method for small to medium-sized samples
  }
}

# Configure and generate simulated data

if (study == "linear") {
  E <- FALSE
  y_bin <- FALSE
} else if (study == "interactions_no_E") {
  E <- FALSE
  y_bin <- FALSE
} else if (study == "interactions_E") {
  E <- TRUE
  y_bin <- FALSE
} else if (study == "complex") {
  E <- TRUE
  y_bin <- FALSE
}
n.snps <- 50
data.list <- list()
scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
N.test <<- 10000
snps <<- TRUE
complete.search <<- TRUE

lin.pred.list <- list()
theoretical.perf <- numeric()
for(i in 1:nrow(scenarios)) {
  sample.size <- scenarios$N[i]
  sample.size <- scenarios$N[i] + N.test
  snr <- scenarios$snr[i]
  if (study == "linear") {
    # Generate SNPs (predictors) using scrime::simulateSNPglm
    data <- prepareSNPData(scrime::simulateSNPglm(sample.size * n_reps, n.snps, list(c(-1), c(-1, 1)),
                                                  list.snp = list(c(1), c(2, 3)),
                                                  maf = 0.25,
                                                  beta0 = -1.0, beta = c(1.5, 1.5)))
    data[,-1] <- data[,-1] - 1 # 0,1,2 SNP encoding
    # Construct linear predictor used for generating y:
    lin.pred <- data$SNP1 + data$SNP2 +
      (data$SNP3 > 0) + (data$SNP4 > 0) +
      (data$SNP5 == 2) + (data$SNP6 == 2)
    noise.var <- var(lin.pred)/snr

    # Change standard deviation to specific values for different signal intensities
    data$y <- lin.pred + rnorm(sample.size * n_reps, 0, sd = sqrt(noise.var))
    theoretical.perf[i] <- eval.performance(lin.pred, data$y)
    data.list[[i]] <- data
    lin.pred.list[[i]] <- lin.pred
  } else if (study == "interactions_no_E") {
    # Generate SNPs (predictors) using scrime::simulateSNPglm
    data <- prepareSNPData(scrime::simulateSNPglm(sample.size * n_reps, n.snps, list(c(-1), c(-1, 1)),
                                                  list.snp = list(c(1), c(2, 3)),
                                                  maf = 0.25,
                                                  beta0 = -1.0, beta = c(1.5, 1.5)))
    data[,-1] <- data[,-1] - 1 # 0,1,2 SNP encoding
    # Construct linear predictor used for generating y:
    quad <- (sqrt(log(1.5)) + sqrt(log(2)))^2 - log(1.5) - log(2)
    lin.pred <- -0.4 + log(1.5) * (data$SNP1 > 0) + log(2) * (data$SNP2 > 0 & data$SNP3 == 0) +
      quad * (data$SNP1 > 0 & data$SNP2 > 0 & data$SNP3 == 0)
    noise.var <- var(lin.pred)/snr

    # Change standard deviation to specific values for different signal intensities
    data$y <- lin.pred + rnorm(sample.size * n_reps, 0, sd = sqrt(noise.var))
    theoretical.perf[i] <- eval.performance(lin.pred, data$y)
    data.list[[i]] <- data
    lin.pred.list[[i]] <- lin.pred
  } else if (study == "interactions_E") {
    # Generate SNPs (predictors) using scrime::simulateSNPglm
    data <- prepareSNPData(scrime::simulateSNPglm(sample.size * n_reps, n.snps, list(c(-1), c(-1, 1)),
                                                  list.snp = list(c(1), c(2, 3)),
                                                  maf = 0.25,
                                                  beta0 = -1.0, beta = c(1.5, 1.5)))
    data[,-1] <- data[,-1] - 1 # 0,1,2 SNP encoding
    # Generate environmental covariable:
    data$E <- rnorm(sample.size * n_reps, 20, 10)
    data$E[data$E < 0] <- 0
    # Construct linear predictor used for generating y:
    lin.pred <- -0.75 + log(2) * (data$SNP1 > 0) +
      log(4) * data$E/20 * (data$SNP2 > 0 & data$SNP3 == 0)
    # SNR = var(truth)/var(noise) <=> var(noise) = var(truth)/SNR
    noise.var <- var(lin.pred)/snr

    # Change standard deviation to specific values for different signal intensities
    data$y <- lin.pred + rnorm(sample.size * n_reps, 0, sd = sqrt(noise.var))
    theoretical.perf[i] <- eval.performance(lin.pred, data$y)
    data.list[[i]] <- data
    lin.pred.list[[i]] <- lin.pred
  } else if (study == "complex") {
    # Generate SNPs (predictors) using scrime::simulateSNPglm
    data <- prepareSNPData(scrime::simulateSNPglm(sample.size * n_reps, n.snps, list(c(-1), c(-1, 1)),
                                                  list.snp = list(c(1), c(2, 3)),
                                                  maf = 0.25,
                                                  beta0 = -1.0, beta = c(1.5, 1.5)))
    data[,-1] <- data[,-1] - 1 # 0,1,2 SNP encoding
    # Generate environmental covariable:
    data$E <- rnorm(sample.size * n_reps, 20, 10)
    data$E[data$E < 0] <- 0
    # Construct linear predictor used for generating y:
    lin.pred <- (data$SNP1 > 0) + 1.5 * data$E/(1.349*10) * (data$SNP2 > 0) -
      (data$SNP3 == 2) - data$SNP4 +
      2 * data$E/(1.349*10) * (data$SNP5 > 0) * data$SNP6 -
      2 * (data$SNP7 > 0) * (data$SNP8 == 0) * data$SNP9
    # SNR = var(truth)/var(noise) <=> var(noise) = var(truth)/SNR
    noise.var <- var(lin.pred)/snr

    # Change standard deviation to specific values for different signal intensities
    data$y <- lin.pred + rnorm(sample.size * n_reps, 0, sd = sqrt(noise.var))
    theoretical.perf[i] <- eval.performance(lin.pred, data$y)
    data.list[[i]] <- data
    lin.pred.list[[i]] <- lin.pred
  }
}

# Get data for one repetition, splitted into train, validation, and test data
# If final evaluation is performed (and no hyperparameter optimization, i.e., hyperopt=FALSE),
# train and validation data sets are combined
getData <- function(i, scen, hyperopt) {
  N.test <- 10000
  sample.size <- scenarios$N[scen] + N.test
  train.part <- scenarios$N[scen]
  val.part <- floor(0.25*train.part)
  train_inds <- (1+sample.size*(i-1)):(train.part+sample.size*(i-1))
  test_inds <- (train.part+1+sample.size*(i-1)):(sample.size*i)
  val_inds <- 1:val.part
  split_data <- data.complete[train_inds,]
  if(hyperopt) {
    split_data_val <- split_data[val_inds,]
    split_data <- split_data[-val_inds,]
  } else {
    split_data_val <- data.complete[test_inds,]
  }
  return(list(train = split_data, val = split_data_val,
              train_inds = train_inds, test_inds = test_inds, val_inds = val_inds))
}

splitSNPs <- function(dat) {
  if(snps) logicDT::splitSNPs(dat) else dat
}

set.seed(NULL)

params_file <- paste(study, "_params.RData", sep="")
scores_file <- paste(study, "_scores.RData", sep="")
terms_file <- paste(study, "_terms.RData", sep="")

##########################################

# Set up parallel computation

library(foreach)
library(doParallel)
cores <- detectCores()
os <- unname(Sys.info()["sysname"])
if(!exists("force_cores")) force_cores <- FALSE
use_cores <- ifelse(os == "Linux" && !force_cores, cores, use_cores)
cl <- makeCluster(use_cores, outfile="")

registerDoParallel(cl)

##########################################

# Individual model fitting and prediction functions

BITS.fun <- function(dummy, gamma = rep(-1, n_reps), lambda = rep(-1, n_reps)) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "E",
                                "snps", "complete.search")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

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

      # Complete search for simulation study and SALIA application
      # Otherwise, hybrid search for Dorothea
      if(complete.search) {
        max.iter <- function(p) -1
      } else {
        max.iter <- function(p) 3*4*p
      }

      # Potentially concatenate environmental data back
      if(E) {
        split_data_bin <- cbind(split_data_bin, Z); split_data_val_bin <- cbind(split_data_val_bin, Z_val)
        Z <- NULL -> Z_val
      }

      if(hyperopt) s <- NULL else s <- lambda[i] # Term inclusion penalty initially not set
      if(hyperopt) {
        model <- BITS::BITS.complete(split_data_bin, split_data[,1],
                                     max.vars = 3,
                                     boosting.iter = 50, learning.rate = 0.1,
                                     lambda = NULL, alpha = 1, nfolds = 0,
                                     term.select = "elnet",
                                     relax = TRUE, gamma2 = 0.5,
                                     max.iter = max.iter, parallel = FALSE,
                                     gmin = function(gmax) 0.001 * gmax, gsteps = 50)
        # Get optimal gamma and lambda penalties
        best.mod <- BITS::get.ideal.model(model, split_data_val_bin, split_data_val[,1], choose = "min")
        g <- best.mod$best.g; s <- best.mod$best.s
        score <- min(best.mod$val.res$score)
      } else {
        model <- BITS::BITS(split_data_bin, split_data[,1],
                            max.vars = 3, gamma = gamma[i],
                            boosting.iter = 50, learning.rate = 0.1,
                            lambda = NULL, alpha = 1, nfolds = 0,
                            term.select = "elnet",
                            relax = TRUE, gamma2 = 0.5,
                            max.iter = max.iter)
        model$s <- s # Set optimal term inclusion penalty obtained in hyperparameter optimization
        # Evaluate model performance
        score <- eval.performance(predict(model, split_data_val_bin), split_data_val[,1])
        g <- gamma[i]
      }
      model$set_vars <- NULL; gc() # Explicitly clean up, as otherwise parallel foreach might not do it
      # Return iteration, score, and implicitly obtained gamma and lambda penalties
      data.frame(Iter = i, Score = score, lambda = s, gamma = g)
    }
  return(list(Score = scores))
}

glinternet.fun <- function(dummy, lambda = rep(-1, n_reps)) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "E", "y_bin",
                                "optimal.lambda", "calcScorePerObservation", "calcScore",
                                "splitSNPs", "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split SNPs into dominant and recessive to not lose any information
      split_data_bin <- splitSNPs(split_data[,-c(1, rem_E)])
      split_data_val_bin <- splitSNPs(split_data_val[,-c(1, rem_E)])
      train.dat <- split_data_bin; val.dat <- split_data_val_bin
      n.levels <- rep(2, ncol(split_data_bin))
      Z <- NULL -> Z_val
      if(E) {
        Z <- split_data[,"E",drop=FALSE]
        Z_val <- split_data_val[,"E",drop=FALSE]
        train.dat <- cbind(train.dat, Z); val.dat <- cbind(val.dat, Z_val)
        n.levels <- c(n.levels, 1)
      }
      s <- NULL
      if(lambda[i] != -1) s <- lambda[i]
      model <- glinternet::glinternet(train.dat, split_data[,1], n.levels, numCores = 1,
                                      family = ifelse(y_bin, "binomial", "gaussian"),
                                      lambda = NULL)
      if(lambda[i] == -1)
        # Obtain optimal lambda penalty for the corresponding validation data set
        s <- optimal.lambda(model, val.dat, split_data_val[,1], choose = "min")$best.s
      else {
        # Get closest lambda in the computed lambda-path to the optimal lambda
        ll <- model$lambda
        ind <- which.min(abs(ll - s))
        s <- ll[ind]
      }
      cat(sprintf("Iteration %d done\n", i), sep = "")
      # Return iteration, score, and implicitly obtained lambda penalty
      data.frame(Iter = i, Score = eval.performance(predict(model, val.dat, type = "response", lambda = s), split_data_val[,1]), lambda = s)
    }
  return(list(Score = scores))
}

# Score calculation methods for optimal lambda derivation
calcScore <- function(preds, y, scoring_rule) {
  # Deviance for binary y, MSE for continuous y
  # Both measures are equivalent to considering the respective log-likelihoods
  # Deviance is used here instead of AUC, as it can be computed observation-wise
  # which would be needed for the 1SE rule
  if(scoring_rule == "deviance") {
    score <- logicDT::calcDev(preds, y)
  } else if(scoring_rule == "mse") {
    score <- logicDT::calcMSE(preds, y)
  }
  score
}
calcScorePerObservation <- function(preds, y, scoring_rule) {
  scores <- vector(mode = "numeric", length = length(preds))
  for(i in 1:length(preds)) scores[i] <- calcScore(preds[i], y[i], scoring_rule)
  scores
}
optimal.lambda <- function(model, X, y, scoring_rule = "auc", choose = "min") {
  if(inherits(model, "glinternet")) {
    lambda <- model$lambda
    preds <- predict(model, X, type="response")
  } else if(inherits(model, "sprinter")) {
    lambda <- model$lambda
    preds <- predict.sprinter(model, X)
  } else if(inherits(model, "pre")) {
    lambda <- model$glmnet.fit$lambda
    preds <- sapply(lambda, function(l) predict(model, X, penalty.par.val = l, type = "response"))
  }
  y_bin <- setequal(unique(y), c(0, 1))
  if(!y_bin) scoring_rule <- "mse"
  per.observation <- !y_bin || scoring_rule == "deviance"
  if(y_bin) y <- as.integer(y)

  val.res <- data.frame()
  for(i in 1:length(lambda)) {
    if(per.observation) {
      scores <- calcScorePerObservation(as.numeric(preds[,i]), y, scoring_rule)
      m <- mean(scores)
      se <- sd(scores)/sqrt(length(scores))
      val.res <- rbind(val.res, data.frame(s = lambda[i], score = m, se = se, score.plus.1se = m + se))
    } else {
      m <- 1 - logicDT::calcAUC(as.numeric(preds[,i]), y)
      val.res <- rbind(val.res, data.frame(s = lambda[i], score = m))
    }
  }

  min.ind <- min(which(val.res$score == min(val.res$score)))
  if(choose == "1se") {
    max.val <- val.res$score.plus.1se[min.ind]
    min.ind <- min(which(val.res$score <= max.val))
  }

  best.s <- lambda[min.ind]
  return(list(val.res = val.res, best.s = best.s))
}

rulefit.fun <- function(subsample.frac, nodesize, lambda = rep(-1, n_reps)) {
  fam <- ifelse(y_bin, "binomial", "gaussian")
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "y_bin",
                                "optimal.lambda", "calcScorePerObservation", "calcScore")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      if(!is.data.frame(split_data)) {
        split_data <- as.data.frame(split_data); split_data_val <- as.data.frame(split_data_val)
      }
      # Set minbucket/minsplit parameters
      # Also, set maxdepth=3 to allow detecting interactions with order up to 3,
      # as for BITS and MARS
      minbucket <- floor(nodesize[i] * nrow(split_data))
      tree.control <- partykit::ctree_control(maxdepth = 3, mtry = Inf,
                                              minbucket = minbucket, minsplit = 2*minbucket,
                                              testtype = "Univariate", alpha = 0.5)
      colnames(split_data)[1] <- "y"
      if(y_bin) split_data$y <- as.factor(split_data$y)
      model <- pre::pre(y ~ ., data = split_data, family = fam,
                        ntrees = 50, tree.control = tree.control,
                        sampfrac = subsample.frac[i], learnrate = 0.1,
                        verbose = TRUE)

      if(lambda[i] == -1)
        # Obtain optimal lambda penalty for the corresponding validation data set
        s <- optimal.lambda(model, split_data_val[,-1], split_data_val[,1], choose = "min")$best.s
      else
        s <- lambda[i]

      pred <- predict(model, split_data_val[,-1], type = "response", penalty.par.val = s)

      # Return iteration, score, and implicitly obtained lambda penalty
      data.frame(Iter = i, Score = eval.performance(pred, split_data_val[,1]),
                 lambda = as.numeric(s))
    }
  return(list(Score = scores))
}

mars.fun <- function(dummy) {
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "y_bin")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      fam2 <- NULL
      if(y_bin) fam2 <- list(family=binomial())
      model <- earth::earth(split_data[,-1], split_data[,1],
                            degree = 3, glm = fam2)
      # Return iteration, and score
      data.frame(Iter = i, Score = eval.performance(predict(model, split_data_val[,-1], type = "response"), split_data_val[,1]))
    }
  return(list(Score = scores))
}

sprinter.fun <- function(dummy, lambda = rep(-1, n_reps)) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  reps <- 1:n_reps
  # Parallel over all repetitions:
  scores <- foreach(i=reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "E", "y_bin",
                                "optimal.lambda", "calcScorePerObservation", "calcScore", "predict.sprinter",
                                "splitSNPs", "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split SNPs into dominant and recessive to not lose any information
      split_data_bin <- splitSNPs(split_data[,-c(1, rem_E)])
      split_data_val_bin <- splitSNPs(split_data_val[,-c(1, rem_E)])
      train.dat <- split_data_bin; val.dat <- split_data_val_bin
      Z <- NULL -> Z_val
      if(E) {
        Z <- split_data[,"E",drop=FALSE]
        Z_val <- split_data_val[,"E",drop=FALSE]
        train.dat <- cbind(train.dat, Z); val.dat <- cbind(val.dat, Z_val)
      }
      s <- NULL
      if(lambda[i] != -1) s <- matrix(lambda[i])
      model <- sprintr::sprinter(train.dat, split_data[,1], square = FALSE,
                                 lambda = s, num_keep = NULL)
      if(lambda[i] == -1)
        # Obtain optimal lambda penalty for the corresponding validation data set
        s <- optimal.lambda(model, val.dat, split_data_val[,1], choose = "min")$best.s
      cat(sprintf("Iteration %d done\n", i), sep = "")
      preds <- predict.sprinter(model, val.dat)
      if(ncol(preds) > 1) {
        lambda.ind <- which(abs(model$lambda - s) < 1e-10)[1]
        preds <- preds[,lambda.ind]
      }
      # Return iteration, score, and implicitly obtained lambda penalty
      data.frame(Iter = i, Score = eval.performance(as.numeric(preds), split_data_val[,1]), lambda = as.numeric(s))
    }
  return(list(Score = scores))
}

# Custom prediction function, as no predict function is exported for the
# original sprinter object in the sprintr package
predict.sprinter <- function(model, X) {
  lambda <- model$lambda
  b0 <- model$a0
  co <- model$coef
  idx <- model$idx
  preds <- matrix(nrow = nrow(X), ncol = length(lambda))
  X <- as.matrix(X)
  for(i in 1:length(lambda)) {
    cb <- cbind(idx, co[,i])
    main <- cb[cb[,1] == 0,,drop=FALSE]
    int <- cb[cb[,1] != 0 & cb[,3] != 0,,drop=FALSE]
    yy <- b0[i] + X %*% main[,3,drop=FALSE]
    dm <- matrix(nrow = nrow(X), ncol=nrow(int))
    if(nrow(int) > 0) {
      for(j in 1:nrow(int)) {
        tmp <- X[,int[j,1:2]]
        dm[,j] <- tmp[,1] * tmp[,2]
      }
      yy <- yy + dm %*% int[,3,drop=FALSE]
    }
    preds[,i] <- yy
  }
  preds
}

rf.fun <- function(nodesize, mtry) {
  splitrule <- ifelse(y_bin, "gini", "variance")
  if(y_bin) {
    form <- factor(y) ~ .
  } else {
    form <- y ~ .
  }
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "y_bin")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val
      
      # Set seed, as random forests model fitting procedure involves
      # random number generation
      set.seed(123+i)

      colnames(split_data)[1] <- "y"
      model <- ranger::ranger(formula = form, data = split_data, mtry = mtry[i], min.node.size = floor(nodesize[i] * nrow(split_data)),
                              num.trees = 2000, num.threads = 1, probability = y_bin, splitrule = splitrule)
      preds <- predict(model, split_data_val[,-1], num.threads = 1)$predictions
      if(y_bin) preds <- preds[,2]
      # Return iteration and score
      data.frame(Iter = i, Score = eval.performance(preds, split_data_val[,1]))
    }
  return(list(Score = scores))
}

lr.fun <- function(ntrees, nleaves) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  lr_type <- ifelse(y_bin, 3, 2)
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=rbind,
                    .export = c("n_reps", "getData", "scenarios", "scen", "study", "data.complete", "hyperopt", "eval.performance", "y_bin", "E",
                                "splitSNPs", "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split SNPs into dominant and recessive to not lose any information
      split_data_bin <- splitSNPs(split_data[,-c(1, rem_E)])
      split_data_val_bin <- splitSNPs(split_data_val[,-c(1, rem_E)])

      # Set seed, as logic regression model fitting procedure involves
      # random number generation
      set.seed(123+i)

      # Setting manually tuned temperature parameters
      if(study == "linear") {
        temps <- c(-1, -3.5)
      } else if(study == "interactions_no_E") {
        temps <- c(-1.5, -4)
      } else if(study == "interactions_E") {
        temps <- c(-1.5, -4)
      } else if(study == "complex") {
        temps <- c(-1, -3.25)
      } else if(study == "salia") {
        temps <- c(1, -0.5)
      } else if(study == "dorothea") {
        temps <- c(-1, -3.25)
      }

      if(!E) {
        model <- LogicReg::logreg(resp = split_data[,1], bin = split_data_bin, type = lr_type, select = 1,
                                  anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                                  ntrees = ntrees[i], nleaves = nleaves[i],
                                  tree.control = LogicReg::logreg.tree.control(opers=2))
        score <- eval.performance(predict(model, newbin = split_data_val_bin), split_data_val[,1])
      } else {
        model <- LogicReg::logreg(resp = split_data[,1], bin = split_data_bin, sep = split_data[,"E",drop=FALSE],
                                  type = lr_type, select = 1,
                                  anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                                  ntrees = ntrees[i], nleaves = nleaves[i],
                                  tree.control = LogicReg::logreg.tree.control(opers=2))
        score <- eval.performance(predict(model, newbin = split_data_val_bin, newsep = as.matrix(split_data_val[,"E",drop=FALSE])), split_data_val[,1])
      }
      if(is.na(score)) {
        score <- 0
      }
      # Return iteration and score
      data.frame(Iter = i, Score = score)
    }
  return(list(Score = scores))
}

###############################

# Interpreting models by extracting resulting terms

# Custom combine function for foreach to combine lists
combine.custom <- function(a, b) {
  # https://stackoverflow.com/a/13433689
  depth <- function(this,thisdepth=0){
    if(!is.list(this)){
      return(thisdepth)
    }else{
      return(max(unlist(lapply(this,depth,thisdepth=thisdepth+1))))
    }
  }

  if(depth(a) > depth(b)) {
    res <- c(a, list(b))
  } else {
    res <- c(list(a), list(b))
  }
  return(res)
}


BITS.terms <- function(hyperparams) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=combine.custom,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "E")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split off environmental variable (if it exists) and obtain modes of
      # inheritance for all SNPs
      Z <- NULL -> Z_val
      if(E) {
        Z <- split_data[,"E",drop=FALSE]
        Z_val <- split_data_val[,"E",drop=FALSE]
      }
      mois <- BITS::MOI(split_data[,-c(1, rem_E)], split_data[,1])
      split_data_bin <- BITS::applyMOI(split_data[,-c(1, rem_E)], mois)
      split_data_val_bin <- BITS::applyMOI(split_data_val[,-c(1, rem_E)], mois)

      # Potentially concatenate environmental data back
      if(E) {
        split_data_bin <- cbind(split_data_bin, Z); split_data_val_bin <- cbind(split_data_val_bin, Z_val)
        Z <- NULL -> Z_val
      }

      gamma <- hyperparams$gamma[hyperparams$Iter == i][1]
      model <- BITS::BITS(split_data_bin, split_data[,1],
                          max.vars = 3, gamma = gamma,
                          boosting.iter = 50, learning.rate = 0.1,
                          lambda = NULL, alpha = 1, nfolds = 0, nlambda = 1000,
                          term.select = "elnet",
                          relax = TRUE, gamma2 = 0.5,
                          max.iter = function(p) -1)
      lm <- model$lin.mod
      terms <- list()
      n.lambda <- length(lm$lambda)
      n.terms.so.far <- integer()
      # For each number of terms (up to 10), get set of induced terms
      for(j in 1:n.lambda) {
        s <- lm$lambda[j]
        model$s <- s
        tmp.terms <- BITS::get.included.vars(model)
        tmp.terms <- unique(tmp.terms)
        n.terms <- nrow(tmp.terms)
        if(!(n.terms %in% n.terms.so.far) && n.terms > 0 && n.terms <= 10) {
          n.terms.so.far <- c(n.terms.so.far, n.terms)
          terms[[n.terms]] <- tmp.terms
        }
      }
      
      # Optimal terms:
      s <- hyperparams$lambda[hyperparams$Iter == i][1]
      model$s <- s
      optimal.terms <- unique(BITS::get.included.vars(model))
      model$set_vars <- NULL; gc() # Clean up due to parallel computation
      # Return list of iteration, terms for all numbers of terms, and set of optimal terms
      list(Iter = i, terms = terms, optimal.terms = optimal.terms)
    }
  return(scores)
}

glinternet.terms <- function(hyperparams) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=combine.custom,
                    .export = c("n_reps", "getData", "scenarios", "scen", "data.complete", "hyperopt", "eval.performance", "E", "y_bin",
                                "optimal.lambda", "calcScorePerObservation", "calcScore",
                                "splitSNPs", "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split SNPs into dominant and recessive to not lose any information
      split_data_bin <- splitSNPs(split_data[,-c(1, rem_E)])
      split_data_val_bin <- splitSNPs(split_data_val[,-c(1, rem_E)])
      train.dat <- split_data_bin; val.dat <- split_data_val_bin
      n.levels <- rep(2, ncol(split_data_bin))

      Z <- NULL -> Z_val
      if(E) {
        Z <- split_data[,"E",drop=FALSE]
        Z_val <- split_data_val[,"E",drop=FALSE]
        train.dat <- cbind(train.dat, Z); val.dat <- cbind(val.dat, Z_val)
        n.levels <- c(n.levels, 1)
      }
      s <- NULL
      model <- glinternet::glinternet(train.dat, split_data[,1], n.levels, numCores = 1,
                                      family = ifelse(y_bin, "binomial", "gaussian"),
                                      lambda = s, nLambda = 1000)
      # For each number of terms (up to 10), get set of induced terms
      lm <- model$activeSet
      terms <- list()
      n.lambda <- length(model$lambda)
      n.terms.so.far <- integer()
      for(j in 1:n.lambda) {
        tmp.terms <- lm[[j]]
        tmp.terms2 <- matrix(nrow=0,ncol=2)
        n.terms <- 0
        
        # Main effects
        if(!(is.null(tmp.terms$cat))) {
          n.terms <- n.terms + nrow(tmp.terms$cat)
          tmp.terms2 <- rbind(tmp.terms2, cbind(tmp.terms$cat, NA))
        }
        if(!(is.null(tmp.terms$cont))) {
          n.terms <- n.terms + nrow(tmp.terms$cont)
          # Adjust variable indices, as cat and cont variable indices start
          # at 1, respectively
          cont.terms <- tmp.terms$cont + ncol(split_data_bin)
          tmp.terms2 <- rbind(tmp.terms2, cbind(cont.terms, NA))
        }
        
        # Interaction effects
        if(!(is.null(tmp.terms$catcat))) {
          n.terms <- n.terms + nrow(tmp.terms$catcat)
          tmp.terms2 <- rbind(tmp.terms2, tmp.terms$catcat)
        }
        if(!(is.null(tmp.terms$catcont))) {
          n.terms <- n.terms + nrow(tmp.terms$catcont)
          catcont.terms <- tmp.terms$catcont
          # Adjust variable indices, as cat and cont variable indices start
          # at 1, respectively
          catcont.terms[,2] <- catcont.terms[,2] + ncol(split_data_bin)
          tmp.terms2 <- rbind(tmp.terms2, catcont.terms)
        }
        
        if(!(n.terms %in% n.terms.so.far) && n.terms > 0 && n.terms <= 10) {
          n.terms.so.far <- c(n.terms.so.far, n.terms)
          terms[[n.terms]] <- tmp.terms2
        }
      }

      # Optimal terms:
      s <- hyperparams$lambda[hyperparams$Iter == i][1]
      ll <- model$lambda
      ind <- which.min(abs(ll - s))
      lm <- model$activeSet
      tmp.terms <- lm[[ind]]
      tmp.terms2 <- matrix(nrow=0,ncol=2)
      if(!(is.null(tmp.terms$cat))) {
        tmp.terms2 <- rbind(tmp.terms2, cbind(tmp.terms$cat, NA))
      }
      if(!(is.null(tmp.terms$cont))) {
        cont.terms <- tmp.terms$cont + ncol(split_data_bin)
        tmp.terms2 <- rbind(tmp.terms2, cbind(cont.terms, NA))
      }
      if(!(is.null(tmp.terms$catcat))) {
        tmp.terms2 <- rbind(tmp.terms2, tmp.terms$catcat)
      }
      if(!(is.null(tmp.terms$catcont))) {
        catcont.terms <- tmp.terms$catcont
        catcont.terms[,2] <- catcont.terms[,2] + ncol(split_data_bin)
        tmp.terms2 <- rbind(tmp.terms2, catcont.terms)
      }
      optimal.terms <- tmp.terms2

      cat(sprintf("Iteration %d done\n", i), sep = "")
      # Return list of iteration, terms for all numbers of terms, and set of optimal terms
      list(Iter = i, terms = terms, optimal.terms = optimal.terms)
    }
  return(scores)
}

lr.terms <- function(hyperparams) {
  rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
  lr_type <- ifelse(y_bin, 3, 2)

  # Parallel over all repetitions:
  scores <- foreach(i=1:n_reps,
                    .combine=combine.custom,
                    .export = c("n_reps", "getData", "scenarios", "scen", "study", "data.complete", "hyperopt", "eval.performance", "y_bin", "E",
                                "splitSNPs", "snps")) %dopar%
    {
      # Get training and evaluation data sets for this repetition:
      dat <- getData(i, scen, hyperopt)
      split_data <- dat$train; split_data_val <- dat$val

      # Split SNPs into dominant and recessive to not lose any information
      split_data_bin <- splitSNPs(split_data[,-c(1, rem_E)])
      split_data_val_bin <- splitSNPs(split_data_val[,-c(1, rem_E)])

      # Set seed, as logic regression model fitting procedure involves
      # random number generation
      set.seed(123+i)

      # Setting manually tuned temperature parameters
      if(study == "linear") {
        temps <- c(-1, -3.5)
      } else if(study == "interactions_no_E") {
        temps <- c(-1.5, -4)
      } else if(study == "interactions_E") {
        temps <- c(-1.5, -4)
      } else if(study == "complex") {
        temps <- c(-1, -3.25)
      } else if(study == "salia") {
        temps <- c(1, -0.5)
      } else if(study == "dorothea") {
        temps <- c(-1, -3.25)
      }

      terms <- list()
      optimal.ntrees <- hyperparams$ntrees[hyperparams$Iter == i][1]
      # For each number of terms (up to maximum number of considered logic trees), get set of induced terms
      for(ntrees in 1:max(hyperparams$ntrees)) {
        nleaves <- hyperparams$nleaves[hyperparams$Iter == i & hyperparams$ntrees == ntrees][1]
        if(!E) {
          model <- LogicReg::logreg(resp = split_data[,1], bin = split_data_bin, type = lr_type, select = 1,
                                    anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                                    ntrees = ntrees, nleaves = nleaves,
                                    tree.control = LogicReg::logreg.tree.control(opers=2))
        } else {
          model <- LogicReg::logreg(resp = split_data[,1], bin = split_data_bin, sep = split_data[,"E",drop=FALSE],
                                    type = lr_type, select = 1,
                                    anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                                    ntrees = ntrees, nleaves = nleaves,
                                    tree.control = LogicReg::logreg.tree.control(opers=2))
        }
        trees <- model$model$trees
        tmp.terms <- matrix(nrow=0, ncol=2)
        for(tree in trees) {
          vars <- tree$trees$knot
          vars <- vars[vars > 0]
          if(length(vars) > ncol(tmp.terms)) {
            tmp.terms <- cbind(tmp.terms, matrix(NA, nrow=nrow(tmp.terms), ncol=length(vars)-ncol(tmp.terms)))
            tmp.terms <- rbind(tmp.terms, matrix(vars, nrow=1))
          } else {
            tmp.terms <- rbind(tmp.terms, matrix(c(vars, rep(NA, ncol(tmp.terms)-length(vars))), nrow=1))
          }
        }
        terms[[length(trees)]] <- tmp.terms

        # Optimal terms:
        if(ntrees == optimal.ntrees) {
          optimal.terms <- tmp.terms
        }
      }
      # Return list of iteration, terms for all numbers of terms, and set of optimal terms
      list(Iter = i, terms = terms, optimal.terms = optimal.terms)
    }
  return(scores)
}

###############################

# Grid search for hyperparameter optimization

grid.search <- function(FUN, bounds) {
  n_params <- ncol(bounds)
  iters <- nrow(bounds)
  
  eval <- NULL

  bounds2 <- as.list(bounds)
  
  for(i in 1:iters) {
    cat(sprintf("Iteration %d/%d (%.0f%%) started", i, iters, i/iters * 100), sep = "")

    current_params <- as.list(bounds[i,,drop=FALSE])

    cat("\r", sprintf("Iteration %d/%d (%.0f%%) |", i, iters, i/iters * 100), sep = "")
    lapply(seq_along(current_params), function(x) cat(sprintf(" %s:", names(current_params)[x]), current_params[[x]]))
    flush.console()
    cat("\n")

    # Call model fitting and evaluation function with currently considered set
    # of hyperparameters over all repetitions
    params_rep <- lapply(current_params, function(x) rep(x, n_reps))
    score <- do.call(FUN, params_rep)
    eval <- rbind(eval, cbind(score$Score, as.data.frame(params_rep)))
  }

  return(list(scoreSummary = eval, bounds = bounds2))
}

###############################

# Explicitly considered hyperparameter settings
# (for hyperparameters that are not implicitly determined such as penalty parameters)

BITS.bounds <- expand.grid(dummy = 0)

glinternet.bounds <- expand.grid(dummy = 0)

rulefit.bounds <- expand.grid(subsample.frac = c(0.5, 0.75, 1),
                              nodesize = c(0.01, 0.05, 0.1))

mars.bounds <- expand.grid(dummy = 0)

sprinter.bounds <- expand.grid(dummy = 0)

mtries <- c(floor(0.5 * floor(sqrt(n.snps+E-1))),
            floor(1 * floor(sqrt(n.snps+E-1))),
            floor(2 * floor(sqrt(n.snps+E-1))))
rf.bounds <- expand.grid(nodesize = c(0.01, 0.05, 0.1),
                         mtry = mtries)

trees_leaves_combs <- expand.grid(nleaves = 1:10, ntrees = 1:5)
trees_leaves_combs <- expand.grid(nleaves = 1:12, ntrees = 1:6)
trees_leaves_combs <- trees_leaves_combs[trees_leaves_combs[,1] >= trees_leaves_combs[,2],]
lr.bounds <- trees_leaves_combs

###############################

# Methods for tuning hyperparameters, evaluating methods, and extracting terms

tune.hyperparameters <- function(method, algo = "") {
  algo <<- algo
  hyperopt <<- TRUE
  method.file.string <- paste(method, ifelse(algo == "", "", "."), algo, sep="")
  for(s in scenes) {
    data.complete <<- data.list[[s]]
    scen <<- s
    result <- grid.search(FUN = get(paste(method, ".fun", sep="")),
                          bounds = get(paste(method, ".bounds", sep="")))
    load(params_file)
    results[[method.file.string]][[s]] <- result
    save(results, file = params_file)
  }
}

eval.method <- function(method, algo = "") {
  algo <<- algo
  hyperopt <<- FALSE
  method.file.string <- paste(method, ifelse(algo == "", "", "."), algo, sep="")
  for(s in scenes) {
    data.complete <<- data.list[[s]]
    scen <<- s
    load(params_file)
    # Get optimal hyperparameters
    hyper <- results[[method.file.string]][[s]]
    hyperparams <- hyper$scoreSummary
    new_params <- NULL
    for(i in 1:n_reps) {
      tmp <- hyperparams[hyperparams$Iter == i,]
      new_params <- rbind(new_params, tmp[which.max(tmp$Score),,drop=FALSE])
    }
    hyperparams <- new_params
    hyperparams <- hyperparams[order(hyperparams$Iter), 3:ncol(hyperparams), drop=FALSE]
    hyperparams <- as.list(hyperparams)
    # Evaluate method using the best hyperparameters
    score <- do.call(get(paste(method, ".fun", sep="")), hyperparams)
    load(scores_file)
    scores[[method.file.string]][[s]] <- score$Score
    save(scores, file = scores_file)
  }
}

term.extract <- function(method, algo = "") {
  algo <<- algo
  hyperopt <<- FALSE
  method.file.string <- paste(method, ifelse(algo == "", "", "."), algo, sep="")
  for(s in scenes) {
    data.complete <<- data.list[[s]]
    scen <<- s
    load(params_file)
    hyper <- results[[method.file.string]][[s]]
    hyperparams <- hyper$scoreSummary
    hyperparams <- hyperparams[order(hyperparams$Score, decreasing = TRUE),,drop=FALSE]
    # Extract fitted terms using the best hyperparameters
    score <- do.call(get(paste(method, ".terms", sep="")), list(hyperparams=hyperparams))
    load(terms_file)
    terms[[method.file.string]][[s]] <- score
    save(terms, file = terms_file)
  }
}
