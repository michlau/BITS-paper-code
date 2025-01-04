# Run hyperparameter optimization and method evaluation for the Dorothea application

# Set number of cores
use_cores <- 32
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE

study <- "dorothea"

scenarios <- 1
scenes <- 1

E <- FALSE

source("../common/prep.R")

n_reps <- 1

loadData <- function(full = TRUE) {
  if(full) load("dorothea_neurips.RData") else load("dorothea_neurips_lasso.RData")
  data <- dat
  resp.ind <<- 1
  snp.offset <<- 2
  data[,resp.ind] <- as.integer(data[,resp.ind])
  n.snps <<- ncol(data) - 1

  data.complete <<- data[,c(resp.ind, snp.offset:(snp.offset + n.snps - 1))]
  if(exists("E.offset")) {
    data.complete$E <<- data[,E.offset]
  }

  E <<- "E" %in% colnames(data.complete)
  y_bin <<- TRUE

  snps <<- FALSE
  complete.search <<- !full

  data.list <<- list()
  data.list[[1]] <<- data.complete
}

getData <- function(i, scen, hyperopt) {
  set.seed(123)
  train_inds <- list(1:1150)
  test_inds <- list(1151:1950)
  val_inds <- list(801:1150)
  N <- nrow(data.complete)
  set.seed(NULL)

  split_data <- data.complete[train_inds[[i]],]
  if(hyperopt) {
    split_data_val <- split_data[val_inds[[i]],]
    split_data <- split_data[-val_inds[[i]],]
  } else {
    split_data_val <- data.complete[test_inds[[i]],]
  }
  return(list(train = split_data, val = split_data_val,
              train_inds = train_inds[[i]], test_inds = test_inds[[i]], val_inds = val_inds[[i]]))
}

set.seed(NULL)

##########################################

# Custom method evaluation functions,
# as Dorothea is a high-dimensional data set (potentially needing different handling)
# and the raw predictions for all test data observations are needed

BITS.fun <- function(dummy, gamma = rep(-1, n_reps), lambda = rep(-1, n_reps)) {
  i <- 1
  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- dat$train; split_data_val <- dat$val

  split_data_bin <- split_data[,-1]
  split_data_val_bin <- split_data_val[,-1]

  # Hybrid search with 4 times as many iterations as a standard greedy search
  # (due to max.vars = 3)
  max.iter.factor <- 3*4
  if(complete.search) {
    max.iter <- function(p) -1
  } else {
    max.iter <- function(p) max.iter.factor*p
  }

  if(hyperopt) s <- NULL else s <- lambda[i]
  if(hyperopt) {
    model <- BITS::BITS.complete(split_data_bin, split_data[,1],
                                 max.vars = 3,
                                 boosting.iter = 100, learning.rate = 0.1,
                                 lambda = NULL, alpha = 1, nfolds = 0,
                                 term.select = "elnet",
                                 relax = TRUE, gamma2 = 0.5,
                                 max.iter = max.iter, parallel = TRUE,
                                 gmin = function(gmax) 0.001 * gmax, gsteps = 50)
    # Get optimal gamma and lambda penalties
    best.mod <- BITS::get.ideal.model(model, split_data_val_bin, split_data_val[,1], choose = "min", metric = "auc")
    g <- best.mod$best.g; s <- best.mod$best.s
    score <- min(best.mod$val.res$score)
  } else {
    model <- BITS::BITS(split_data_bin, split_data[,1],
                        max.vars = 3, gamma = gamma[i],
                        boosting.iter = 100, learning.rate = 0.1,
                        lambda = NULL, alpha = 1, nfolds = 0,
                        term.select = "elnet",
                        relax = TRUE, gamma2 = 0.5,
                        max.iter = max.iter)
    model$s <- s # Set optimal term inclusion penalty obtained in hyperparameter optimization
    # Get test data predictions
    preds <- predict(model, split_data_val_bin)
    score <- eval.performance(preds, split_data_val[,1])
    g <- gamma[i]
  }
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score, lambda = s, gamma = g)))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score, lambda = s, gamma = g),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

glinternet.fun <- function(dummy, lambda = rep(-1, n_reps)) {
  i <- 1
  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- as.matrix(dat$train); split_data_val <- as.matrix(dat$val)

  split_data_bin <- split_data[,-1]
  split_data_val_bin <- split_data_val[,-1]
  n.levels <- rep(2, ncol(split_data_bin))
  Z <- NULL -> Z_val
  s <- NULL
  if(lambda[i] != -1) s <- lambda[i]
  # For Dorothea, set screenLimit to 1000. Otherwise, memory issues might occur
  # due to high dimensionality.
  screenLimit <- 1000
  model <- glinternet::glinternet(split_data_bin, split_data[,1], n.levels, numCores = use_cores,
                                  family = "binomial",
                                  lambda = NULL,
                                  screenLimit = screenLimit)
  if(lambda[i] == -1)
    s <- optimal.lambda(model, split_data_val_bin, split_data_val[,1], choose = "min")$best.s
  else {
    ll <- model$lambda
    ind <- which.min(abs(ll - s))
    s <- ll[ind]
  }
  cat(sprintf("Iteration %d done\n", i), sep = "")
  # Get test data predictions
  preds <- predict(model, split_data_val_bin, type = "response", lambda = s)
  score <- eval.performance(preds, as.integer(split_data_val[,1]))
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score, lambda = s)))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score, lambda = s),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

rulefit.fun <- function(subsample.frac, nodesize, lambda = rep(-1, n_reps)) {
  fam <- "binomial"
  i <- 1
  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- dat$train; split_data_val <- dat$val

  # Set minbucket/minsplit parameters
  # Also, set maxdepth=3 to allow detecting interactions with order up to 3,
  # as for BITS and MARS
  minbucket <- floor(nodesize[i] * nrow(split_data))
  tree.control <- partykit::ctree_control(maxdepth = 3, mtry = Inf,
                                          minbucket = minbucket, minsplit = 2*minbucket,
                                          testtype = "Univariate", alpha = 0.5)
  colnames(split_data)[1] <- "y"
  split_data$y <- as.factor(split_data$y)
  model <- pre::pre(y ~ ., data = split_data, family = fam,
                    ntrees = 50, tree.control = tree.control,
                    sampfrac = subsample.frac[i], learnrate = 0.1,
                    verbose = TRUE)

  if(lambda[i] == -1)
    s <- optimal.lambda(model, split_data_val[,-1], split_data_val[,1], choose = "min")$best.s
  else
    s <- lambda[i]

  # Get test data predictions
  preds <- predict(model, split_data_val[,-1], type = "response", penalty.par.val = s)
  score <- eval.performance(preds, split_data_val[,1])
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score, lambda = as.numeric(s))))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score, lambda = as.numeric(s)),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

mars.fun <- function(dummy) {
  i <- 1
  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- dat$train; split_data_val <- dat$val

  fam2 <- list(family=binomial())
  model <- earth::earth(split_data[,-1], split_data[,1],
                        degree = 3, glm = fam2)
  # Get test data predictions
  preds <- predict(model, split_data_val[,-1], type = "response")
  score <- eval.performance(preds, split_data_val[,1])
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score)))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

rf.fun <- function(nodesize, mtry) {
  i <- 1
  splitrule <- "gini"
  form <- factor(y) ~ .

  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- dat$train; split_data_val <- dat$val

  colnames(split_data)[1] <- "y"
  model <- ranger::ranger(formula = form, data = split_data, mtry = mtry[i], min.node.size = floor(nodesize[i] * nrow(split_data)),
                          num.trees = 2000, probability = TRUE, splitrule = splitrule)
  # Get test data predictions
  preds <- predict(model, split_data_val[,-1])$predictions
  preds <- preds[,2]
  score <- eval.performance(preds, split_data_val[,1])
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score)))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

lr.fun <- function(ntrees, nleaves) {
  i <- 1
  lr_type <- 3

  # Get training and evaluation data sets:
  dat <- getData(i, scen, hyperopt)
  split_data <- dat$train; split_data_val <- dat$val

  split_data_bin <- split_data[,-1]
  split_data_val_bin <- split_data_val[,-1]

  set.seed(123+i)

  temps <- c(2, 0)

  model <- LogicReg::logreg(resp = split_data[,1], bin = split_data_bin, type = lr_type, select = 1,
                            anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                            ntrees = ntrees[i], nleaves = nleaves[i],
                            tree.control = LogicReg::logreg.tree.control(opers=2))
  # Get test data predictions
  preds <- predict(model, newbin = split_data_val_bin)
  score <- eval.performance(preds, split_data_val[,1])
  if(is.na(score)) {
    score <- 0
  }
  if(hyperopt) {
    return(list(Score = data.frame(Iter = i, Score = score)))
  } else {
    # Also return the full model with test data predictions for the final evaluation
    return(list(Score = list(dat = data.frame(Iter = i, Score = score),
                             model = model, preds = preds, y_val = split_data_val[,1])))
  }
}

##########################################

algo <- ""
# algo <- "lasso" # For ablation analysis of BITS and glinternet using the lasso selection as well
method <- c("BITS", "glinternet", "rulefit", "mars", "rf", "lr")

for(m in method) {
  full <- (m %in% c("BITS", "glinternet")) && (algo == "")
  loadData(full)

  tune.hyperparameters(m, algo = algo)
  eval.method(m, algo = algo)
}

##########################################

stopCluster(cl)
quit("no")

##########################################

