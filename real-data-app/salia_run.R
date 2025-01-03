# Run hyperparameter optimization and method evaluation for the SALIA application

# Set number of cores
use_cores <- 32
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE

study <- "salia"

scenarios <- 1
scenes <- 1

load("salia_data.RData")
resp.ind <- 3
snp.offset <- match(TRUE, startsWith(colnames(data), "rs"))
data <- data[!is.na(data[,resp.ind]),]
y_bin <- setequal(unique(data[,resp.ind]), c(0, 1))
if(y_bin) data[,resp.ind] <- as.integer(data[,resp.ind])
n.snps <- sum(startsWith(colnames(data), "rs"))

data.complete <- data[,c(resp.ind, snp.offset:(snp.offset + n.snps - 1))]
if(exists("E.offset")) {
  data.complete$E <- data[,E.offset]
}

E <- "E" %in% colnames(data.complete)
y_bin <- setequal(unique(data.complete[,1]), c(0, 1))

source("../common/prep.R")

n_reps <- 100

data.list <- list()
data.list[[1]] <- data.complete

getData <- function(i, scen, hyperopt) {
  set.seed(123)
  train_inds <- list()
  test_inds <- list()
  val_inds <- list()
  N <- nrow(data.complete)
  for(j in 1:n_reps) {
    train_inds[[j]] <- sample(1:N, floor(0.7 * N))
    test_inds[[j]] <- -train_inds[[j]]
    val_inds[[j]] <- sample(1:length(train_inds[[j]]), floor(0.25 * length(train_inds[[j]])))
  }
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

algo <- ""
method <- c("BITS", "glinternet", "rulefit", "mars", "rf", "lr")

for(m in method) {
  tune.hyperparameters(m, algo = algo)
  eval.method(m, algo = algo)
}

stopCluster(cl)

# Quit script here on the computer cluster
quit("no")
