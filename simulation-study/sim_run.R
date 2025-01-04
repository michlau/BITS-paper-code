# Run hyperparameter optimization, method evaluation and term extraction
# Preferably on a compute cluster with parallel jobs per simulation setting
# (Provide simulation study index (1,...,4) and simulation study setting (1,...,9)
#  as command line arguments)

# Set number of cores
use_cores <- 32
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE

studies <- c("linear", "interactions_no_E", "interactions_E", "complex")

args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  study <- studies[as.numeric(args[1])]
  scenes <- as.numeric(args[2])
}
cat(study, scenes); cat("\n")

algo <- ""
method <- c("BITS", "glinternet", "sprinter", "rulefit", "mars", "rf", "lr")

source("../common/prep.R")

for(m in method) {
  cat(m); cat("\n")

  tune.hyperparameters(m, algo = algo)
  eval.method(m, algo = algo)
  term.extract(m, algo = algo)
}

stopCluster(cl)

# Quit script here on the computer cluster
quit("no")

###########################################
