# Measure time of model fitting and evaluation for each considered method
# Execute this script for each method (in parallel) and supply corresponding
# method index (1,...,8) as command-line argument

# Only use 1 core to obtain comparable model fitting times
use_cores <- 1
# Force this number of cores
force_cores <- TRUE

# Only for complex scenario and medium (5) setting
study <- "complex"
scenes <- 5
times_file <- "times.RData"

algo <- ""

# BITS twice: once using complete search, once using hybrid search
methods <- c("BITS", "BITS", "glinternet", "rulefit", "mars", "sprinter", "lr", "rf")
args <- commandArgs(trailingOnly = TRUE)
if(length(args) > 0) {
  method <- methods[as.numeric(args[1])]
}

source("../common/prep.R")

for(m in method) {
  m.name <- m
  if(as.numeric(args[1]) > 1 && m == "BITS") {
    complete.search <<- FALSE
    m.name <- "BITS.hybrid"
  }

  cat(sprintf("%s, %s", m, m.name)); cat("\n")

  tic <- Sys.time()
  eval.method(m, algo = algo)
  toc <- Sys.time()
  eval.time <- as.numeric(toc-tic, units="secs")/100 # 100 reps
  load(times_file)
  eval.times[[m.name]] <- eval.time
  save(tune.times, eval.times, file = times_file)
}

stopCluster(cl)
quit("no")

###########################

# Show results

load(times_file)
eval.times
lapply(eval.times, function(x) round(x, 1))


