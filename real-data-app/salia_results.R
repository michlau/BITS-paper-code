# Process SALIA results here locally

# Set number of cores
use_cores <- 1
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE

# Load SALIA data (as in salia_run.R):

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

####################################

# Create predictive performance plot:

load(scores_file)

metric <- "Area under the Curve"

complete.eval <- NULL
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "BITS", Performance = scores$BITS[[1]]$Score))
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "glinternet", Performance = scores$glinternet[[1]]$Score))
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "MARS", Performance = scores$mars[[1]]$Score))
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "RuleFit", Performance = scores$rulefit[[1]]$Score))
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "Random Forests", Performance = scores$rf[[1]]$Score))
complete.eval <- rbind(complete.eval, data.frame(Algorithm = "Logic Regression", Performance = scores$lr[[1]]$Score))
complete.eval$Algorithm <- factor(complete.eval$Algorithm,
                                  levels = unique(complete.eval$Algorithm),
                                  ordered = TRUE)

p <- ggplot(complete.eval, aes(x=Algorithm, y=Performance, fill=Algorithm)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_x_discrete("Algorithm", waiver()) +
  scale_y_continuous(metric) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(8, "Set1")[c(1,2,4,5,6,8)]) +
  theme_bw() +
  theme(legend.position="bottom")
ggsave(plot = p, filename = "salia_results.pdf", width = 9, height = 5.73)

##########################################

# Model interpretation:

load(params_file)

# BITS
mois <- BITS::MOI(data.complete[,-1], data.complete[,1], fast = FALSE)
split_data_bin <- BITS::applyMOI(data.complete[,-1], mois)
snp.names <- colnames(data.complete[,-1])

params <- results$BITS[[1]]$scoreSummary
lambda <- median(params$lambda)
gamma <- median(params$gamma)

model <- BITS::BITS(split_data_bin, data.complete[,1],
                    max.vars = 3,
                    gamma = gamma,
                    boosting.iter = 50, learning.rate = 0.1,
                    lambda = NULL, alpha = 1, nfolds = 0,
                    term.select = "elnet",
                    relax = TRUE, gamma2 = 0.5,
                    max.iter = function(x) -1)
model$s <- lambda
vars <- BITS::get.included.vars(model)
vars
nrow(vars) # Number of terms
sum(is.na(vars[,2])) # Number of main effects
sum(!is.na(vars[,2]) & is.na(vars[,3])) # Number of second order interaction effects
sum(!is.na(vars[,3])) # Number of third order interaction effects
print(model, X.names=snp.names, MOI=mois, SNP=TRUE, latex=TRUE)

# glinternet interpretation
params <- results$glinternet[[1]]$scoreSummary
lambda <- median(params$lambda)
train.dat <- splitSNPs(data.complete[,-1])
n.levels <- rep(2, ncol(train.dat))
model <- glinternet::glinternet(train.dat, data.complete[,1], n.levels, numCores = 4,
                                family = "binomial")
ll <- model$lambda
ind <- which.min(abs(ll - lambda))
s <- ll[ind]
lm <- model$activeSet
tmp.terms <- lm[[ind]]
tmp.terms2 <- matrix(nrow=0,ncol=2)
if(!(is.null(tmp.terms$cat))) {
  tmp.terms2 <- rbind(tmp.terms2, cbind(tmp.terms$cat, NA))
}
if(!(is.null(tmp.terms$catcat))) {
  tmp.terms2 <- rbind(tmp.terms2, tmp.terms$catcat)
}
optimal.terms <- tmp.terms2
optimal.terms

# Logic regression interpretation
hyperparams <- results$lr[[1]]$scoreSummary
new_params <- NULL
for(i in 1:100) {
  tmp <- hyperparams[hyperparams$Iter == i,]
  new_params <- rbind(new_params, tmp[which.max(tmp$Score),,drop=FALSE])
}
params <- as.list(apply(new_params, 2, median))
temps <- c(1, -0.5)
split_data_bin <- splitSNPs(data.complete[,-1])
set.seed(123)
model <- LogicReg::logreg(resp = data.complete[,1], bin = split_data_bin, type = 3, select = 1,
                          anneal.control = LogicReg::logreg.anneal.control(temps[1], temps[2], 500000, update = -1),
                          ntrees = params$ntrees, nleaves = params$nleaves,
                          tree.control = LogicReg::logreg.tree.control(opers=2))
model
plot(model)
