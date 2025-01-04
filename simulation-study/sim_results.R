# Continue here locally processing the obtained results from sim_run.R
# by constructing plots

# Set number of cores
use_cores <- 1
# Force this number of cores, otherwise use all available cores
force_cores <- TRUE

plot.dat <- NULL; scenes <- 1:9
studies <- c("linear", "interactions_no_E", "interactions_E", "complex")
# For every study, assemble the complete.eval data gather it in plot.dat:
for(study in studies) {
  source("../common/prep.R")
  
  # Get true model performance for comparison
  true.scores <- list()
  for(s in scenes) {
    data.complete <- data.list[[s]]
    lin.pred <- lin.pred.list[[s]]
    true.scores[[s]] <- rep(0.0, n_reps)
    for(i in 1:n_reps) {
      split_data_val <- getData(i, s, FALSE)$val
      true.scores[[s]][i] <- eval.performance(lin.pred[getData(i, s, FALSE)$test_inds], split_data_val[,1])
    }
  }
  
  ##########################################
  
  # Collect plot data
  
  load(scores_file)
  
  metric <- ifelse(!y_bin, "1-NMSE", "Area under the Curve")
  
  complete.eval <- NULL
  plot.scenes <- scenes
  for(s in plot.scenes) {
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "BITS", Performance = scores$BITS[[s]]$Score, s = s))
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "glinternet", Performance = scores$glinternet[[s]]$Score, s = s))
    if(!y_bin) complete.eval <- rbind(complete.eval, data.frame(Algorithm = "sprinter", Performance = scores$sprinter[[s]]$Score, s = s))
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "MARS", Performance = scores$mars[[s]]$Score, s = s))
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "RuleFit", Performance = scores$rulefit[[s]]$Score, s = s))
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "Random Forests", Performance = scores$rf[[s]]$Score, s = s))
    complete.eval <- rbind(complete.eval, data.frame(Algorithm = "Logic Regression", Performance = scores$lr[[s]]$Score, s = s))
    if(exists("true.scores")) {
      complete.eval <- rbind(complete.eval, data.frame(Algorithm = "True Model", Performance = true.scores[[s]], s = s))
    }
  }
  complete.eval$Algorithm <- factor(complete.eval$Algorithm,
                                    levels = unique(complete.eval$Algorithm),
                                    ordered = TRUE)
  
  ######################################
  
  plot.dat <- rbind(plot.dat, cbind(complete.eval, scenarios[complete.eval$s,], Study = study))
  
}

study.string <- function(study) {
  study[study == "linear"] <- "Linear"
  study[study == "interactions_no_E"] <- "Interactions I"
  study[study == "interactions_E"] <- "Interactions II"
  study[study == "complex"] <- "Complex"
  study
}

plot.dat$snr <- factor(plot.dat$snr,
                       levels = unique(plot.dat$snr),
                       ordered = TRUE)
plot.dat$N <- paste("N =", plot.dat$N)
plot.dat$N <- factor(plot.dat$N,
                     levels = unique(plot.dat$N),
                     ordered = TRUE)
plot.dat$Study <- study.string(plot.dat$Study)
plot.dat$Study <- factor(plot.dat$Study, levels = unique(plot.dat$Study))

library(ggplot2)

# Full plot depicting all scenarios and settings
p <- ggplot(data = plot.dat, aes(x=.data[["snr"]], y=.data[["Performance"]], fill=Algorithm)) +
  stat_boxplot(geom = "errorbar", lwd=.1) +
  geom_boxplot(lwd=.1, outlier.size=.5) +
  scale_x_discrete("SNR [Signal-to-Noise Ratio]", waiver()) +
  scale_y_continuous(latex2exp::TeX("Predictive Performance [$R^2$]")) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(8, "Set1")[c(1:6, 8, 7)]) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(rows = vars(Study), cols = vars(N))
ggsave(plot = p, filename = "results.pdf", width = 8, height = 12)

# Medium settings for main text
plot.dat.medium <- plot.dat[plot.dat$snr == 1 & plot.dat$N == "N = 1000",]
p <- ggplot(plot.dat.medium, aes(x=Algorithm, y=Performance, fill=Algorithm)) +
  stat_boxplot(geom = "errorbar", lwd=.1) +
  geom_boxplot(lwd=.1, outlier.size=.5) +
  scale_x_discrete("Algorithm", waiver()) +
  scale_y_continuous(latex2exp::TeX("Predictive Performance [$R^2$]")) +
  scale_fill_manual(values=RColorBrewer::brewer.pal(8, "Set1")[c(1:6, 8, 7)]) +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  facet_wrap(~ Study, nrow=1, scales="free")
ggsave(plot = p, filename = "results_medium.pdf", width = 8, height = 3)

##########################################

# Hyperparameter plot

# Using complex study, medium setting (5), and first replication

i <- 1; scen <- 5; hyperopt <- TRUE
data.complete <- data.list[[scen]]

dat <- getData(i, scen, hyperopt)
split_data <- dat$train; split_data_val <- dat$val

rem_E <- ifelse(E, ncol(data.complete), ncol(data.complete)+1)
Z <- NULL -> Z_val
fast <- ifelse(ncol(split_data)-1 <= 1000, FALSE, TRUE)
mois <- BITS::MOI(split_data[,-c(1, rem_E)], split_data[,1], fast = fast)
split_data_bin <- BITS::applyMOI(split_data[,-c(1, rem_E)], mois)
split_data_val_bin <- BITS::applyMOI(split_data_val[,-c(1, rem_E)], mois)

if(E) {
  split_data_bin <- cbind(split_data_bin, Z); split_data_val_bin <- cbind(split_data_val_bin, Z_val)
  Z <- NULL -> Z_val
}

model <- BITS::BITS.complete(split_data_bin, split_data[,1],
                             max.vars = 3,
                             boosting.iter = 50, learning.rate = 0.1,
                             lambda = NULL, alpha = 1, nfolds = 0,
                             term.select = "elnet",
                             relax = TRUE, gamma2 = 0.5,
                             max.iter = function(p) -1, parallel = FALSE,
                             gmin = function(gmax) 0.001 * gmax, gsteps = 50)

result <- val.res$val.res; class(result ) <- c("val.performance", "data.frame")
p <- plot(result) + ylab("Score [Mean Squared Error]")
p
ggsave(p, filename = "val.res.plot.pdf", width = 6, height = 3)
