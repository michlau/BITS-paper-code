# Local term identification analysis

# Function to judge how many terms have been correctly identified
howManyCorrect <- function(identified.vars, correct.vars) {
  if(is.null(identified.vars) || nrow(identified.vars) == 0 || is.null(correct.vars) || nrow(correct.vars) == 0) return(0)

  identified.vars <- abs(identified.vars)

  # Ignore environmental covariable in judgment to not put glinternet and
  # logic regression at a disadvantage, as these methods either cannot identify
  # three-way interactions (glinternet) or cannot identify interactions with
  # continuous predictors (logic regression)
  identified.vars[identified.vars == 51] <- NA
  correct.vars[correct.vars == 51] <- NA

  if(ncol(identified.vars) > 1) identified.vars <- unique(t(apply(identified.vars, 1, sort.int, method = "quick", na.last = TRUE)))
  if(ncol(identified.vars) > ncol(correct.vars)) {
    correct.vars <- cbind(correct.vars, matrix(NA, nrow=nrow(correct.vars),
                                               ncol=ncol(identified.vars)-ncol(correct.vars)))
  } else if(ncol(identified.vars) < ncol(correct.vars)) {
    identified.vars <- cbind(identified.vars, matrix(NA, nrow=nrow(identified.vars),
                           ncol=ncol(correct.vars)-ncol(identified.vars)))
  }
  main.disj <- rbind(identified.vars, correct.vars)
  # Count duplicates in concatenated matrix of terms, i.e., term matches
  nrow(main.disj[duplicated(main.disj),,drop=FALSE])
}

# Compute false discovery rates for each method, simulation scenario, setting,
# and number of identified terms
errors <- NULL
for(method in c("BITS", "glinternet", "lr")) {
  for(study in c("linear", "interactions_no_E", "interactions_E", "complex")) {
    terms_file <- paste(study, "_terms.RData", sep="")
    load(terms_file)

    # Data-generating models:
    if(study == "linear")
      correct.vars <- matrix(1:6, ncol=1,byrow=TRUE)
    if(study == "interactions_no_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,NA,
                               1,2,3), ncol=3,byrow=TRUE)
    if(study == "interactions_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,51),ncol=3,byrow=TRUE)
    if(study == "complex")
      correct.vars <- matrix(c(1,NA,NA,
                               3,NA,NA,
                               4,NA,NA,
                               7,8,9,
                               2,51,NA,
                               5,6,51),ncol=3,byrow=TRUE)

    for(s in 1:9) {
      tmp <- terms[[method]][[s]]
      for(k in 1:nrow(correct.vars)) {
        for(i in 1:100) {
          if(length(tmp[[i]]$terms) < k) next
          identified.vars <- tmp[[i]]$terms[[k]]
          if(is.null(identified.vars)) next
          # glinternet and logic regression use both binary SNP variables
          # (SNPD and SNPR) as input
          # Map them to the same SNP index to not put these methods at a disadvantage
          if(method %in% c("glinternet", "lr")) identified.vars <- ceiling(identified.vars/2)
          n.correct.terms <- howManyCorrect(identified.vars, correct.vars)
          errors <- rbind(errors, data.frame(method = method, study = study, s = s, k = k, i = i,
                                             fdr = 1-n.correct.terms/k))
        }
      }
    }
  }
}

se <- function(x) sd(x) / sqrt(length(x))
study.string <- function(study) {
  study[study == "linear"] <- "Linear"
  study[study == "interactions_no_E"] <- "Interactions I"
  study[study == "interactions_E"] <- "Interactions II"
  study[study == "complex"] <- "Complex"
  study
}

errors2 <- aggregate(fdr ~ method + study + s + k, data=errors, FUN = function(x) {
  cbind(mean(x), se(x))
})
errors2$se <- errors2$fdr[,2]; errors2$fdr <- errors2$fdr[,1]
errors2$study <- study.string(errors2$study)
errors2$study <- factor(errors2$study, levels = c("Linear", "Interactions I", "Interactions II", "Complex"))
errors2$method[errors2$method == "lr"] <- "Logic Regression"
colnames(errors2)[1] <- "Algorithm"

scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
errors2$setting <- paste("N = ", scenarios[errors2$s, "N"], "\nSNR = ", scenarios[errors2$s, "snr"], sep="")
errors2$setting <- factor(errors2$setting,
                          levels = unique(errors2$setting),
                          ordered = TRUE)

# Medium setting
errors5 <- errors2[errors2$s == 5,]

library(ggplot2)
p <- ggplot(data = errors5, aes(x=.data[["k"]], y=.data[["fdr"]], fill=Algorithm, color=Algorithm)) +
  geom_line() +
  # Asymptotic 95% confidence interval:
  geom_ribbon(aes(ymin = pmax(0,fdr-1.96*se), ymax = pmin(1,fdr+1.96*se)), alpha = 0.2, color=NA) +
  scale_x_continuous("Number of Identified Terms", breaks = 1:6) +
  scale_y_continuous(latex2exp::TeX("False Discovery Rate")) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(study), scales = "free")
p
ggsave(plot = p, filename = "fdr5.pdf", width = 8, height = 3)

# All settings:
p <- ggplot(data = errors2, aes(x=.data[["k"]], y=.data[["fdr"]], fill=Algorithm, color=Algorithm)) +
  geom_line() +
  geom_ribbon(aes(ymin = pmax(0,fdr-1.96*se), ymax = pmin(1,fdr+1.96*se)), alpha = 0.2, color=NA) +
  scale_x_continuous("Number of Identified Terms", breaks = 1:6) +
  scale_y_continuous(latex2exp::TeX("False Discovery Rate")) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position="bottom")
p + facet_grid(cols = vars(study), rows = vars(setting), scales = "free")
ggsave(plot = p + facet_grid(cols = vars(study), rows = vars(setting), scales = "free"), filename = "fdr.pdf", width = 8, height = 12)

#########################################

# False discovery and false negative rates of optimal models

errors <- NULL
for(method in c("BITS", "glinternet", "lr")) {
  for(study in c("linear", "interactions_no_E", "interactions_E", "complex")) {
    terms_file <- paste(study, "_terms.RData", sep="")
    load(terms_file)

    if(study == "linear")
      correct.vars <- matrix(1:6, ncol=1,byrow=TRUE)
    if(study == "interactions_no_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,NA,
                               1,2,3), ncol=3,byrow=TRUE)
    if(study == "interactions_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,51),ncol=3,byrow=TRUE)
    if(study == "complex")
      correct.vars <- matrix(c(1,NA,NA,
                               3,NA,NA,
                               4,NA,NA,
                               7,8,9,
                               2,51,NA,
                               5,6,51),ncol=3,byrow=TRUE)

    for(s in 1:9) {
      tmp <- terms[[method]][[s]]
      for(i in 1:100) {
        identified.vars <- tmp[[i]]$optimal.terms
        if(is.null(identified.vars)) next
        if(method %in% c("glinternet", "lr")) identified.vars <- ceiling(identified.vars/2)
        n.correct.terms <- howManyCorrect(identified.vars, correct.vars)
        errors <- rbind(errors, data.frame(method = method, study = study, s = s, i = i,
                                           fdr = 1-n.correct.terms/nrow(identified.vars),
                                           fnr = 1-n.correct.terms/nrow(correct.vars),
                                           n.terms = nrow(identified.vars)))
      }
    }
  }
}

errors2 <- aggregate(cbind(fdr, fnr) ~ method + study + s, data=errors, FUN = function(x) {
  cbind(mean(x), se(x))
})
errors2$fdr.se <- errors2$fdr[,2]; errors2$fdr <- errors2$fdr[,1]
errors2$fnr.se <- errors2$fnr[,2]; errors2$fnr <- errors2$fnr[,1]
n.terms <- aggregate(n.terms ~ method + study + s, data=errors, FUN=median)
errors2 <- merge(errors2, n.terms, by=c("method", "study", "s"), sort=FALSE)
mean.df <- reshape2::melt(cbind(errors2[,1:3], FDR=errors2$fdr, FNR=errors2$fnr), measure.vars=c("FDR", "FNR"))
se.df <- reshape2::melt(cbind(errors2[,1:3], FDR=errors2$fdr.se, FNR=errors2$fnr.se), measure.vars=c("FDR", "FNR"))
n.terms.df <- reshape2::melt(cbind(errors2[,1:3], FDR=errors2$n.terms, FNR=errors2$n.terms), measure.vars=c("FDR", "FNR"))
errors2 <- merge(mean.df, se.df, by=c("method", "study", "s", "variable"), sort=FALSE)
errors2 <- merge(errors2, n.terms.df, by=c("method", "study", "s", "variable"), sort=FALSE)
colnames(errors2)[4:7] <- c("Error Rate Type", "mean", "se", "n.terms")
errors2$study <- study.string(errors2$study)
errors2$study <- factor(errors2$study, levels = c("Linear", "Interactions I", "Interactions II", "Complex"))
errors2$method[errors2$method == "lr"] <- "Logic Regression"
colnames(errors2)[1] <- "Algorithm"

scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
errors2$snr <- factor(scenarios[errors2$s, "snr"],
                      levels = unique(scenarios[errors2$s, "snr"]),
                      ordered = TRUE)
errors2$N <- paste("N =", scenarios[errors2$s, "N"])
errors2$N <- factor(errors2$N,
                    levels = unique(errors2$N),
                    ordered = TRUE)

library(ggh4x)

p <- ggplot(data = errors2, aes(x=interaction(.data[["Error Rate Type"]], snr, sep = "!"), y=.data[["mean"]], fill=Algorithm, color=Algorithm)) +
  geom_point(position = position_dodge(width=.75)) +
  geom_errorbar(aes(ymin = pmax(0,mean-1.96*se), ymax = pmin(1,mean+1.96*se)),
                width=0.75,
                position = position_dodge(width=.75)) +
  geom_text(aes(label=n.terms, y = pmin(1,mean+1.96*se), vjust = -0.5),
            position = position_dodge(width=.75)) +
  scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = "SNR [Signal-to-Noise Ratio]") +
  scale_y_continuous(latex2exp::TeX("False Discovery Rate / False Negative Rate"),
                     expand = expansion(mult = 0.2),
                     breaks = seq(0, 1, by=0.25)) +
  geom_vline(xintercept = c(2.5,4.5)) +
  geom_vline(xintercept = c(1.5,3.5,5.5), linetype="dashed", alpha=.5) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position="bottom") +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0))) +
  facet_grid(cols = vars(N), rows = vars(study), scales = "free")
p
ggsave(plot = p, filename = "optimal_terms.pdf", width = 8, height = 8)

#########################################

# Number of identified terms

n.terms <- NULL
for(method in c("BITS", "glinternet", "lr")) {
  for(study in c("linear", "interactions_no_E", "interactions_E", "complex")) {
    terms_file <- paste(study, "_terms.RData", sep="")
    load(terms_file)

    for(s in 1:9) {
      tmp <- terms[[method]][[s]]
      for(i in 1:100) {
        identified.vars <- tmp[[i]]$optimal.terms
        if(is.null(identified.vars)) next
        n.terms <- rbind(n.terms, data.frame(method = method, study = study, s = s, i = i,
                                             n.terms = nrow(identified.vars)))
      }
    }
  }
}

n.terms$study <- study.string(n.terms$study)
n.terms$study <- factor(n.terms$study, levels = c("Linear", "Interactions I", "Interactions II", "Complex"))
n.terms$method[n.terms$method == "lr"] <- "Logic Regression"
colnames(n.terms)[1] <- "Algorithm"
scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
n.terms$snr <- factor(scenarios[n.terms$s, "snr"],
                      levels = unique(scenarios[n.terms$s, "snr"]),
                      ordered = TRUE)
n.terms$N <- paste("N =", scenarios[n.terms$s, "N"])
n.terms$N <- factor(n.terms$N,
                    levels = unique(n.terms$N),
                    ordered = TRUE)

p <- ggplot(data = n.terms, aes(x=snr, y=.data[["n.terms"]], fill=Algorithm)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot() +
  scale_x_discrete("SNR [Signal-to-Noise Ratio]", waiver()) +
  scale_y_continuous(latex2exp::TeX("Number of Identified Terms")) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  coord_cartesian(ylim=c(0, 25)) +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(N), rows = vars(study), scales = "free")
p
ggsave(plot = p, filename = "n_terms.pdf", width = 8, height = 8)

#########################################

# Main effects vs. interaction effects:
# False discovery and false negative rates

errors <- NULL
for(method in c("BITS", "glinternet", "lr")) {
  for(study in c("complex")) {
    terms_file <- paste(study, "_terms.RData", sep="")
    load(terms_file)

    if(study == "linear")
      correct.vars <- matrix(1:6, ncol=1,byrow=TRUE)
    if(study == "interactions_no_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,NA,
                               1,2,3), ncol=3,byrow=TRUE)
    if(study == "interactions_E")
      correct.vars <- matrix(c(1,NA,NA,
                               2,3,51),ncol=3,byrow=TRUE)
    if(study == "complex")
      correct.vars <- matrix(c(1,NA,NA,
                               3,NA,NA,
                               4,NA,NA,
                               7,8,9,
                               2,51,NA,
                               5,6,51),ncol=3,byrow=TRUE)

    for(s in 1:9) {
      tmp <- terms[[method]][[s]]
      for(i in 1:100) {
        identified.vars <- tmp[[i]]$optimal.terms
        if(is.null(identified.vars)) next
        if(method %in% c("glinternet", "lr")) identified.vars <- ceiling(identified.vars/2)

        n.vars.identified <- rowSums(!is.na(identified.vars)); n.vars.correct <- rowSums(!is.na(correct.vars))
        main.identified <- identified.vars[n.vars.identified == 1,,drop=FALSE]; main.correct <- correct.vars[n.vars.correct == 1,,drop=FALSE]
        int.identified <- identified.vars[n.vars.identified > 1,,drop=FALSE]; int.correct <- correct.vars[n.vars.correct > 1,,drop=FALSE]

        # For FDR computation, compute number of correctly identified terms with respect to
        # the full true terms matrix, respectively, to not put logic regression at a disadvantage
        # due to not being able to identify interactions with E
        # as, e.g., in identifying SNP2 from the true term SNP2 * E
        # (since E is ignored in the judgment of correct terms)
        main.n.correct <- howManyCorrect(main.identified, correct.vars)
        int.n.correct <- howManyCorrect(int.identified, correct.vars)
        main.fdr <- ifelse(nrow(main.identified) > 0, (nrow(main.identified)-main.n.correct)/nrow(main.identified), 0)
        int.fdr <- ifelse(nrow(int.identified) > 0, (nrow(int.identified)-int.n.correct)/nrow(int.identified), 0)

        # Similarly to above, for FNR comparison, compute number of correctly identified terms with respect to
        # the full identified terms matrix, respectively, to not put logic regression at a disadvantage
        # due to not being able to identify interactions with E
        main.n.correct <- howManyCorrect(identified.vars, main.correct)
        int.n.correct <- howManyCorrect(identified.vars, int.correct)
        main.fnr <- (nrow(main.correct)-main.n.correct)/nrow(main.correct)
        int.fnr <- (nrow(int.correct)-int.n.correct)/nrow(int.correct)

        errors <- rbind(errors, data.frame(method = method, study = study, s = s, i = i,
                                           main.fdr = main.fdr,
                                           int.fdr = int.fdr,
                                           main.fnr = main.fnr,
                                           int.fnr = int.fnr,
                                           main.n.terms = nrow(main.identified),
                                           int.n.terms = nrow(int.identified)))
      }
    }
  }
}

errors2 <- aggregate(cbind(main.fdr, int.fdr, main.fnr, int.fnr) ~ method + study + s, data=errors, FUN = function(x) {
  cbind(mean(x), se(x))
})
errors2$main.fdr.se <- errors2$main.fdr[,2]; errors2$main.fdr <- errors2$main.fdr[,1]
errors2$int.fdr.se <- errors2$int.fdr[,2]; errors2$int.fdr <- errors2$int.fdr[,1]
errors2$main.fnr.se <- errors2$main.fnr[,2]; errors2$main.fnr <- errors2$main.fnr[,1]
errors2$int.fnr.se <- errors2$int.fnr[,2]; errors2$int.fnr <- errors2$int.fnr[,1]
n.terms <- aggregate(cbind(main.n.terms, int.n.terms) ~ method + study + s, data=errors, FUN=median)
errors2 <- merge(errors2, n.terms, by=c("method", "study", "s"), sort=FALSE)
mean.df <- reshape2::melt(cbind(errors2[,1:3], main.fdr=errors2$main.fdr, int.fdr=errors2$int.fdr, main.fnr=errors2$main.fnr, int.fnr=errors2$int.fnr), measure.vars=c("main.fdr", "int.fdr", "main.fnr", "int.fnr"))
se.df <- reshape2::melt(cbind(errors2[,1:3], main.fdr=errors2$main.fdr.se, int.fdr=errors2$int.fdr.se, main.fnr=errors2$main.fnr.se, int.fnr=errors2$int.fnr.se), measure.vars=c("main.fdr", "int.fdr", "main.fnr", "int.fnr"))
n.terms.df <- reshape2::melt(cbind(errors2[,1:3], main.fdr=errors2$main.n.terms, int.fdr=errors2$int.n.terms, main.fnr=errors2$main.n.terms, int.fnr=errors2$int.n.terms), measure.vars=c("main.fdr", "int.fdr", "main.fnr", "int.fnr"))
errors2 <- merge(mean.df, se.df, by=c("method", "study", "s", "variable"), sort=FALSE)
errors2 <- merge(errors2, n.terms.df, by=c("method", "study", "s", "variable"), sort=FALSE)
colnames(errors2)[5:7] <- c("mean", "se", "n.terms")
splitted.names <- strsplit(as.character(errors2$variable), ".", fixed=TRUE)
splitted.names <- matrix(unlist(splitted.names), ncol=2, byrow=TRUE)
splitted.names[splitted.names == "main"] <- "Main Effects"
splitted.names[splitted.names == "int"] <- "Interaction Effects"
splitted.names[splitted.names == "fdr"] <- "FDR"
splitted.names[splitted.names == "fnr"] <- "FNR"
errors2[["Term Type"]] <- factor(splitted.names[,1],
                                 levels = unique(splitted.names[,1]),
                                 ordered = TRUE)
errors2[["Error Rate Type"]] <- factor(splitted.names[,2],
                                       levels = unique(splitted.names[,2]),
                                       ordered = TRUE)
errors2$method[errors2$method == "lr"] <- "Logic Regression"
colnames(errors2)[1] <- "Algorithm"
scenarios <- expand.grid(snr = c(0.5, 1, 2), N = c(500, 1000, 2000))
errors2$snr <- paste("SNR =", scenarios[errors2$s, "snr"])
errors2$snr <- factor(errors2$snr,
                      levels = unique(errors2$snr),
                      ordered = TRUE)
errors2$N <- paste("N =", scenarios[errors2$s, "N"])
errors2$N <- factor(errors2$N,
                    levels = unique(errors2$N),
                    ordered = TRUE)

p <- ggplot(data = errors2, aes(x=interaction(.data[["Error Rate Type"]], .data[["Term Type"]], sep = "!"), y=.data[["mean"]], fill=Algorithm, color=Algorithm)) +
  geom_point(position = position_dodge(width=.75)) +
  geom_errorbar(aes(ymin = pmax(0,mean-1.96*se), ymax = pmin(1,mean+1.96*se)),
                width=0.5,
                position = position_dodge(width=.75)) +
  geom_text(aes(label=n.terms, y = pmin(1,mean+1.96*se), vjust = -0.5),
            position = position_dodge(width=.75)) +
  scale_x_discrete(guide = guide_axis_nested(delim = "!"), name = NULL) +
  scale_y_continuous(latex2exp::TeX("False Discovery Rate / False Negative Rate"),
                     expand = expansion(mult = 0.15),
                     breaks = seq(0, 1, by=0.25)) +
  geom_vline(xintercept = 2.5) +
  geom_vline(xintercept = c(1.5,3.5), linetype="dashed", alpha=.5) +
  scale_fill_brewer(palette="Set1") +
  scale_color_brewer(palette="Set1") +
  theme_bw() +
  theme(legend.position="bottom") +
  facet_grid(cols = vars(snr), rows = vars(N))
p
ggsave(plot = p, filename = "fdr_main_int.pdf", width = 8, height = 8)
