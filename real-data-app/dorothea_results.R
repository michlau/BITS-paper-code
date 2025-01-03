# Process Dorothea results here locally

# Write test data predictions to .predict file, as required by the evaluation
# system on https://codalab.lisn.upsaclay.fr/competitions/7363

# These need to be zipped (per method) with the dummy predictions for the other
# data sets from the competition created below and submitted to the website linked above

load("dorothea_scores.RData")

preds.BITS <- scores$BITS[[1]]$preds
write.table(preds.BITS, file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

preds.glinternet <- scores$glinternet[[1]]$preds
write.table(preds.glinternet, file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

preds.rulefit <- scores$rulefit[[1]]$preds
write.table(preds.rulefit, file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

preds.mars <- scores$mars[[1]]$preds
write.table(preds.mars, file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

preds.rf <- scores$rf[[1]]$preds
write.table(preds.rf, file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

preds.lr <- scores$lr[[1]]$preds
write.table(round(preds.lr, digits=5), file = "dorothea_test.predict", col.names = FALSE, row.names = FALSE)

# Create dummy prediction files for the other data sets of the competition:
preds.arcene <- runif(700)
write.table(preds.arcene, file = "arcene_test.predict", col.names = FALSE, row.names = FALSE)

preds.gisette <- runif(6500)
write.table(preds.gisette, file = "gisette_test.predict", col.names = FALSE, row.names = FALSE)

preds.dexter <- runif(2000)
write.table(preds.dexter, file = "dexter_test.predict", col.names = FALSE, row.names = FALSE)

preds.madelon <- runif(1800)
write.table(preds.madelon, file = "madelon_test.predict", col.names = FALSE, row.names = FALSE)

###########################################

# Model interpretation

# BITS:
model <- scores$BITS[[1]]$model
vars <- BITS::get.included.vars(model)
vars
nrow(vars) # Number of terms
sum(is.na(vars[,2])) # Number of main effects
sum(!is.na(vars[,2])) # Number of second order interaction effects
print(model, latex=TRUE)

# glinternet:
model <- scores$glinternet[[1]]$model
# As in glinternet.terms (common/prep.R) to obtain the contained terms:
load("dorothea_params.RData")
s <- results$glinternet[[1]]$scoreSummary$lambda
ll <- model$lambda
ind <- which.min(abs(ll - s))
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

# Random Forests:
model <- scores$rf[[1]]$model
# Distribution of number of leaves per tree:
summary(as.numeric(lapply(model$forest$terminal.class.counts, length)))
