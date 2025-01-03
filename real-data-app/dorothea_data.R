# Dorothea data can be downloaded here:
# https://doi.org/10.24432/C5NK6X

# Dorothea predictions can be evaluated here:
# https://codalab.lisn.upsaclay.fr/competitions/7363

# Load raw Dorothea data:

mat <- matrix(0, nrow = 1950, ncol = 100000)

train.raw <- readLines("DOROTHEA/DOROTHEA/dorothea_train.data")
for(i in 1:length(train.raw)) {
  ones <- as.integer(unlist(strsplit(train.raw[i], " ")))
  mat[i, ones] <- 1
}

val.raw <- readLines("DOROTHEA/DOROTHEA/dorothea_valid.data")
for(i in 1:length(val.raw)) {
  ones <- as.integer(unlist(strsplit(val.raw[i], " ")))
  mat[i+length(train.raw), ones] <- 1
}

test.raw <- readLines("DOROTHEA/DOROTHEA/dorothea_test.data")
for(i in 1:length(test.raw)) {
  ones <- as.integer(unlist(strsplit(test.raw[i], " ")))
  mat[i+length(train.raw)+length(val.raw), ones] <- 1
}

y <- integer()

train.labels <- readLines("DOROTHEA/DOROTHEA/dorothea_train.labels")
train.labels <- as.integer(train.labels == "1")
y <- c(y, train.labels)

val.labels <- readLines("DOROTHEA/dorothea_valid.labels")
val.labels <- as.integer(val.labels == "1")
y <- c(y, val.labels)

y2 <- c(y, rep(NA, length(test.raw)))
dat <- data.frame(y = y2, mat)
train <- 1:800
val <- 801:1150
test <- 1151:1950

save(dat, train, val, test, file = "dorothea_neurips.RData")

###################################################

# Select most important predictors using lasso:

load("dorothea_neurips.RData")

model <- glmnet::glmnet(as.matrix(dat[train,-1]), dat$y[train], family = "binomial", alpha = 1,
                     lambda.min.ratio = 1e-4)
apply(model$beta, 2, function(x) sum(x != 0))
vars <- names(which(model$beta[,ncol(model$beta)] != 0))

dat <- dat[,c("y", vars)]
save(dat, train, val, test, file = "dorothea_neurips_lasso.RData")
