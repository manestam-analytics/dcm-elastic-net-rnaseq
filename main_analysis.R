## Script 1 – Primary Objective
## Elastic Net Model for Predicting DCM from RNA-seq Data
##
## This script reproduces the primary analysis reported in the thesis.

rm(list = ls())
graphics.off()

## -----------------------------
## Load data and packages
## -----------------------------

Data <- read.csv("Data.csv", header = FALSE)

library(glmnet)
library(caret)
library(pROC)
library(dplyr)

## -----------------------------
## Data preparation
## -----------------------------

patient_ids <- as.character(unlist(Data[1, -1]))

diagnosis_raw <- as.character(unlist(Data[2, -1]))

y <- factor(
  ifelse(diagnosis_raw == "Dilated cardiomyopathy (DCM)", "DCM", "Control"),
  levels = c("Control", "DCM")
)

print(table(y))

gene_names <- as.character(Data[3:nrow(Data), 1])

expr <- Data[3:nrow(Data), -1]
expr <- apply(expr, 2, as.numeric)

rownames(expr) <- gene_names
colnames(expr) <- patient_ids

## Transpose so rows = patients and columns = genes
X <- t(expr)

## Basic checks
stopifnot(nrow(X) == length(y))
stopifnot(all(rownames(X) == patient_ids))

## Remove near-zero variance genes
nzv <- nearZeroVar(X)
if (length(nzv) > 0) {
  X <- X[, -nzv]
}

## -----------------------------
## Train/test split
## -----------------------------

set.seed(123)

train_idx <- createDataPartition(
  y,
  p = 0.7,
  list = FALSE
)

X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]

y_train <- y[train_idx]
y_test  <- y[-train_idx]

## Check for sample overlap
any(rownames(X_train) %in% rownames(X_test))

## -----------------------------
## Elastic net model
## -----------------------------

set.seed(123)

cv_fit <- cv.glmnet(
  x = as.matrix(X_train),
  y = y_train,
  family = "binomial",
  alpha = 0.5,
  standardize = TRUE,
  type.measure = "auc"
)

## -----------------------------
## Prediction
## -----------------------------

pred_prob <- predict(
  cv_fit,
  newx = as.matrix(X_test),
  s = "lambda.min",
  type = "response"
)[, 1]

## -----------------------------
## ROC AUC
## -----------------------------

roc_obj <- roc(
  response = y_test,
  predictor = pred_prob,
  levels = c("Control", "DCM"),
  direction = "<"
)

auc_value <- auc(roc_obj)
print(auc_value)

## Bootstrap confidence interval for AUC
set.seed(123)

auc_ci <- ci.auc(
  roc_obj,
  method = "bootstrap",
  boot.n = 1000
)

print(auc_ci)

## Save AUC results
auc_results <- data.frame(
  AUC = as.numeric(auc_value),
  CI_lower = as.numeric(auc_ci[1]),
  CI_upper = as.numeric(auc_ci[3])
)

write.csv(
  auc_results,
  "AUC_results_test_set.csv",
  row.names = FALSE
)

## -----------------------------
## ROC curve
## -----------------------------

png(
  "ROC_curve_test_set.png",
  width = 2400,
  height = 2400,
  res = 300
)

plot(
  roc_obj,
  col = "#1B6CA8",
  lwd = 3,
  legacy.axes = TRUE,
  main = paste0(
    "ROC Curve - Test Set\n",
    "AUC = ", round(auc_value, 4),
    " (95% CI: ",
    round(auc_ci[1], 4),
    " - ",
    round(auc_ci[3], 4),
    ")"
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2,
  col = "gray"
)

dev.off()

## -----------------------------
## Calibration analysis
## -----------------------------

calibration_data <- data.frame(
  patient_id = rownames(X_test),
  observed = as.numeric(y_test == "DCM"),
  predicted_probability = pred_prob
)

write.csv(
  calibration_data,
  "calibration_predictions_test_set.csv",
  row.names = FALSE
)

## Calibration slope
cal_model <- glm(
  observed ~ qlogis(predicted_probability),
  data = calibration_data,
  family = binomial
)

summary(cal_model)

## Calibration plot using deciles
cal_data <- data.frame(
  pred = pred_prob,
  obs = as.numeric(y_test == "DCM")
)

cal_data$bin <- cut(
  cal_data$pred,
  breaks = quantile(cal_data$pred, probs = seq(0, 1, 0.1)),
  include.lowest = TRUE
)

cal_plot <- cal_data %>%
  group_by(bin) %>%
  summarise(
    mean_pred = mean(pred),
    mean_obs = mean(obs),
    .groups = "drop"
  )

write.csv(
  cal_plot,
  "calibration_deciles_test_set.csv",
  row.names = FALSE
)
