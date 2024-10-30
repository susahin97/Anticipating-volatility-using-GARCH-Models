################################################################################
## Variable Selection 
################################################################################
## Libraries
library(tseries) # for adf.test
library(readr) # for read_csv
library(dplyr) # for data manipulation
library(forecast)
library(corrplot)
library(glmnet)
library(caret)
library(ggcorrplot)
library(GGally)
library(ggplot2)
library(lars)
library(leaps)
library(MASS)
library(gridExtra)
library(patchwork)
library(caTools)
library(caret)
library(broom)
library(olsrr)
library(cowplot)
library(pROC)

## Data
str(selection)
skimr::skim(selection)

################################################################################
## feature selection for WTI returns 
################################################################################
## calculate returns ----------------------------------------------------------
selection$wti_r <- 100 * c(NA, diff(log(selection$WTI)))
selection$brent_r <- 100 * c(NA, diff(log(selection$Brent)))
selection$baid_r <- 100 * c(NA, diff(log(selection$BAID)))
selection$dxy_r <- 100 * c(NA, diff(log(selection$DXY)))
selection$sp500_r <- 100 * c(NA, diff(log(selection$SP500)))
selection$ism_r <- 100 * c(NA, diff(log(selection$ISM)))
selection$steel_r <- 100 * c(NA, diff(log(selection$steel)))
selection$surplus_r <- 100 * c(NA, diff(log(selection$Opec_surplus_capacity)))
selection$production_r <- 100 * c(NA, diff(log(selection$crude_production)))
selection$import_r <- 100 * c(NA, diff(log(selection$crude_import)))
selection$export_r <- 100 * c(NA, diff(log(selection$crude_export)))
selection$kilian_r <- 100 * c(NA, diff(log(selection$Kilian)))

# remove NA from the returns
selection <- selection[-1,]
selection$kilian_r <- na_locf(selection$kilian_r)

## stationary check ------------------------------------------------------------
# Iterate through each column except for the Date column
for (col_name in names(selection)[-1]) { 
  cat("Performing ADF Test on:", col_name, "\n")
  
  # Extract the column by name
  column_data <- selection[[col_name]]
  
  # Perform the ADF test
  adf_test_result <- adf.test(column_data, alternative = "stationary")
  
  # Print the results
  print(adf_test_result)
  cat("\n") # Just for better readability of output
}

## Not Stationary -------------------------
# WTI_r: p-value = 0.01  - stationary
# Brent_r : p-value = 0.01 - stationary
# SP500_r: p-value = 0.01 - stationary
# Dollar_index_r: p-value = 0.01 - stationary
# BAID_r: p-value = 0.01 - Stationary
# ISM_r: p-value = 0.01 - Stationary
# steel_r: p-value = 0.01 - Stationary
# Opec_surpls_capacity: p-value = 0.01 - Stationary
# crude_production: p-value = 0.01 -  stationary
# crude_import: p-value = 0.01 - stationary
# crude_export: p-value = 0.01 - stationary 
# Kilian: p-value = 0.2981 - Not stationary


# List of non-stationary variables identified from ADF test results
non_stationary_vars <- c("Kilian")

# Apply first-order differencing and update the dataframe
for (var in non_stationary_vars) {
  # Apply differencing
  differenced_series <- diff(selection[[var]], differences = 1)
  
  # Correctly update the dataframe with the differenced series
  selection[[var]] <- c(NA, differenced_series)
}

# After updating all variables, remove the first row from the dataframe
selection <- selection[-1,]


# Reiterate through each variable for ADF test
for (var in non_stationary_vars) { 
  cat("Rechecking ADF Test for differenced:", var, "\n")
  
  # Perform the ADF test on the differenced data
  adf_test_result <- adf.test(selection[[var]], alternative = "stationary")
  
  # Print the ADF test results
  print(adf_test_result)
  cat("\n") 
} # all variables are now stationary


## Corrplot for WTI ------------------------------------------------------------
# Exclude the Brent
selected_vars <- selection[, !(names(selection) %in% c("Brent", "Date","WTI",
                                                       "Kilian","SP500","DXY",
                                                       "BAID", "ISM", "steel",
                                                       "Opec_surplus_capacity",
                                                       "crude_production",
                                                       "crude_import", "crude_export", "brent_r"))]

# Calculate the correlation matrix
corr_matrix <- cor(selected_vars, use = "complete.obs")

corrplot(corr_matrix, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45,
         addCoef.col = "black", # Add correlation coefficients
         type = "upper",        # Show upper half only
         diag = FALSE) 
# wti/ ism : 0.45
# opec_surplus/production : -0.4

# Split the data into training and test sets -----------------------------------
set.seed(123) # For reproducibility
sample_size <- floor(0.7 * nrow(selected_vars))
training_index <- sample(seq_len(nrow(selected_vars)), size = sample_size) 
train_data <- selected_vars[training_index, ]
test_data <- selected_vars[-training_index, ]

# Standardize the variables
scalingData <- preProcess(train_data[,1:10],method = c("center", "scale")) 
train_data[,1:10] <- predict(scalingData,train_data[,1:10]) 
test_data[,1:10] <- predict(scalingData,test_data[,1:10])

x_train <- model.matrix(wti_r~.-1,data=train_data)
y_train <- train_data$wti_r

x_test <- model.matrix(wti_r~.-1, data=test_data)
y_test <- test_data$wti_r

## correlation plot for train set ----------------
corr_matrix <- cor(train_data, use = "complete.obs")

wti <- corrplot(corr_matrix, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45,
         addCoef.col = "black", # Add correlation coefficients
         type = "upper",        # Show upper half only
         diag = FALSE) 
# Pairs Plot -------------------------------------------------------------------
ggpairs(train_data, lower = list(continuous = wrap("points", size = 0.5)),
        upper = list(continuous = wrap("cor", digits = 2, size = 3)))


## Box Plot of vars ------------------------------------------------------------
# Scale the data
box <- as.data.frame(train_data[, -11])

# Create boxplots
ggplot(stack(box), aes(x = ind, y = values)) +
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90)) +
  labs(title = "Boxplots of Scaled Variables", x = "Variable", y = "Scaled Value")

outliers <- boxplot(train_data[,-11], plot = FALSE)$out
if (length(outliers) > 0) {
  print("Potential outliers found:")
  print(outliers)
} else {
  print("No potential outliers found.")
}

##### Collinearity #############################################################
X = as.matrix(train_data)
e = eigen(t(X) %*% X)$val
max(e)/min(e) 

model.ols <- lm(wti_r~., data=train_data) 

## TOLERANCE & VARIANCE INFLATION FACTOR (VIF)
ols_vif_tol(model.ols)
vif_results <- ols_vif_tol(model.ols)
variables <- vif_results[vif_results$Tolerance < 0.1 & vif_results$VIF > 10, "Variables"]
variables # THERE ARE NO HIGHLY COLLINEAR VARIABLES


## lINEAR REGRESSION ###########################################################
model <- lm(y_train ~ x_train - 1) 
summary(model)

# DXY_r, SP500_r, ISM_r and steel_r are significant

model_2 <- lm(y_train ~ dxy_r + sp500_r + ism_r + steel_r - 1, data=train_data)
summary(model_2)

# Making predictions
predictions_train <- predict(model_2, newdata=train_data)
predictions_test <- predict(model_2, newdata=test_data)

# Calculate MAE for training and testing sets
mae_train <- mean(abs(predictions_train - y_train))
mae_test <- mean(abs(predictions_test - y_test))

print(paste("MAE Train: ", mae_train))
print(paste("MAE Test: ", mae_test))

### Regularization Methods #####################################################
## 1) RIDGE ####################################################################
fit_ridge <- glmnet(x_train,
                    y_train,
                    family = "gaussian",
                    type.measure = "mse",
                    alpha = 0, 
                    lambda = NULL)

sprintf("Max Lambda: %.3f", max(fit_ridge$lambda))
sprintf("Min Lambda: %.3f", min(fit_ridge$lambda))

# Example of min coeff
coef(fit_ridge)[c("baid_r", "sp500_r"),1]

# Example of max coeff
coef(fit_ridge)[c("baid_r", "sp500_r"),100]

# cv.glmnet 
set.seed(123)
lambdas_to_try <- 10^seq(4, -4, length = 100)
ridge_cv = cv.glmnet(x_train,
                     y_train,
                     alpha=0,
                     family="gaussian",
                     type.measure = "mse",
                     lambda=lambdas_to_try)
ridge_cv

ridge_cv$lambda.min 
ridge_cv$lambda.1se 

plot(fit_ridge, xvar = "lambda", label = TRUE)
abline(v = log(ridge_cv$lambda.min),col="blue", lty = 5)
abline(v = log(ridge_cv$lambda.1se),col="red", lty = 5)
title("Ridge Regression: Log Lambda vs. Coefficient Shrinkage")

# coefficients with lambda.min and lambda.1se
round(cbind(
  coef(ridge_cv, s = 'lambda.min'),
  coef(ridge_cv, s = 'lambda.1se')),   
  digits = 3 )

# Plot cross-validation results
plot(ridge_cv)
abline( h = ridge_cv$cvup[ridge_cv$index[1]], lty = 5, col="red")
title("Ridge Regression: Log Lambda vs. Misclassification Error")


# Ridge regression with lambda.min
fit.ridge.min.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 0, 
                               family="gaussian", 
                               type.measure= "mse",
                               lambda = ridge_cv$lambda.min)
coef.ridge.min <- coef(fit.ridge.min.lambda)

# Ridge regression with lambda.1se
fit.ridge.1se.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 0,
                               family="gaussian",
                               type.measure= "mse",
                               lambda = ridge_cv$lambda.1se)
coef.ridge.1se <- coef(fit.ridge.1se.lambda)


# Generate predictions for test data using lambda.min
predictions.ridge.min <- predict(fit.ridge.min.lambda, newx = x_train, type="response")

postResample(pred = predictions.ridge.min, obs = y_train)
#    RMSE      Rsquared    MAE 
# 0.9010129   0.2496606   0.6362903 

# Generate predictions for test data using lambda.1se
predictions.ridge.1se <- predict(fit.ridge.1se.lambda, newx = x_train, type="response")
postResample(pred = predictions.ridge.1se, obs = y_train)

# --------------------------------------------------------------------------------
## FINAL MODEL
ridge.model <- glmnet(x_train,
                      y_train,
                      alpha = 0,
                      family = "gaussian",
                      type.measure = "mse",
                      lambda = ridge_cv$lambda.min)

predictions.ridge.model <- predict(ridge.model, newx = x_train, type = "response")


ridge.perf <- postResample(pred = predictions.ridge.model, obs = y_train)
ridge.perf
## Variable Importance ---------------------------------------------------------
var_importance <- varImp(ridge.model, ridge_cv$lambda.min)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle("Variable Importance Ridge Regularization")

# gives the same results as linear regression.

# plot of the coefficients
plot(coef.ridge.min, xlab = "Coefficient Index", ylab = "Coefficient Value", main = "Ridge Regression Coefficients")
abline(h = 0,col="red", lty = 2)  # Add a dashed line at 0

# Add labels for observation points
text(seq_along(coef.ridge.min), coef.ridge.min, labels = colnames(x_train), pos = 4, xpd = TRUE)


##### Lasso ####################################################################
fit_lasso <- glmnet(x_train,
                    y_train,
                    alpha = 1,
                    family="gaussian",
                    type.measure = "mse",
                    lambda=NULL)

sprintf("Max Lambda: %.3f", max(fit_lasso$lambda))
sprintf("Min Lambda: %.3f", min(fit_lasso$lambda))



lambdas_to_try <- 10^seq(4, -4, length = 100)
lasso_cv = cv.glmnet(x_train, 
                     y_train,
                     alpha=1, 
                     family="gaussian", 
                     type.measure="mse",
                     lambda = lambdas_to_try )
lasso_cv

plot(fit_lasso, xvar = "lambda", label = TRUE)
abline(v = log(lasso_cv$lambda.min),col="blue", lty = 5)
abline(v = log(lasso_cv$lambda.1se),col="red", lty = 5)

lasso_cv$lambda.min 
lasso_cv$lambda.1se 

# coefficients with lambda.min and lambda.1se
round(cbind(
  coef(lasso_cv, s = 'lambda.min'),
  coef(lasso_cv, s = 'lambda.1se')),   
  digits = 3 )


# Plot cross-validation results
plot(lasso_cv)
abline( h = lasso_cv$cvup[lasso_cv$index[1]], lty = 4 )

## Lasso Fit -------------------------------------------------------------------
# Lasso regression with lambda.min
fit.lasso.min.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 1, 
                               family="gaussian", 
                               type.measure= "mse",
                               lambda = lasso_cv$lambda.min)
coef.lasso.min <- coef(fit.lasso.min.lambda)

# Lasso regression with lambda.1se
fit.lasso.1se.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 1,
                               family="gaussian",
                               type.measure="mse",
                               lambda = lasso_cv$lambda.min)
coef.lasso.1se <- coef(fit.lasso.1se.lambda)

# Generate predictions for test data using lambda.min
predictions.lasso.min <- predict(fit.lasso.min.lambda, newx = x_train, type = "response")

postResample(pred = predictions.lasso.min, obs = y_train)

# Generate predictions for test data using lambda.1se
predictions.lasso.1se <- predict(fit.lasso.1se.lambda, newx = x_train, type="response")

postResample(pred = predictions.ridge.min, obs = y_train)


## final model-------------------------------------------------------------------------------
lasso.model <- glmnet(x_train,
                      y_train, 
                      family="gaussian",
                      alpha = 1, lambda = lasso_cv$lambda.min)

coef.lasso <- coef(lasso.model)

predictions.lasso.model <- predict(lasso.model, newx = x_train, type="response")
lasso.perf <- postResample(pred = predictions.lasso.model, obs = y_train)

var_importance <- varImp(lasso.model, lasso_cv$lambda.min)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle(" WTI Variable Importance Lasso Regularization")
#Lasso also gives the same results.

## plotting coefficiennts --------------------------------------------------------

plot(coef.lasso, xlab = "Coefficient Index", ylab = "Coefficient Value", main = "Lasso Regression Coefficients")
abline(h = 0,col="red", lty = 2) 
text(seq_along(coef.ridge.min), coef.ridge.min, labels = colnames(x_train), pos = 4, xpd = TRUE)

## random forest ---------------------------------------------------------------
# Load the randomForest package
library(randomForest)

# Set seed for reproducibility
set.seed(123)

# Train the Random Forest model
rf_model <- randomForest(x_train, y_train, ntree=500, importance=TRUE)

# Predict using the trained model on the test data
predictions <- predict(rf_model, x_test)

# You can also evaluate the model's performance on the test set. For a regression problem, you might use RMSE (Root Mean Squared Error)
actual_vs_predicted <- data.frame(Actual = y_test, Predicted = predictions)
rmse <- sqrt(mean((actual_vs_predicted$Actual - actual_vs_predicted$Predicted)^2))
print(paste("RMSE on test data:", rmse))

# Optionally, print the model's importance
print(importance(rf_model))
var_importance <- varImp(rf_model)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle(" WTI Variable Importance Random Forest")

varImpPlot(rf_model)
################################################################################
## feature selection for Brent retuns 
################################################################################
## Corrplot for brent returns ------------------------------------------------------------
# Exclude the Brent
selected_vars <- selection[, !(names(selection) %in% c("Brent", "Date","WTI",
                                                       "Kilian","SP500","DXY",
                                                       "BAID", "ISM", "steel",
                                                       "Opec_surplus_capacity",
                                                       "crude_production",
                                                       "crude_import", "crude_export", "wti_r"))]

# Calculate the correlation matrix
corr_matrix <- cor(selected_vars, use = "complete.obs")

corrplot(corr_matrix, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45,
         addCoef.col = "black", # Add correlation coefficients
         type = "upper",        # Show upper half only
         diag = FALSE) 

# ism_r / brent_r: 0.42
# surplus_r / production_r : -0.4

# Split the data into training and test sets -----------------------------------
set.seed(123) # For reproducibility
sample_size <- floor(0.7 * nrow(selected_vars))
training_index <- sample(seq_len(nrow(selected_vars)), size = sample_size) 
train_data <- selected_vars[training_index, ]
test_data <- selected_vars[-training_index, ]

# Standardize the variables
scalingData <- preProcess(train_data[,1:10],method = c("center", "scale")) 
train_data[,1:10] <- predict(scalingData,train_data[,1:10]) 
test_data[,1:10] <- predict(scalingData,test_data[,1:10])

x_train <- model.matrix(brent_r~.-1,data=train_data)
y_train <- train_data$brent_r

x_test <- model.matrix(brent_r~.-1, data=test_data)
y_test <- test_data$brent_r

## correlation plot for train set ----------------
# Calculate the correlation matrix
corr_matrix <- cor(train_data, use = "complete.obs")

brent <- corrplot(corr_matrix, method = "shade", shade.col = NA, tl.col = "black", tl.srt = 45,
         addCoef.col = "black", # Add correlation coefficients
         type = "upper",        # Show upper half only
         diag = FALSE) 
# brent retuns - DXY_r : -0.34


# Pairs Plot -------------------------------------------------------------------
ggpairs(train_data, lower = list(continuous = wrap("points", size = 0.5)),
        upper = list(continuous = wrap("cor", digits = 2, size = 3)))


## Box Plot of vars ------------------------------------------------------------
# Scale the data
box <- as.data.frame(train_data[, -11])

# Create boxplots
ggplot(stack(box), aes(x = ind, y = values)) +
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90)) +
  labs(title = "Boxplots of Scaled Variables", x = "Variable", y = "Scaled Value")

outliers <- boxplot(train_data[,-11], plot = FALSE)$out
if (length(outliers) > 0) {
  print("Potential outliers found:")
  print(outliers)
} else {
  print("No potential outliers found.")
}

##### Collinearity #############################################################
X = as.matrix(train_data)
e = eigen(t(X) %*% X)$val
max(e)/min(e)  # 9548.968

model.ols <- lm(brent_r~., data=train_data) 

## TOLERANCE & VARIANCE INFLATION FACTOR (VIF)
ols_vif_tol(model.ols)
vif_results <- ols_vif_tol(model.ols)
variables <- vif_results[vif_results$Tolerance < 0.1 & vif_results$VIF > 10, "Variables"]
variables # THERE ARE NO HIGHLY COLLINEAR VARIABLES


## lINEAR REGRESSION ###########################################################
model <- lm(y_train ~ x_train - 1) 
summary(model)

# DXY, ISM and sp500 are significant

model_2 <- lm(y_train ~ dxy_r +  sp500_r + ism_r -1, data=train_data)
summary(model_2)

# Making predictions
predictions_train <- predict(model_2, newdata=train_data)
predictions_test <- predict(model_2, newdata=test_data)

# Calculate MAE for training and testing sets
mae_train <- mean(abs(predictions_train - y_train))
mae_test <- mean(abs(predictions_test - y_test))

print(paste("MAE Train: ", mae_train))
print(paste("MAE Test: ", mae_test))

### Regularization Methods #####################################################
## 1) RIDGE ####################################################################
fit_ridge <- glmnet(x_train,
                    y_train,
                    family = "gaussian",
                    type.measure = "mse",
                    alpha = 0, 
                    lambda = NULL)

sprintf("Max Lambda: %.3f", max(fit_ridge$lambda))
sprintf("Min Lambda: %.3f", min(fit_ridge$lambda))

# Example of min coeff
coef(fit_ridge)[c("dxy_r", "sp500_r"),1]

# Example of max coeff
coef(fit_ridge)[c("dxy_r", "sp500_r"),100]

# cv.glmnet 
set.seed(123)
lambdas_to_try <- 10^seq(4, -4, length = 100)
ridge_cv = cv.glmnet(x_train,
                     y_train,
                     alpha=0,
                     family="gaussian",
                     type.measure = "mse",
                     lambda=lambdas_to_try)
ridge_cv

ridge_cv$lambda.min 
ridge_cv$lambda.1se 

plot(fit_ridge, xvar = "lambda", label = TRUE)
abline(v = log(ridge_cv$lambda.min),col="blue", lty = 5)
abline(v = log(ridge_cv$lambda.1se),col="red", lty = 5)
title("Ridge Regression: Log Lambda vs. Coefficient Shrinkage")

# coefficients with lambda.min and lambda.1se
round(cbind(
  coef(ridge_cv, s = 'lambda.min'),
  coef(ridge_cv, s = 'lambda.1se')),   
  digits = 3 )

# Plot cross-validation results
plot(ridge_cv)
abline( h = ridge_cv$cvup[ridge_cv$index[1]], lty = 5, col="red")
title("Ridge Regression: Log Lambda vs. Misclassification Error")


# Ridge regression with lambda.min
fit.ridge.min.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 0, 
                               family="gaussian", 
                               type.measure= "mse",
                               lambda = ridge_cv$lambda.min)
coef.ridge.min <- coef(fit.ridge.min.lambda)

# Ridge regression with lambda.1se
fit.ridge.1se.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 0,
                               family="gaussian",
                               type.measure= "mse",
                               lambda = ridge_cv$lambda.1se)
coef.ridge.1se <- coef(fit.ridge.1se.lambda)


# Generate predictions for test data using lambda.min
predictions.ridge.min <- predict(fit.ridge.min.lambda, newx = x_train, type="response")

postResample(pred = predictions.ridge.min, obs = y_train)
#    RMSE      Rsquared    MAE 
# 0.8856237    0.2487546  0.6503031  

# Generate predictions for test data using lambda.1se
predictions.ridge.1se <- predict(fit.ridge.1se.lambda, newx = x_train, type="response")
postResample(pred = predictions.ridge.1se, obs = y_train)

# --------------------------------------------------------------------------------

## FINAL MODEL
ridge.model <- glmnet(x_train,
                      y_train,
                      alpha = 0,
                      family = "gaussian",
                      type.measure = "mse",
                      lambda = ridge_cv$lambda.min)

predictions.ridge.model <- predict(ridge.model, newx = x_train, type = "response")


ridge.perf <- postResample(pred = predictions.ridge.model, obs = y_train)
ridge.perf

## Variable Importance ---------------------------------------------------------
var_importance <- varImp(ridge.model, ridge_cv$lambda.min)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle("Variable Importance Ridge Regularization")

# DXY, ISM, SP500

# plot of the coefficients
plot(coef.ridge.min, xlab = "Coefficient Index", ylab = "Coefficient Value", main = "Ridge Regression Coefficients")
abline(h = 0,col="red", lty = 2)  # Add a dashed line at 0

# Add labels for observation points
text(seq_along(coef.ridge.min), coef.ridge.min, labels = colnames(x_train), pos = 4, xpd = TRUE)

##### Lasso ####################################################################
fit_lasso <- glmnet(x_train,
                    y_train,
                    alpha = 1,
                    family="gaussian",
                    type.measure = "mse",
                    lambda=NULL)

sprintf("Max Lambda: %.3f", max(fit_lasso$lambda))
sprintf("Min Lambda: %.3f", min(fit_lasso$lambda))



lambdas_to_try <- 10^seq(4, -4, length = 100)
lasso_cv = cv.glmnet(x_train, 
                     y_train,
                     alpha=1, 
                     family="gaussian", 
                     type.measure="mse",
                     lambda = lambdas_to_try )
lasso_cv

plot(fit_lasso, xvar = "lambda", label = TRUE)
abline(v = log(lasso_cv$lambda.min),col="blue", lty = 5)


lasso_cv$lambda.min 
lasso_cv$lambda.1se 

# coefficients with lambda.min and lambda.1se
round(cbind(
  coef(lasso_cv, s = 'lambda.min'),
  coef(lasso_cv, s = 'lambda.1se')),   
  digits = 3 )


# Plot cross-validation results
plot(lasso_cv)
abline( h = lasso_cv$cvup[ridge_cv$index[1]], lty = 4 )

## Lasso Fit -------------------------------------------------------------------
# Lasso regression with lambda.min
fit.lasso.min.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 1, 
                               family="gaussian", 
                               type.measure= "mse",
                               lambda = lasso_cv$lambda.min)
coef.lasso.min <- coef(fit.lasso.min.lambda)

# Lasso regression with lambda.1se
fit.lasso.1se.lambda <- glmnet(x_train,
                               y_train,
                               alpha = 1,
                               family="gaussian",
                               type.measure="mse",
                               lambda = lasso_cv$lambda.min)
coef.lasso.1se <- coef(fit.lasso.1se.lambda)

# Generate predictions for test data using lambda.min
predictions.lasso.min <- predict(fit.lasso.min.lambda, newx = x_train, type = "response")

postResample(pred = predictions.lasso.min, obs = y_train)

# Generate predictions for test data using lambda.1se
predictions.lasso.1se <- predict(fit.lasso.1se.lambda, newx = x_train, type="response")

postResample(pred = predictions.ridge.min, obs = y_train)

## final model-------------------------------------------------------------------------------
lasso.model <- glmnet(x_train,
                      y_train, 
                      family="gaussian",
                      alpha = 1, lambda = lasso_cv$lambda.min)

coef.lasso <- coef(lasso.model)

predictions.lasso.model <- predict(lasso.model, newx = x_train, type="response")
lasso.perf <- postResample(pred = predictions.lasso.model, obs = y_train)

var_importance <- varImp(lasso.model, lasso_cv$lambda.min)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle(" Brent Variable Importance Lasso Regularization")

# sp500, DXY, ISM + import_r and steel_r

## plotting coefficiennts --------------------------------------------------------
plot(coef.lasso, xlab = "Coefficient Index", ylab = "Coefficient Value", main = "Lasso Regression Coefficients")
abline(h = 0,col="red", lty = 2) 
text(seq_along(coef.lasso.min), coef.lasso.min, labels = colnames(x_train), pos = 4, xpd = TRUE)

## random forest ---------------------------------------------------------------
# Load the randomForest package
library(randomForest)

# Set seed for reproducibility
set.seed(123)

# Train the Random Forest model
rf_model <- randomForest(x_train, y_train, ntree=500, importance=TRUE)

# Predict using the trained model on the test data
predictions <- predict(rf_model, x_test)

# You can also evaluate the model's performance on the test set. For a regression problem, you might use RMSE (Root Mean Squared Error)
actual_vs_predicted <- data.frame(Actual = y_test, Predicted = predictions)
rmse <- sqrt(mean((actual_vs_predicted$Actual - actual_vs_predicted$Predicted)^2))
print(paste("RMSE on test data:", rmse))

# Optionally, print the model's importance
print(importance(rf_model))
var_importance <- varImp(rf_model)

ggplot(var_importance, aes(x = Overall, y = reorder(rownames(var_importance), -Overall))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(x = "Importance", y = "Variable") +
  ggtitle(" WTI Variable Importance Random Forest")

varImpPlot(rf_model)
