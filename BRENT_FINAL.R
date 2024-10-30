################################################################################
## BRENT 
################################################################################
# LIBRARIES --------------------------------------------------------------------
library(nortsTest)
library(fGarch)
library(rugarch)
library(quantmod)
library(xts)
library(tseries)
library(urca)
library(ggplot2)
library(zoo)
library(QRM)
library(fracdiff)
library(vrtest)
library(forecast)
library(lmtest)
library(sandwich)
library(pracma)
library(gridExtra)
library(readr)
library(garma)
library(tidyverse)
library(lubridate)
library(imputeTS)
library(strucchange)
library(arfima)
library(readxl)
library(parallel)
library(doParallel)
library(MTS)
library(PerformanceAnalytics)
library(segMGarch)

## Download the Data -----------------------------------------------------------
SP500 <- read_csv("SP500.csv",col_types = cols(Date = col_date(format = "%d/%m/%Y")))
Dollar_index <- read_csv("Dollar_index.csv",col_types = cols(Date = col_date(format = "%d/%m/%Y")))
Brent_price <- read_csv("Brent_price.csv", col_types = cols(Date = col_date(format = "%d/%m/%Y")))
WTI_price <- read_csv("wti_price.csv", col_types = cols(Date = col_date(format = "%d/%m/%Y")))
OPEC_meetings <- read_excel("OPEC meetings.xlsx")

# Convert POSIXct to Date in OPEC_meetings
OPEC_meetings$Date <- as.Date(OPEC_meetings$Date)

## Data Exploration ------------------------------------------------------------
ggplot(Brent_price, aes(x = Date, y = Brent)) +
  geom_line(color = "blue") + 
  theme_minimal() +  
  labs(
    title = "Daily Brent Crude Oil Closing Price",
    x = "Date",
    y = "Price"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),  
    legend.position = "topright",  
    legend.title = element_blank() 
  ) +
  guides(
    color = guide_legend("Brent Daily Prices") 
  )

qqnorm(Brent_price$Brent,main="Brent SPOT PRICES")
qqline(Brent_price$Brent)

## Obtain Returns --------------------------------------------------------------
# Calculate log returns according to the formula
Brent_price$Returns <- 100 * c(NA, diff(log(Brent_price$Brent)))

# Plotting Brent Returns
Brent_return_p <- ggplot(Brent_price, aes(x = Date, y = Returns)) +
  geom_line(color = "blue", na.rm = TRUE) +
  theme_minimal() +
  labs(title = "Daily Brent Crude Oil Returns", x = "Date", y = "Returns") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(-30, 30))  

# Plotting histogram and distribution ------------------------------------------
par(mfrow = c(1, 2))

# Calculate mean and standard deviation of log returns
meanreturn <- mean(na.omit(Brent_price$Returns), na.rm = TRUE)
sdreturn <- sd(na.omit(Brent_price$Returns), na.rm = TRUE)

# Plot the histogram and the normal distribution
hist(Brent_price$Returns, breaks = 40, freq = FALSE, main = "Log Return Histogram", xlab = "Log Returns", col = "lightblue")
curve(dnorm(x, mean = meanreturn, sd = sdreturn), add = TRUE, col = "red", lwd = 2)

# Plotting the density of log returns
plot(density(na.omit(Brent_price$Returns)), main = "Log Return Empirical Distribution", 
     xlab = "Log Returns", ylab = "Density", col = "blue")
curve(dnorm(x, mean = meanreturn, sd = sdreturn), 
      from = min(na.omit(Brent_price$Returns)), 
      to = max(na.omit(Brent_price$Returns)), 
      add = TRUE, col = "red", lwd = 2)

# Plotting the density of log returns
plot(density(na.omit(Brent_price$Returns)), 
     main = "Log Return Empirical Distribution", 
     xlab = "Log Returns", 
     ylab = "Density", 
     col = "blue")

# Plotting squared log returns --------------------------------------------------
Brent_price$SquaredReturns <- Brent_price$Returns^2

# Plotting the squared log returns
ggplot(Brent_price, aes(x = Date, y = SquaredReturns)) +
  geom_line(color = "red", na.rm = TRUE) +  
  theme_minimal() +  
  labs(
    title = "Daily Squared Log Returns of Brent Crude Oil",
    x = "Date",
    y = "Squared Log Returns"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5) 
  )

# Descriptive Stats-------------------------------------------------------------
# Remove NA values before calculations
Brent_returns <- na.omit(Brent_price$Returns)

# Calculate mean
mean_returns <- mean(Brent_returns)

# Calculate variance
variance_returns <- var(Brent_returns)

# Calculate standard deviation
std_dev_returns <- sd(Brent_returns)

# Calculate skewness
skewness_returns <- skewness(Brent_returns)

# Calculate kurtosis
kurtosis_returns <- kurtosis(Brent_returns)

# Print the results
print(paste("Mean: ", mean_returns))
print(paste("Variance: ", variance_returns))
print(paste("Standard Deviation: ", std_dev_returns))
print(paste("Skewness: ", skewness_returns))
print(paste("Kurtosis: ", kurtosis_returns))

# Normality test ---------------------------------------------------------------
qqnorm(Brent_price$Returns,main="Brent RETURNS")
qqline(Brent_price$Returns)

# Jarque-Bera Test
jarque.bera.test(Brent_returns)

# Ljung-Box Test
Box.test(Brent_returns, type = "Ljung-Box")
Box.test(Brent_returns, lag = 15, type = "Ljung-Box")

# Testing for Stationary -------------------------------------------------------
# Augmented Dickey-Fuller Test
adf.test(Brent_returns, alternative = "stationary")

# Phillips-Perron Test
ur.pp(Brent_returns, type = "Z-alpha")

result <- ur.pp(Brent_returns, type='Z-tau')
cat("Test Statistic:", result@teststat, "\n")
cat("Critical Values (for example): 1%:", result@cval[3], "\n")

# KPSS Test
kpss.test(Brent_returns)

## Heteroskedasticity test -----------------------------------------------------
## ARCH LM test ----------------------------------------------------------------
Brent_returns_ts <- ts(Brent_returns) 

arch.test(Brent_returns_ts, arch = "box", alpha = 0.05, lag.max = 2) # presence of heterocasdicity.

acf_plot <- ggAcf(Brent_returns_ts) + ggtitle("ACF of Brent Returns")
pacf_plot <- ggPacf(Brent_returns_ts) + ggtitle("Pacf of Brent Returns")

grid.arrange(acf_plot, pacf_plot, ncol = 2)

fit <- auto.arima(Brent_returns_ts)
summary(fit)
checkresiduals(fit)

p_brent <- ggplot(Brent_price, aes(x = Date, y = Brent)) +
  geom_line(color = "black") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Decision == 'maintain',], aes(xintercept = as.Date(Date), color = 'Maintain'), linetype = "dotted") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Decision == 'cut',], aes(xintercept = as.Date(Date), color = 'Decrease'), linetype = "dotted") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Decision == 'hike',], aes(xintercept = as.Date(Date), color = 'Increase'), linetype = "dotted") +
  scale_y_continuous(limits = c(0, 150), name = "Spot price/Dollars per Barrel") +
  scale_color_manual(values = c("Maintain" = "gray", "Decrease" = "red", "Increase" = "#4CBB17")) +
  labs(title = "Spot Price of Brent", color = "Decision Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))

##------------------------------------------------------------------------------
## Long memory tests: Table 4 - Full data
## -----------------------------------------------------------------------------
# Calculate the Hurst exponent
hurst_result <- hurstexp(Brent_returns_ts)

## Hurt exponent through regression --------------------------------------------
L <- length(Brent_returns)

# Minimum length of sub series = around 8
vLs <- round(L/2^(0:floor(log2(L/8)))) 

# Output container
df_rs <- data.frame(RS=NA, Ls = vLs)

# Each the length of sub series (Ls)
for(k in 1:length(vLs)) {
  Ls <- vLs[k]
  Ns <- round(L/Ls, 0)
  rs <- rep(NA, Ns)
  
  for(i in 1:Ns) {
    x_is <- Brent_returns[(Ls*(i-1)+1):(Ls*i)]
    x_is <- x_is[!is.na(x_is)]
    ms <- mean(x_is)
    ss <- sd(x_is)
    ds <- x_is - ms
    zs <- cumsum(ds)
    
    rs[i] <- (max(zs) - min(zs))/ss
  }
  
  df_rs[k,] <- c(mean(rs, na.rm = TRUE), Ls)
}

df_rs

#-------------------------------------------------
# Estimation of Hurst exponent

log_fit <- lm(log(RS) ~ log(Ls), data = df_rs)
summary(log_fit)
hurst_coefficient <- log_fit$coefficients[2]

# Print the Hurst exponent
print(hurst_coefficient) # persistent long-term trend 

## Lo's Modified R/S test------------------------------------------------------
# Test statistic (q=1)
L <- length(Brent_returns_ts)
vLs <- round(L/2^(0:floor(log2(L/8)))) 
df_rs <- data.frame(QT=NA, Ls=vLs)

for (k in 1:length(vLs)) {
  Ls <- vLs[k]
  Ns <- round(L/Ls)
  qt <- rep(NA, Ns)
  
  for (i in 1:Ns) {
    x_is <- Brent_returns[(Ls*(i-1)+1):(Ls*i)]
    x_is <- x_is[!is.na(x_is)]
    ms <- mean(x_is)
    T <- length(x_is)
    zs <- cumsum(x_is - ms)
    R_T <- max(zs) - min(zs)
    
    # Use NeweyWest function from sandwich package for HAC estimator
    nw_est <- NeweyWest(lm(x_is ~ 1), lag = 1, prewhite = FALSE)
    S_T_sq <- sum(diag(nw_est))
    
    # Compute Q_T
    qt[i] <- R_T / sqrt(S_T_sq)
  }
  
  df_rs[k,] <- c(mean(qt, na.rm=TRUE), Ls)
}

# Estimation of the Hurst exponent via regression
log_fit <- lm(log(QT) ~ log(Ls), data=df_rs)
summary(log_fit)

# Test statistic (q=5)
for (k in 1:length(vLs)) {
  Ls <- vLs[k]
  Ns <- round(L/Ls)
  qt <- rep(NA, Ns)
  
  for (i in 1:Ns) {
    x_is <- Brent_returns[(Ls*(i-1)+1):(Ls*i)]
    x_is <- x_is[!is.na(x_is)]
    ms <- mean(x_is)
    T <- length(x_is)
    zs <- cumsum(x_is - ms)
    R_T <- max(zs) - min(zs)
    
    # Use NeweyWest function from sandwich package for HAC estimator
    nw_est <- NeweyWest(lm(x_is ~ 1), lag = 5, prewhite = FALSE)
    S_T_sq <- sum(diag(nw_est))
    
    # Compute Q_T
    qt[i] <- R_T / sqrt(S_T_sq)
  }
  
  df_rs[k,] <- c(mean(qt, na.rm=TRUE), Ls)
}

# Estimation of the Hurst exponent via regression
log_fit <- lm(log(QT) ~ log(Ls), data=df_rs)
summary(log_fit)

## GPH Test --------------------------------------------------------------------
gph_result<- fdGPH(Brent_returns_ts, bandw.exp = 0.45)
print(gph_result)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))


## alpha = 0.50 ----------------------------------------------------------------
gph_result <- fdGPH(Brent_returns_ts, bandw.exp = 0.5)
print(gph_result)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))


gph_result <- fdGPH(Brent_returns_ts, bandw.exp = 0.55)
print(gph_result)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))


## GSP Test --------------------------------------------------------------------
# Compute the periodogram
periodogram <- function(log_returns) {
  T <- length(log_returns)
  frequencies <- (2 * pi * (1:(T/2))) / T
  fft_values <- fft(log_returns)[1:(T/2)]
  I_w <- (1 / (2 * pi * T)) * Mod(fft_values)^2
  data.frame(Frequency = frequencies, Periodogram = I_w)
}

# Function R(H)
R_H <- function(H, w, I_w, m) {
  term1 <- log((1/m) * sum(I_w[1:m] / w[1:m]^(2*H)))
  term2 <- (2*H - 1) * (1/m) * sum(log(w[1:m]))
  term1 - term2
}

# Calculate the periodogram from the log returns
log_returns <- na.omit(Brent_returns) # Make sure to remove NA values
peri <- periodogram(log_returns)

# Search for the H that minimizes R(H)
H_values <- seq(0.01, 0.99, length.out = 200)
m <- length(log_returns) / 2
R_values <- sapply(H_values, function(H) R_H(H, peri$Frequency, peri$Periodogram, m))

# Find the H value that minimizes R(H)
min_index <- which.min(R_values)
H_min <- H_values[min_index]

# Print the estimated Hurst exponent
print(H_min)

ggbr_semipara(
  Brent_returns_ts,
  alpha = 0.5,
  method = "gsp",
  min_freq = 0,
  max_freq = 0.5
)

################################################################################
## Long Memory Test - Squared Returns 
################################################################################
Brent_price$SquaredReturns <- Brent_price$Returns^2
Brent_squared_returns <- na.omit(Brent_price$SquaredReturns)
Brent_squared_returns_ts <- ts(Brent_squared_returns)

## Hurst Exponent --------------------------------------------------------------
## Hurt exponent through regression --------------------------------------------
L <- length(Brent_squared_returns)

vLs <- round(L/2^(0:floor(log2(L/8))))

df_rs <- data.frame(RS=NA, Ls = vLs)

for(k in 1:length(vLs)) {
  Ls <- vLs[k]
  Ns <- round(L/Ls, 0)
  rs <- rep(NA, Ns)
  for(i in 1:Ns) {
    x_is <- Brent_squared_returns[(Ls*(i-1)+1):(Ls*i)]
    x_is <- x_is[!is.na(x_is)]
    ms <- mean(x_is)
    ss <- sd(x_is)
    ds <- x_is - ms
    zs <- cumsum(ds)
    
    rs[i] <- (max(zs) - min(zs))/ss
  }
  
  # Save Ls and average of RS pairs
  df_rs[k,] <- c(mean(rs, na.rm = TRUE), Ls)
}

df_rs

log_fit <- lm(log(RS) ~ log(Ls), data = df_rs)
summary(log_fit)

## Lo's modified R/S test ------------------------------------------------------
L <- length(Brent_squared_returns_ts)
vLs <- round(L/2^(0:floor(log2(L/8)))) 
df_rs <- data.frame(QT=NA, Ls=vLs)

for (k in 1:length(vLs)) {
  Ls <- vLs[k]
  Ns <- round(L/Ls)
  qt <- rep(NA, Ns)
  
  for (i in 1:Ns) {
    x_is <- Brent_squared_returns[(Ls*(i-1)+1):(Ls*i)]
    x_is <- x_is[!is.na(x_is)]
    ms <- mean(x_is)
    T <- length(x_is)
    zs <- cumsum(x_is - ms)
    R_T <- max(zs) - min(zs)
    
    # Use NeweyWest function from sandwich package for HAC estimator
    nw_est <- NeweyWest(lm(x_is ~ 1), lag = 5, prewhite = FALSE) # change lag=1  to get the results for q=1
    S_T_sq <- sum(diag(nw_est))
    
    # Compute Q_T
    qt[i] <- R_T / sqrt(S_T_sq)
  }
  
  df_rs[k,] <- c(mean(qt, na.rm=TRUE), Ls)
}

# Estimation of the Hurst exponent via regression
log_fit <- lm(log(QT) ~ log(Ls), data=df_rs)
summary(log_fit)

## GPH test --------------------------------------------------------------------
gph_result<- fdGPH(Brent_squared_returns_ts, bandw.exp = 0.45)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))

## alpha = 0.5 -----------------------------------------------------------------
gph_result <- fdGPH(Brent_squared_returns_ts, bandw.exp = 0.5)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))

## alpha = 0.55 ----------------------------------------------------------------
gph_result <- fdGPH(Brent_squared_returns_ts, bandw.exp = 0.55)
d_hat <- gph_result$d
sd_as <- gph_result$sd.as
t_statistic <- d_hat / sd_as

critical_value <- qnorm(1 - 0.01/2)

is_significant <- abs(t_statistic) > critical_value

print(paste("T-statistic:", t_statistic))
print(paste("Critical value at 1% significance level:", critical_value))
print(paste("Is the estimate statistically significant at the 1% level?:", is_significant))

## Merge datas -----------------------------------------------------------------
Brent_price$Returns <- 100 * c(NA, diff(log(Brent_price$Brent)))
merged_df <- merge(Brent_price, Dollar_index, by="Date", all.x=TRUE)
merged_df <- merge(merged_df, SP500, by="Date", all.x=TRUE)
merged_df <- merge(merged_df, OPEC_meetings, by="Date", all.x=TRUE)
sum(!is.na(merged_df$Decision)) 

merged_df$dxy_r <- 100 * c(NA, diff(log(merged_df$DXY)))
merged_df$SP500_r <- 100 * c(NA, diff(log(merged_df$SP500)))


## Creating dummy variables for OPEC Meetings ----------------------------------
# Create binary columns for each decision type
merged_df <- merged_df %>%
  mutate(cut = as.integer(Decision == "cut"),
         maintain = as.integer(Decision == "maintain"),
         hike = as.integer(Decision == "hike"))

## Handle Missing values in the merged df --------------------------------------
skimr::skim(merged_df)

# Handle NA values in OPEC Cols
cols_to_fill <- c("cut", 
                  "maintain",
                  "hike")

merged_df[cols_to_fill] <- na.fill(merged_df[cols_to_fill], 0)
merged_df$Meeting[is.na(merged_df$Meeting)] <- "No Meeting"
merged_df$Decision[is.na(merged_df$Decision)] <- "No Decision"

# Handle NA value in SP500
merged_df$SP500 <- na_locf(merged_df$SP500)

# Handle NA value in DXY
merged_df$DXY <- na_locf(merged_df$DXY)


# Handdle NA value for Brent returns
merged_df$Returns <- na_locf(merged_df$Returns)

merged_df$dxy_r <- na_locf(merged_df$dxy_r)
merged_df$SP500_r <- na_locf(merged_df$SP500_r)


# Remove ISM, DXY and SP500
merged_df <- merged_df[, !(names(merged_df) %in% c("SP500", "DXY"))]


# Convert your dataframe to an xts object
brent_xts <- xts(merged_df$Returns, order.by = merged_df$Date)

# Extracting the returns data from the xts object as a numeric vector
returns_data <- coredata(brent_xts)

brent_returns <- returns_data  

# Get the ACF values without plotting
acf_values <- acf(brent_returns, lag.max = 30, plot = FALSE)

lag_1_autocorr <- acf_values$acf[2] 
lag_2_autocorr <- acf_values$acf[3] 
lag_3_autocorr <- acf_values$acf[4] 
lag_4_autocorr <- acf_values$acf[5] 
lag_5_autocorr <- acf_values$acf[6] 

print(paste("Autocorrelation at lag 1:", lag_1_autocorr)) # not significant
print(paste("Autocorrelation at lag 2:", lag_2_autocorr)) # not significant
print(paste("Autocorrelation at lag 3:", lag_3_autocorr)) # not significant
print(paste("Autocorrelation at lag 4:", lag_4_autocorr)) # significant
print(paste("Autocorrelation at lag 5:", lag_5_autocorr)) # not significant

# Create lagged returns
merged_df <- merged_df %>%
  mutate(
    rt_lag4 = lag(Returns, 4))

## MODIFYING DUMMY VARIABLE FOR REGRESSION #####################################
modify_all_dummies <- function(data, params) {
  # Helper function to modify a single dummy type
  modify_dummy <- function(data, dummy_col, s1, s2, s3) {
    # Step A1: Shift
    shifted <- dplyr::lag(data[[dummy_col]], n = s1, default = 0)
    
    # Step A2: Extend
    n <- nrow(data)
    extended <- rep(0, n)
    for (i in (1 + s1):n) {
      if (shifted[i] == 1) {
        indices <- seq(i, min(i + s2 - 1, n))
        extended[indices] <- 1
      }
    }
    
    # Step A3: Apply Decay
    adjoining <- rep(0, n)
    for (i in 1:n) {
      if (extended[i] == 1) {
        j <- 1
        decay_value <- s3
        while ((i + j) <= n && decay_value >= 0.1) {
          if (extended[i + j] == 0) {
            adjoining[i + j] <- decay_value
          }
          decay_value <- decay_value * s3
          j <- j + 1
        }
      }
    }
    
    data[[paste0(dummy_col, "_modified")]] <- pmax(extended, adjoining)
    return(data)
  }
  
  data <- modify_dummy(data, "cut", params$s1_cut, params$s2_cut, params$s3_cut)
  data <- modify_dummy(data, "maintain", params$s1_maintain, params$s2_maintain, params$s3_maintain)
  data <- modify_dummy(data, "hike", params$s1_hike, params$s2_hike, params$s3_hike)
  
  return(data)
}

## Choosing the optimal parameters for modified dummy variables ----------------
# Define parameter ranges
s1_range <- 1:5
s2_range <- 1:5
s3_range <- seq(0.1, 0.9, by = 0.1)

# Function to generate a random set of parameters
generate_random_params <- function() {
  list(
    s1_cut = sample(s1_range, 1),
    s2_cut = sample(s2_range, 1),
    s3_cut = sample(s3_range, 1),
    s1_maintain = sample(s1_range, 1),
    s2_maintain = sample(s2_range, 1),
    s3_maintain = sample(s3_range, 1),
    s1_hike = sample(s1_range, 1),
    s2_hike = sample(s2_range, 1),
    s3_hike = sample(s3_range, 1)
  )
}

# Number of random samples
n_samples <- 30000

# Placeholder for the best AIC and parameters
best_aic <- Inf
best_params <- NULL

# Use a loop to sample and evaluate parameters
for (i in 1:n_samples) {
  params <- generate_random_params()
  
  # Modify your data based on the sampled parameters
  modified_data <- modify_all_dummies(merged_df, params)
  
  # Fit the model
  model <- lm(Returns ~ rt_lag4 + cut_modified + maintain_modified + hike_modified, data = modified_data)
  current_aic <- AIC(model)
  
  # Check if this is the best AIC
  if (current_aic < best_aic) {
    best_aic <- current_aic
    best_params <- params
  }
}

# Print the best AIC and parameters
cat("Best AIC:", best_aic, "\nBest Parameters:\n")
print(best_params)

# Optimized parameters
optimized_params <- list(
  s1_cut = 1,
  s2_cut = 2,
  s3_cut = 0.1,
  s1_maintain = 1,
  s2_maintain = 3,
  s3_maintain = 0.5,
  s1_hike = 2,
  s2_hike = 1,
  s3_hike = 0.1
)

# Apply the modify_all_dummies function to your data with the optimized parameters
modified_data <- modify_all_dummies(merged_df, optimized_params)

# Create a new dataframe that contains the dates and modified dummy variables for plotting
plot_data <- modified_data %>%
  select(Date, cut_modified, maintain_modified, hike_modified) %>%
  gather(key = "DecisionType", value = "Value", -Date)

# Plot the modified dummy variables for each decision type
ggplot(plot_data, aes(x = Date, y = Value, color = DecisionType)) +
  geom_line() +
  facet_wrap(~DecisionType, scales = "free_y") +
  theme_minimal() +
  labs(title = "Modified Dummy Variables Over Time",
       x = "Date",
       y = "Value",
       color = "Decision Type")

## Analysing the modified dummy variables --------------------------------------
# create the linear regression model
model <- lm(Returns ~ hike_modified + maintain_modified + cut_modified, data = modified_data)

# Summary of the model to see coefficients and statistics
summary(model)

#  test for autocorrelation in residuals, use the Durbin-Watson test
dwtest(model)

# check for heteroskedasticity, use the Breusch-Pagan test
bptest(model)


## Pure regression models (expected return only, without considering the volatility)
# A. Regression with unmodified dummy variables --------------------------------
model <- lm(Returns ~  rt_lag4 + cut + hike + maintain + SP500_r + dxy_r, data = modified_data)

# Summary of the model to check coefficients and model statistics
summary(model)
AIC(model)

# Assuming 'model' is your fitted lm model object
residuals <- residuals(model)

# Plotting the residuals
par(mfrow=c(2,1)) # Setting up the plotting area to have 2 rows and 1 column

# Plot 1: Residuals
plot(residuals, type = 'l', main = "Residuals of the Regression Model", ylab = "Residuals", xlab = "Observation")
abline(h = 0, col = "red") 

# Plot 2: Squared Residuals
plot(residuals^2, type = 'l', main = "Squared Residuals of the Regression Model", ylab = "Squared Residuals", xlab = "Observation")
abline(h = 0, col = "red") 

# Correlogram of Residuals
acf(residuals, main="ACF of Regression Residuals", lag.max=20)

# Correlogram of Squared Residuals
acf(residuals^2, main="ACF of Squared Regression Residuals", lag.max=20)

par(mfrow=c(1,1))


## B. Regression with optimally modified dummy variables -----------------------
model2 <- lm(Returns ~  rt_lag4 + cut_modified + hike_modified + maintain_modified , data = modified_data)
model2 <- lm(Returns ~  rt_lag4 + cut_modified + hike_modified + maintain_modified +SP500_r + dxy_r , data = modified_data)
# Summary of the model to check coefficients and model statistics
summary(model2)
AIC(model2)

# Calculate residuals from the model
residuals_model <- residuals(model2)

# Plot the ACF of the residuals
acf(residuals_model, main = "ACF of Regression Residuals")

# Plot the ACF of the squared residuals
acf(residuals_model^2, main = "ACF of Squared Regression Residuals")


# Running the regression model
model3 <- lm(Returns ~ rt_lag4 + cut_modified + hike_modified +
               maintain_modified + SP500_r +dxy_r, data = modified_data)

# Summary of the model to check coefficients and model statistics
summary(model3)
AIC(model3)
#everything is significant other than mantain_modified and ism returns.

model4 <- lm(Returns ~ rt_lag4 + cut_modified + hike_modified +
               SP500_r +dxy_r, data = modified_data)
summary(model4)
AIC(model4)
residuals_model <- residuals(model4)


## GARCH MODELS ################################################################

# DETERMINING MEAN EQUATION ----------------------------------------------------
InSampleData  = brent_xts[1:4448]
OutSampleData = brent_xts[4448:5559]

covariates <- merged_df[, c("cut", "maintain", "hike")]
covariates_matrix <- data.matrix(covariates)
covariates_matrix_aligned <- covariates_matrix[1:4448, ]

auto_fit <- auto.arima(InSampleData,xreg=covariates_matrix_aligned, stationary = TRUE, approximation = F, stepwise = TRUE, ic = "aic")
summary(auto_fit)
checkresiduals(auto_fit)

## Arma(2,0) mode with covariates ----------------------------------------------
Arima(InSampleData, xreg=covariates_matrix_aligned, order=c(2,0,0))

## Arma(0,2) mode with covariates ----------------------------------------------
Arima(InSampleData, xreg=covariates_matrix_aligned, order=c(0,0,2))

## Arma(1,2) mode with covariates ----------------------------------------------
Arima(InSampleData, xreg=covariates_matrix_aligned, order=c(1,0,2))

## Arma(2,1) mode with covariates ----------------------------------------------
Arima(InSampleData, xreg=covariates_matrix_aligned, order=c(2,0,1))
# lowest aic

# A) GARCH without dummy variables ---------------------------------------------
## In sample ----------------------
## Normal dist --------
external_regressor <- as.matrix(merged_df[, c("SP500_r")])

## ARFIMA FIGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal Dist---------
norm_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "norm")

# Fit the GARCH(1,1) model using the in-sample data and reserving the last T_test observations for out-of-sample testing
norm_garch_fit <- ugarchfit(spec = norm_spec, data = InSampleData)
show(norm_garch_fit)

## student t Dist---------
std_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                       mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                       distribution.model = "std")

# Fit the GARCH(1,1) model using the in-sample data and reserving the last T_test observations for out-of-sample testing
dtf_garch_fit <- ugarchfit(spec = std_spec, data = InSampleData)
show(dtf_garch_fit)

## skewews student t Dist---------
sstd_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "sstd")

# Fit the GARCH(1,1) model using the in-sample data and reserving the last T_test observations for out-of-sample testing
sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = InSampleData)
show(sstd_garch_fit)

## EGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
norm_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = InSampleData)

show(norm_garch_fit)

## std dist --------
std_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                       mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
                       distribution.model = "std")

std_garch_fit <- ugarchfit(spec = std_spec, data = InSampleData)

show(std_garch_fit)

## sstd dist --------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = InSampleData)

show(sstd_garch_fit)

## Out- Sample -------------------------
T <- nrow(merged_df)-1
T
T_train <- round(0.80*T)
T_train
T_test <- T - T_train
T_test
dates_out_of_sample <- tail(index(merged_df), T_test)
dates_all <- index(merged_df)
dates_in_sample <- dates_all[1:T_train]

## ARFIMA FIGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal Dist---------
norm_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = Brent_xts, out.sample = T_test)
show(norm_garch_fit)

## student t Dist---------
std_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                       mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                       distribution.model = "std")


std_garch_fit <- ugarchfit(spec = std_spec, data = Brent_xts, out.sample = T_test)
show(dtf_garch_fit)

## skewed student t Dist---------
sstd_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = Brent_xts, out.sample = T_test)
show(sstd_garch_fit)

## EGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal dist --------
norm_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "norm")


norm_garch_fit <- ugarchfit(spec = norm_spec, data = brent_xts, out.sample = T_test)
show(norm_garch_fit)

## std dist --------
std_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                       mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
                       distribution.model = "std")

std_garch_fit <- ugarchfit(spec = std_spec, data = brent_xts, out.sample = T_test)

show(std_garch_fit)

## sstd dist --------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                        mean.model = list(armaOrder = c(2, 1), include.mean = TRUE, arfima=T),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = brent_xts, out.sample = T_test)

show(sstd_garch_fit)

## B. GARCH with unmodified dummy variables ####################################
dummies <- as.matrix(merged_df[, c("cut","hike","maintain", "SP500_r")])
dummies <- as.matrix(merged_df[, c("hike","maintain","SP500_r")])

## In-sample --------------------------------
## ARFIMA FIGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

## Normal Dist---------
norm_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = InSampleData)
show(norm_garch_fit)

## student t Dist---------
std_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                       mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                       distribution.model = "std")


std_garch_fit <- ugarchfit(spec = std_spec, data = InSampleData)
show(std_garch_fit)

## skewed student t Dist---------
sstd_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = InSampleData)
show(sstd_garch_fit)

## EGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal dist ----------
norm_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = InSampleData)

show(norm_garch_fit)

## std dist --------
std_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                       mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
                       distribution.model = "std")

std_garch_fit <- ugarchfit(spec = std_spec, data = InSampleData)

show(std_garch_fit)

## sstd dist --------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2,1), include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = InSampleData)

show(sstd_garch_fit)

## Out - Sample ---------------------------------------
## ARFIMA FIGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal Dist---------
norm_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = Brent_xts, out.sample = T_test)
show(norm_garch_fit)

## student t Dist---------
std_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                       mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                       distribution.model = "std")


std_garch_fit <- ugarchfit(spec = std_spec, data = Brent_xts, out.sample = T_test)
show(dtf_garch_fit)

## skewed student t Dist---------
sstd_spec <- ugarchspec(variance.model = list(model = "fiGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1), arfima = T, include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = Brent_xts, out.sample = T_test)
show(sstd_garch_fit)


## EGARCH Model >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Normal dist --------
norm_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "norm")

norm_garch_fit <- ugarchfit(spec = norm_spec, data = Brent_xts, out.sample = T_test)

show(norm_garch_fit)

## std dist --------
std_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                       mean.model = list(armaOrder = c(2, 1), include.mean = TRUE),
                       distribution.model = "std")

std_garch_fit <- ugarchfit(spec = std_spec, data = Brent_xts, out.sample = T_test)

show(std_garch_fit)

## sstd dist --------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2,1), include.mean = TRUE),
                        distribution.model = "sstd")

sstd_garch_fit <- ugarchfit(spec = sstd_spec, data = Brent_xts, out.sample = T_test)

show(sstd_garch_fit)

## C. GARCH with modified dummy variables ######################################
## Annoucement day sc before ----------------------------------------------------
merged_df$announcement <- with(merged_df, as.integer(hike | maintain | cut))
modify_all_dummies_reverse <- function(data, params) {
  modify_dummy_reverse <- function(data, dummy_col, s1, s2, s3) {
    
    shifted <- lead(data[[dummy_col]], n = s1, default = 0)
    
    
    n <- nrow(data)
    extended <- rep(0, n)
    for (i in 1:(n - s1)) {
      if (shifted[i] == 1) {
        indices <- seq(max(1, i - s2 + 1), i)
        extended[indices] <- 1
      }
    }
    
    adjoining <- rep(0, n)
    for (i in n:1) {
      if (extended[i] == 1) {
        j <- 1
        decay_value <- s3
        while ((i - j) >= 1 && decay_value >= 0.1) {
          if (extended[i - j] == 0) {
            adjoining[i - j] <- decay_value
          }
          decay_value <- decay_value * s3
          j <- j + 1
        }
      }
    }
    
    data[[paste0(dummy_col, "_anticipated")]] <- pmax(extended, adjoining)
    return(data)
  }
  
  data$announcement <- as.integer((data$hike | data$maintain | data$cut) > 0)
  
  data <- modify_dummy_reverse(data, "announcement", params$s1_announcement, params$s2_announcement, params$s3_announcement)
  
  return(data)
}

best_params <- list(
  s1_announcement = 0,
  s2_announcement = 3,
  s3_announcement = 0.0
)
# Assuming 'merged_df' is your original data frame and 'modify_all_dummies_reverse' is your function
modified_data <- modify_all_dummies_reverse(merged_df, best_params)

external_regressor <- as.matrix(modified_data[, c("announcement_anticipated")])
spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                   mean.model = list(armaOrder = c(0, 1),include.mean = TRUE, arfima=T),
                   distribution.model = "sstd")

garch_fit <- ugarchfit(spec = spec, data = InSampleData)
garch_fit <- ugarchfit(spec = spec, data = brent_xts, out.sample = T_test)
show(garch_fit)

## Announcement day sc after ---------------------------------------------------
modify_announcement_dummies <- function(data, params) {
  # Helper function to modify the announcement dummy with anticipation and decay effects
  modify_announcement_dummy <- function(data, dummy_col, s1, s2, s3) {
    # Step A1: Shift
    shifted <- dplyr::lag(data[[dummy_col]], n = s1, default = 0)
    
    # Step A2: Extend
    n <- nrow(data)
    extended <- rep(0, n)
    for (i in (1 + s1):n) {
      if (shifted[i] == 1) {
        indices <- seq(i, min(i + s2 - 1, n))
        extended[indices] <- 1
      }
    }
    
    # Step A3: Apply Decay
    adjoining <- rep(0, n)
    for (i in 1:n) {
      if (extended[i] == 1) {
        j <- 1
        decay_value <- s3
        while ((i + j) <= n && decay_value >= 0.1) {
          if (extended[i + j] == 0) {
            adjoining[i + j] <- decay_value
          }
          decay_value <- decay_value * s3
          j <- j + 1
        }
      }
    }
    
    data[[paste0(dummy_col, "_modified")]] <- pmax(extended, adjoining)
    return(data)
  }
  
  # Modify the announcement dummy according to provided parameters
  data <- modify_announcement_dummy(data, "announcement", params$s1_announcement, params$s2_announcement, params$s3_announcement)
  
  return(data)
}

best_params <- list(
  s1_announcement = 0,
  s2_announcement = 3,
  s3_announcement = 0.0
)

# Assuming 'merged_df' is your original data frame and 'modify_all_dummies_reverse' is your function
modified_data <- modify_announcement_dummies(merged_df, best_params)

external_regressor <- as.matrix(modified_data[, c("announcement_modified")])
spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = external_regressor),
                   mean.model = list(armaOrder = c(0, 1),include.mean = TRUE, arfima=T),
                   distribution.model = "sstd")

garch_fit <- ugarchfit(spec = spec, data = InSampleData)
garch_fit <- ugarchfit(spec = spec, data = brent_xts, out.sample = T_test)
show(garch_fit)

## MODIFYING DUMMY VARIABLE FOR GARCH ##########################################
modify_all_dummies_reverse <- function(data, params) {
  # Helper function to modify a single dummy type with anticipation effect
  modify_dummy_reverse <- function(data, dummy_col, s1, s2, s3) {
    # Step A1: Shift forward (anticipation effect)
    shifted <- lead(data[[dummy_col]], n = s1, default = 0)
    
    # Step A2: Extend forward
    n <- nrow(data)
    extended <- rep(0, n)
    for (i in 1:(n - s1)) {
      if (shifted[i] == 1) {
        indices <- seq(max(1, i - s2 + 1), i)
        extended[indices] <- 1
      }
    }
    
    # Step A3: Apply Decay in reverse (anticipation building up to the event)
    adjoining <- rep(0, n)
    for (i in n:1) {
      if (extended[i] == 1) {
        j <- 1
        decay_value <- s3
        while ((i - j) >= 1 && decay_value >= 0.1) {
          if (extended[i - j] == 0) {
            adjoining[i - j] <- decay_value
          }
          decay_value <- decay_value * s3
          j <- j + 1
        }
      }
    }
    
    data[[paste0(dummy_col, "_anticipated")]] <- pmax(extended, adjoining)
    return(data)
  }
  
  # Apply modification function to each dummy type with respective parameters
  data <- modify_dummy_reverse(data, "hike", params$s1_hike, params$s2_hike, params$s3_hike)
  data <- modify_dummy_reverse(data, "maintain", params$s1_maintain, params$s2_maintain, params$s3_maintain)
  return(data)
}


# Function to fit GARCH model and calculate AIC
fit_garch_and_get_aic <- function(data, params) {
  modified_data <- modify_all_dummies_reverse(data, params)
  modified_data <- na.omit(modified_data)  # Removing rows with NAs to avoid errors in model fitting
  
  external_regressors <- as.matrix(modified_data[, c("hike_anticipated", "maintain_anticipated" ,"cut", "SP500_r")])
  spec <- ugarchspec(mean.model = list(armaOrder = c(0,1), include.mean = TRUE, arfima=T),
                     variance.model = list(model = "eGARCH", garchOrder = c(1, 1), external.regressors = external_regressors),
                     distribution.model = "sstd")
  fit <- tryCatch({
    ugarchfit(spec = spec, data = InSampleData, solver = 'hybrid')
  }, error = function(e) {
    NULL  # Return NULL on error
  })
  
  if (!is.null(fit)) {
    aic <- infocriteria(fit)[1]
    return(list(aic = aic, params = params))  # Return both AIC and params
  } else {
    return(list(aic = Inf, params = params))  # Return Inf AIC on error
  }
}


# Setup for parallel computing
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

# Example parameter ranges and sample size for the demonstration
s1_range <- 0:3
s2_range <- 1:3
s3_range <- seq(0.1, 0.9, by = 0.1)
n_samples <- 300 # Number of random samples for testing

# Generating random parameter combinations
generate_random_params <- function(n) {
  lapply(1:n, function(x) {
    list(
      s1_maintain = sample(s1_range, 1),
      s2_maintain = sample(s2_range, 1),
      s3_maintain = sample(s3_range, 1),
      s1_hike = sample(s1_range, 1),
      s2_hike = sample(s2_range, 1),
      s3_hike = sample(s3_range, 1)
    )
  })
}

random_params <- generate_random_params(n_samples)

# Export necessary objects to the cluster
clusterExport(cl, varlist = c("modify_all_dummies_reverse", "fit_garch_and_get_aic", "random_params", "merged_df"))

# Performing optimization
results <- foreach(param = random_params, .packages = c("rugarch", "dplyr")) %dopar% {
  fit_garch_and_get_aic(merged_df, param)
}


# Stop the parallel cluster
stopCluster(cl)

# Extracting and processing results
aic_values <- sapply(results, `[[`, "aic")
min_aic_index <- which.min(aic_values)
min_aic <- aic_values[min_aic_index]
best_params <- results[[min_aic_index]]$params

# Output the best AIC and the corresponding parameters
cat("Best AIC:", min_aic, "\nBest Parameters:\n")
print(best_params)

optimized_params <- list(
  s1_maintain = 0,
  s2_maintain = 1,
  s3_maintain = 0.4,
  s1_hike = 3,
  s2_hike = 1,
  s3_hike = 0.8
)


modified_data <- modify_all_dummies_reverse(merged_df, optimized_params)

dummies <- as.matrix(modified_data[,c("hike_anticipated","maintain_anticipated", "cut","SP500_r")])

## sstd ----------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1,1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "sstd")

fit <- ugarchfit(spec = sstd_spec, data =InSampleData, solver = 'hybrid')
fit <- ugarchfit(spec = sstd_spec, data = brent_xts, out.sample = T_test, solver = 'hybrid')

show(fit)
# Interpretation and Diagnostics -----------
err <- fit@fit[["residuals"]] 
# Diagnostic plots
par(mfrow = c(2, 3))
plot(err, main = "Residuals")
hist(err, breaks = 50, main = "Histogram of Residuals")
qqnorm(err); qqline(err, col = "red", main = "Q-Q Plot of Residuals")
acf(err, main = "ACF of Residuals")
pacf(err, main = "PACF of Residuals")

# Statistical tests
Box.test(err, lag = 16, type = "Ljung-Box")
adf.test(err)
jarque.bera.test(err)

# squared residuals
squared_err <- err^2
squared_err_ts <- ts(squared_err)

plot(squared_err, main = "Squared Residuals")
hist(squared_err, breaks = 50, main = "Histogram of Squared Residuals")
acf(squared_err, main = "ACF of Squared Residuals")
pacf(squared_err, main = "PACF of Squared Residuals")

Box.test(squared_err, lag = 20, type = "Ljung-Box")
arch.test(squared_err_ts, arch = "box", alpha = 0.05, lag.max = 2)
adf.test(squared_err)
jarque.bera.test(squared_err)
## Producing plots  ----------------
quantile(brent_xts , 0.05)
qplot(brent_xts , geom = 'histogram') + geom_histogram(fill = 'lightblue' , bins = 30) +
  geom_histogram(aes(brent_xts[brent_xts < quantile(brent_xts , 0.05)]) , fill = 'red' , bins = 30) +
  labs(x = 'Daily Returns')


fitdist(distribution = 'sstd' , x =brent_xts)$pars

cat("For a = 0.05 the quantile value of skewed t-distribution is: " ,
    qdist(distribution = 'sstd' , shape = 3.56514510 , p = 0.05) , "\n" , "\n" , 
    "For a = 0.01 the quantile value of skewed t-distribution is: " , 
    qdist(distribution = 'sstd' , shape = 3.56514510 , p = 0.01) , sep = "")

qplot(y = brent_xts , x = 1:5559 , geom = 'point') + geom_point(colour = 'lightgrey' , size = 2) + 
  geom_line(aes(y = fit@fit$sigma*(-1.463507) , x = 1:5559) , colour = 'red') +
  geom_hline(yintercept = sd(brent_xts)*qnorm(0.05) , colour = 'darkgreen' , size = 1.2) + theme_light() + 
  labs(x = '' , y = 'Daily Returns' , title = 'Value at Risk Comparison')

control <- list(maxiter = 10000, trace = TRUE, algorithm = "hybrid", reltol = 1e-8)

# Retry ugarchroll with the same specification but adjusted control settings
model.roll <- ugarchroll(spec = sstd_spec, data = brent_xts[1:4448], n.start = 3948, refit.window = "moving",
                         calculate.VaR = TRUE, VaR.alpha = c(0.05, 0.01))

# For 99% confidence level
report(model.roll, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)

# For 95% confidence level
report(model.roll, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)

roll = ugarchroll(spec = sstd_spec , data = brent_xts[1:4448] , n.start = 3948,
                  refit.window = 'moving',calculate.VaR = TRUE,
                  VaR.alpha = c(0.05, 0.01))

plot(roll, which = 4)

## VAR Backtesting #############################################################
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(0, 1), include.mean = TRUE, arfima=T),
                        distribution.model = "sstd")


model.roll <- ugarchroll(sstd_spec, data = brent_xts[1:4448], refit.window = "moving",
                         calculate.VaR = TRUE, VaR.alpha = c(0.05, 0.01))

# For 99% confidence level
report(model.roll, type = "VaR", VaR.alpha = 0.01, conf.level = 0.95)

# For 95% confidence level
report(roll, type = "VaR", VaR.alpha = 0.05, conf.level = 0.95)

# Extracting VaR and actual values
if ("forecast" %in% slotNames(model.roll) && !is.null(model.roll@forecast$VaR)) {
  EbacktestVaR_5 <- zoo(model.roll@forecast$VaR[, "alpha(5%)"], order.by = as.Date(rownames(model.roll@forecast$VaR)))
  EbacktestVaR_1 <- zoo(model.roll@forecast$VaR[, "alpha(1%)"], order.by = as.Date(rownames(model.roll@forecast$VaR)))
  Ebackactual <- zoo(model.roll@forecast$VaR[, "realized"], order.by = as.Date(rownames(model.roll@forecast$VaR)))
  
  # Plotting the backtest results
  plot(Ebackactual, type = "b", main = " Brent VaR Backtesting for 5% and 1%", xlab = "Date", ylab = "Return/VaR in percent",
       col = "black", lwd = 2, xaxt = "n")  # disable default x-axis to customize
  lines(EbacktestVaR_5, col = "red", lwd = 2, lty = "dashed")
  lines(EbacktestVaR_1, col = "blue", lwd = 2, lty = "dashed")
  
  # Color points where actual returns exceed VaR thresholds
  points(index(Ebackactual[Ebackactual < EbacktestVaR_5]), Ebackactual[Ebackactual < EbacktestVaR_5], col = "red", pch = 19, cex = 0.8)
  points(index(Ebackactual[Ebackactual < EbacktestVaR_1]), Ebackactual[Ebackactual < EbacktestVaR_1], col = "blue", pch = 19, cex = 0.8)
  
  # Customizing the x-axis with formatted dates
  axis.Date(1, at = seq(min(index(Ebackactual)), max(index(Ebackactual)), by = "months"),
            labels = format(seq(min(index(Ebackactual)), max(index(Ebackactual)), by = "months"), "%b %Y"),
            cex.axis = 0.7, las = 2)  # 'las = 2' will make labels perpendicular to the axis
  
  legend("topleft", legend = c("Actual Returns", "VaR 5%", "VaR 1%", "Breaches 5%", "Breaches 1%"),
         col = c("black", "red", "blue", "red", "blue"), lty = c(1, 2, 2, 0, 0), pch = c(NA, NA, NA, 19, 19), 
         lwd = 2, bty = "n", cex = 0.75)  # Adjust legend text size with cex, position to top left
} else {
  cat("No VaR forecast data available. Check model output.\n")
}

## Expected Shortfall ----------------------------------------------------------
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1, 1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "sstd")


fit = ugarchfit(sstd_spec, data = brent_xts[1:4448])
spec2 = sstd_spec
setfixed(spec2)<-as.list(coef(fit))
filt = ugarchfilter(spec2, brent_xts[4449:5559])
actual = brent_xts[4449:5559]

# 95% quantile
VaR = fitted(filt) + sigma(filt)*qdist("sstd", p=0.05, mu = 0, sigma = 1, 
                                       skew  = coef(fit)["skew"], shape=coef(fit)["shape"])
# calculate ES
f = function(x) qdist("sstd", p=x, mu = 0, sigma = 1, 
                      skew  = coef(fit)["skew"], shape=coef(fit)["shape"])
ES = fitted(filt) + sigma(filt)*integrate(f, 0, 0.05)$value/0.05
print(ESTest(0.05, actual, ES, VaR, boot = TRUE))

# 99% quantile ------------------------
# Calculate the VaR for the 99% quantile
VaR_99 = fitted(filt) + sigma(filt) * qdist("sstd", p=0.01, mu = 0, sigma = 1, 
                                            skew  = coef(fit)["skew"], shape=coef(fit)["shape"])

f_99 = function(x) qdist("sstd", p=x, mu = 0, sigma = 1, 
                         skew  = coef(fit)["skew"], shape=coef(fit)["shape"])
ES_99 = fitted(filt) + sigma(filt) * integrate(f_99, 0, 0.01)$value / 0.01

print(ESTest(0.01, actual, ES_99, VaR_99, boot = TRUE))

## Forecasting -----------------------------------------------------------------
## Forecast for dummy model and modified dummy model and compare.
sstd_spec <- ugarchspec(variance.model = list(model = "eGARCH", garchOrder = c(1,1),external.regressors = dummies),
                        mean.model = list(armaOrder = c(2, 1),include.mean = TRUE),
                        distribution.model = "sstd")

fit <- ugarchfit(spec = sstd_spec, data = brent_xts, out.sample = T_test, solver = 'hybrid')

forecasts <- ugarchforecast(fit, n.ahead = 500, n.roll = T_test-1)

# Extract the forecasted values
forecasted_values <- forecasts@forecast$seriesFor[1, ]

forecast_dates <- index(brent_xts)[4449:5558]

# Create the xts object for forecasted log returns
forecast_log_returns <- xts(forecasted_values, order.by = forecast_dates)

volatility <-forecasts@forecast$sigmaFor[1, ]
forecast_volatility <- xts(volatility, order.by = forecast_dates)
fit4<-fitted(fit)
fit4[1:1500]
fit4[1501:2500]
fit4[2502:3000]

actual_log_returns <- brent_xts[forecast_dates]

accuracy_forecast <- accuracy(as.numeric(coredata(forecast_log_returns)), as.numeric(coredata(actual_log_returns)))

plot(forecasts)

ggplot(Brent_price, aes(x = Date, y = Brent)) +
  geom_line(color = "blue") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Meeting == 'Extraordinary OPEC Conference' & OPEC_meetings$Decision == 'maintain',], aes(xintercept = as.Date(Date), color = 'Maintain'), linetype = "dotted") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Meeting == 'Extraordinary OPEC Conference' & OPEC_meetings$Decision == 'cut',], aes(xintercept = as.Date(Date), color = 'Decrease'), linetype = "dotted") +
  geom_vline(data = OPEC_meetings[OPEC_meetings$Meeting == 'Extraordinary OPEC Conference' & OPEC_meetings$Decision == 'hike',], aes(xintercept = as.Date(Date), color = 'Increase'), linetype = "dotted") +
  scale_y_continuous(limits = c(0, 150), name = "Spot price/Dollars per Barrel") +
  scale_color_manual(values = c("Maintain" = "darkgray", "Decrease" = "red", "Increase" = "#4CBB17")) +
  labs(title = "Spot Price of Brent - Extraordinary OPEC Conference Decisions", color = "Decision Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5))




