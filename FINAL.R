# Author: Cooper Lawrenz (With AI ofc)

# Load libraries
suppressMessages({
  library(dplyr); library(ggplot2); library(lubridate); library(xtable)
  library(urca); library(tseries); library(stargazer); library(scales)
  library(sandwich); library(lmtest); library(car); library(reshape2); library(tidyr)
})

# Data preparation
dta <- read.csv("/Users/cooperlawrenz/Documents/LSE Work/Dissertation/Index/Code/Data/Full.csv") #if reproducing, you will have to change this
dta$Date <- as.Date(paste(dta$Year, dta$Month, dta$Day, sep="-"), format="%Y-%m-%d")
dta$Yield <- (3 / dta$Consols) * 100

analysis_data <- dta %>%
  filter(Date >= as.Date("1753-08-01") & Date <= as.Date("1844-11-30")) %>%
  arrange(Date) %>%
  mutate(Market_debt_growth = c(NA, diff(log(Market_debt)) * 100),
         dYield = c(NA, diff(Yield)),
         Month_factor = as.factor(Month)) %>%
  slice(-1) %>%
  select(Date, dYield, Market_debt_growth, Shock, `Gold.Standard`, 
         `Partly.Expected.War`, Consols, Yield, Month_factor) %>%
  filter(complete.cases(.))

print(paste("Data:", nrow(analysis_data), "observations"))

# ADF Unit Root Tests
quick_adf <- function(data, variables) {
  results <- data.frame(Variable = character(), Level = character(), 
                        Test_Stat = numeric(), P_Value = character(), 
                        Decision = character(), stringsAsFactors = FALSE)
  
  for (var in variables) {
    series <- data[[var]][!is.na(data[[var]])]
    
    # Test levels and differences
    for (type in c("Levels", "First Diff")) {
      test_series <- if(type == "Levels") series else diff(series)
      if(length(test_series) > 2) {
        tryCatch({
          adf_test <- ur.df(test_series, type = "drift", lags = 1)
          adf_sum <- summary(adf_test)
          results <- rbind(results, data.frame(
            Variable = var, Level = type,
            Test_Stat = round(adf_sum@teststat[1], 3),
            P_Value = ifelse(adf_sum@teststat[1] < adf_sum@cval[1,2], "< 0.05", "> 0.10"),
            Decision = ifelse(adf_sum@teststat[1] < adf_sum@cval[1,2], "Stationary", "Unit Root")
          ))
        }, error = function(e) cat("Error with", var, type, "\n"))
      }
    }
  }
  return(results)
}

adf_results <- quick_adf(analysis_data, c("Yield", "dYield", "Market_debt_growth", "Shock"))
print(xtable(adf_results, caption = "ADF Unit Root Test Results", label = "tab:adf_tests"), 
      type = "latex", caption.placement = "top", include.rownames = FALSE)

# Optimal Lags and Granger Causality
find_optimal_lags <- function(data, vars, max_lags = 6) {
  results <- data.frame(From = character(), To = character(), 
                        AIC_Lags = numeric(), BIC_Lags = numeric(), 
                        HIC_Lags = numeric(), stringsAsFactors = FALSE)
  
  clean_data <- data[, vars, drop = FALSE][complete.cases(data[, vars, drop = FALSE]), ]
  
  for (from_var in vars) {
    for (to_var in vars[vars != from_var]) {
      best_lags <- c(AIC = 1, BIC = 1, HIC = 1)
      best_vals <- c(AIC = Inf, BIC = Inf, HIC = Inf)
      
      for (p in 1:max_lags) {
        temp_df <- clean_data
        for (lag_i in 1:p) {
          temp_df[paste0(from_var, "_lag", lag_i)] <- c(rep(NA, lag_i), 
                                                        temp_df[[from_var]][1:(nrow(temp_df)-lag_i)])
          temp_df[paste0(to_var, "_lag", lag_i)] <- c(rep(NA, lag_i), 
                                                      temp_df[[to_var]][1:(nrow(temp_df)-lag_i)])
        }
        
        temp_df <- temp_df[(p+1):nrow(temp_df), ][complete.cases(temp_df[(p+1):nrow(temp_df), ]), ]
        if (nrow(temp_df) < 30) next
        
        formula_str <- paste0(to_var, " ~ ", 
                              paste0(to_var, "_lag", 1:p, collapse = " + "), " + ",
                              paste0(from_var, "_lag", 1:p, collapse = " + "))
        
        tryCatch({
          model <- lm(as.formula(formula_str), data = temp_df)
          aic_val <- AIC(model)
          bic_val <- BIC(model)
          hic_val <- -2 * logLik(model) + 2 * length(coef(model)) * log(log(nrow(temp_df)))
          
          if (aic_val < best_vals["AIC"]) { best_vals["AIC"] <- aic_val; best_lags["AIC"] <- p }
          if (bic_val < best_vals["BIC"]) { best_vals["BIC"] <- bic_val; best_lags["BIC"] <- p }
          if (hic_val < best_vals["HIC"]) { best_vals["HIC"] <- hic_val; best_lags["HIC"] <- p }
        }, error = function(e) NULL)
      }
      
      results <- rbind(results, data.frame(From = from_var, To = to_var,
                                           AIC_Lags = best_lags["AIC"], 
                                           BIC_Lags = best_lags["BIC"], 
                                           HIC_Lags = best_lags["HIC"]))
    }
  }
  return(results)
}

granger_causality <- function(data, vars, lag_results) {
  results <- data.frame(From = character(), To = character(), Lags_Used = numeric(),
                        F_Stat = numeric(), P_Value = numeric(), Significant = character(),
                        stringsAsFactors = FALSE)
  
  clean_data <- data[, vars, drop = FALSE][complete.cases(data[, vars, drop = FALSE]), ]
  
  for (i in 1:nrow(lag_results)) {
    from_var <- lag_results$From[i]
    to_var <- lag_results$To[i]
    optimal_p <- lag_results$AIC_Lags[i]
    
    temp_df <- clean_data
    for (lag_j in 1:optimal_p) {
      temp_df[paste0(from_var, "_lag", lag_j)] <- c(rep(NA, lag_j), 
                                                    temp_df[[from_var]][1:(nrow(temp_df)-lag_j)])
      temp_df[paste0(to_var, "_lag", lag_j)] <- c(rep(NA, lag_j), 
                                                  temp_df[[to_var]][1:(nrow(temp_df)-lag_j)])
    }
    
    temp_df <- temp_df[(optimal_p+1):nrow(temp_df), ][complete.cases(temp_df[(optimal_p+1):nrow(temp_df), ]), ]
    if (nrow(temp_df) < 30) next
    
    restricted_formula <- paste0(to_var, " ~ ", paste0(to_var, "_lag", 1:optimal_p, collapse = " + "))
    unrestricted_formula <- paste0(restricted_formula, " + ", paste0(from_var, "_lag", 1:optimal_p, collapse = " + "))
    
    tryCatch({
      restricted_model <- lm(as.formula(restricted_formula), data = temp_df)
      unrestricted_model <- lm(as.formula(unrestricted_formula), data = temp_df)
      f_test <- anova(restricted_model, unrestricted_model)
      
      results <- rbind(results, data.frame(
        From = from_var, To = to_var, Lags_Used = optimal_p,
        F_Stat = round(f_test$F[2], 3), P_Value = round(f_test$`Pr(>F)`[2], 3),
        Significant = ifelse(f_test$`Pr(>F)`[2] < 0.05, "Yes", "No")))
    }, error = function(e) NULL)
  }
  return(results)
}

# Run Granger analysis
key_vars <- c("dYield", "Market_debt_growth", "Shock", "Partly.Expected.War")
optimal_lags_combined <- find_optimal_lags(analysis_data, key_vars)
granger_results <- granger_causality(analysis_data, key_vars, optimal_lags_combined)

# Create Granger matrix
granger_matrix <- matrix("-", length(key_vars), length(key_vars), 
                         dimnames = list(key_vars, key_vars))
diag(granger_matrix) <- "-"

for (i in 1:nrow(granger_results)) {
  from_idx <- which(key_vars == granger_results$From[i])
  to_idx <- which(key_vars == granger_results$To[i])
  granger_matrix[from_idx, to_idx] <- paste0(granger_results$P_Value[i],
                                             ifelse(granger_results$Significant[i] == "Yes", "*", ""))
}

# Print tables
invisible(lapply(list(
  list(optimal_lags_combined, "Optimal Lag Selection Results", "tab:optimal_lags"),
  list(granger_matrix, "Granger Causality P-Values (* p < 0.05)", "tab:granger_matrix"),
  list(granger_results, "Detailed Granger Causality Results", "tab:granger_detailed")
), function(x) print(xtable(x[[1]], caption = x[[2]], label = x[[3]]), 
                     type = "latex", caption.placement = "top", include.rownames = FALSE)))

# Summary Statistics
summary_vars <- c("dYield", "Market_debt_growth", "Shock", "Partly.Expected.War")
summary_stats <- data.frame(
  Variable = summary_vars,
  Mean = round(sapply(summary_vars, function(x) mean(analysis_data[[x]], na.rm=T)), 3),
  SD = round(sapply(summary_vars, function(x) sd(analysis_data[[x]], na.rm=T)), 3),
  Min = round(sapply(summary_vars, function(x) min(analysis_data[[x]], na.rm=T)), 3),
  Max = round(sapply(summary_vars, function(x) max(analysis_data[[x]], na.rm=T)), 3),
  Obs = sapply(summary_vars, function(x) sum(!is.na(analysis_data[[x]])))
)

print(xtable(summary_stats, caption = "Summary Statistics", label = "tab:summary"), 
      type = "latex", caption.placement = "top", include.rownames = FALSE)

# Seasonality test
data_vector <- na.omit(analysis_data$Market_debt_growth)
seasons <- rep(1:12, length.out = length(data_vector))
print(round(kruskal.test(data_vector, factor(seasons))$p.value, 4))

# Traditional Linear Regression
regression_data <- analysis_data %>%
  mutate(dYield_L1 = lag(dYield, 1),
         Market_debt_growth_L1 = lag(Market_debt_growth, 1),
         Shock_L1 = lag(Shock, 1),
         Partly_Expected_War_L1 = lag(`Partly.Expected.War`, 1)) %>%
  slice(-1) %>% filter(complete.cases(.))

models <- list(
  lm(dYield ~ Market_debt_growth + Market_debt_growth_L1, data = regression_data),
  lm(dYield ~ Market_debt_growth + Market_debt_growth_L1 + dYield_L1, data = regression_data),
  lm(dYield ~ Market_debt_growth + Market_debt_growth_L1 + dYield_L1 + Month_factor, data = regression_data),
  lm(dYield ~ Market_debt_growth + Market_debt_growth_L1 + dYield_L1 + Shock + Shock_L1 + 
       `Partly.Expected.War` + Partly_Expected_War_L1 + Month_factor, data = regression_data),
  lm(dYield ~ Market_debt_growth + Market_debt_growth_L1 + dYield_L1 + I(Market_debt_growth * `Gold.Standard`) + 
       `Gold.Standard` + Shock + Shock_L1 + `Partly.Expected.War` + Partly_Expected_War_L1 + 
       Month_factor, data = regression_data)
)

stargazer(models, type = "latex", title = "Debt Growth Impact on Bond Yield Changes",
          label = "tab:debt_yield_models", column.labels = c("Simple", "AR", "Time", "Shocks", "Full"),
          omit = "Month", omit.labels = "Month Fixed Effects", digits = 4, header = FALSE)

# Historical Yields Plot
plot_data <- dta %>%
  filter(Date >= as.Date("1753-08-01") & Date <= as.Date("1844-11-30")) %>%
  arrange(Date) %>%
  mutate(Yield = (3 / Consols) * 100) %>%
  select(Date, Yield, Gold.Standard = `Gold.Standard`, War = `Partly.Expected.War`) %>%
  filter(complete.cases(.))

war_rectangles <- data.frame(
  xmin = as.Date(c("1756-05-01", "1775-04-01", "1792-04-01", "1803-05-01", "1812-06-01")),
  xmax = as.Date(c("1763-02-01", "1783-09-01", "1802-03-01", "1814-04-01", "1815-06-01")),
  war_label = c("Seven Years' War", "American Revolutionary War", "French Revolutionary Wars", 
                "Napoleonic Wars", "War of 1812")
) %>% filter(xmax >= min(plot_data$Date) & xmin <= max(plot_data$Date)) %>%
  mutate(xmin = pmax(xmin, min(plot_data$Date)), xmax = pmin(xmax, max(plot_data$Date)))

yield_plot <- ggplot(plot_data, aes(x = Date, y = Yield)) +
  geom_rect(data = war_rectangles, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, 
                                       fill = war_label), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(color = "black", linewidth = 0.8) +
  scale_fill_manual(name = "War Periods", 
                    values = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "inside", legend.position.inside = c(0.02, 0.95),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  scale_x_date(date_breaks = "10 years", date_labels = "%Y") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(title = "British Government Bond Yields, 1753–1844", 
       subtitle = "3% Consols during five major conflicts", x = "Year", y = "Yield")

print(yield_plot)
ggsave("historical_yields_enhanced.png", yield_plot, width = 12, height = 7, dpi = 300)

# Local Projections Function (handles both standard and NW errors)
run_lp <- function(data, shock_var, gold_standard = FALSE, horizon = 12, newey_west = FALSE) {
  temp_data <- data %>%
    mutate(dYield_L1 = lag(dYield, 1), shock_L1 = lag(.data[[shock_var]], 1),
           Market_debt_growth_L1 = lag(Market_debt_growth, 1),
           Market_debt_growth_L2 = lag(Market_debt_growth, 2),
           Market_debt_growth_L3 = lag(Market_debt_growth, 3),
           Market_debt_growth_L4 = lag(Market_debt_growth, 4))
  
  results <- list()
  lag_terms <- "dYield_L1 + shock_L1 + Market_debt_growth_L1 + Market_debt_growth_L2 + Market_debt_growth_L3 + Market_debt_growth_L4"
  
  for (h in 0:horizon) {
    temp_data$y_lead <- lead(temp_data$dYield, h)
    clean_data <- temp_data %>% slice(5:(n()-h)) %>% filter(complete.cases(.))
    if (nrow(clean_data) < 20) next
    
    formula_str <- if (gold_standard) {
      paste("y_lead ~", shock_var, "+ Gold.Standard + I(", shock_var, 
            " * Gold.Standard) + Market_debt_growth +", lag_terms, "+ Month_factor")
    } else {
      paste("y_lead ~", shock_var, "+ Market_debt_growth +", lag_terms, "+ Month_factor")
    }
    
    model <- lm(as.formula(formula_str), data = clean_data)
    coef_val <- coef(model)[shock_var]
    
    # Calculate standard errors
    if (newey_west) {
      T_obs <- nrow(clean_data)
      lag_trunc <- max(1, floor(4 * (T_obs/100)^(2/9)))
      vcov_nw <- NeweyWest(model, lag = lag_trunc)
      shock_idx <- which(names(coef(model)) == shock_var)
      se_val <- sqrt(vcov_nw[shock_idx, shock_idx])
    } else {
      se_val <- summary(model)$coefficients[shock_var, "Std. Error"]
    }
    
    results[[h+1]] <- data.frame(horizon = h, coefficient = coef_val,
                                 lower_ci = coef_val - 1.96 * se_val,
                                 upper_ci = coef_val + 1.96 * se_val)
  }
  return(do.call(rbind, results))
}

# Run all LP specifications
lp_specs <- expand.grid(
  shock = c("Shock", "Partly.Expected.War"),
  gold = c(FALSE, TRUE),
  nw = c(FALSE, TRUE),
  stringsAsFactors = FALSE
)

lp_results <- list()
for (i in 1:nrow(lp_specs)) {
  spec_name <- paste0(
    ifelse(lp_specs$shock[i] == "Shock", "All Shocks", "War Shocks"),
    ifelse(lp_specs$gold[i], " + Gold", ""),
    ifelse(lp_specs$nw[i], " (NW)", "")
  )
  lp_results[[spec_name]] <- run_lp(analysis_data, lp_specs$shock[i], 
                                    lp_specs$gold[i], newey_west = lp_specs$nw[i])
  lp_results[[spec_name]]$model <- spec_name
}

all_results <- do.call(rbind, lp_results)
all_results$model <- factor(all_results$model)

# Create plots for both standard and NW errors
create_lp_plot <- function(results_data, title_suffix = "", color = "darkblue", fill = "lightblue") {
  ggplot(results_data, aes(x = horizon, y = coefficient)) +
    geom_hline(yintercept = 0, linetype = "solid", color = "black", size = 0.5) +
    geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = fill, alpha = 0.3) +
    geom_line(color = color, size = 0.8) +
    facet_wrap(~ model, ncol = 2, scales = "free_y") +
    labs(x = "Horizon (months)", y = "Impact on Yield Change (basis points)",
         title = paste("Local Projections: Impact of Shocks on Government Bond Yields", title_suffix)) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          strip.text = element_text(size = 10, face = "bold")) +
    scale_x_continuous(breaks = seq(0, 12, 2)) + coord_cartesian(xlim = c(0, 12))
}

# Standard errors plot
std_results <- all_results[!grepl("\\(NW\\)", all_results$model), ]
main_plot <- create_lp_plot(std_results)
print(main_plot)
ggsave("main_plot.png", main_plot, width = 12, height = 9, dpi = 300)

# Newey-West errors plot
nw_results <- all_results[grepl("\\(NW\\)", all_results$model), ]
main_plot_nw <- create_lp_plot(nw_results, "(NW Errors)", "darkred", "lightcoral")
print(main_plot_nw)
ggsave("main_plot_nw.png", main_plot_nw, width = 12, height = 9, dpi = 300)

# Robustness Analysis (Condensed)
run_lp_robust <- function(data, shock_var, outcome_var = "dYield", lag_structure = "medium", 
                          controls = "standard", sample_period = "full", gold_standard = FALSE, 
                          horizon = 12) {
  # Apply sample filters
  if (sample_period == "early") data <- data %>% filter(Date <= as.Date("1800-01-01"))
  else if (sample_period == "late") data <- data %>% filter(Date >= as.Date("1800-01-01"))
  else if (sample_period == "no_war") data <- data %>% filter(`Partly.Expected.War` == 0)
  
  # Define lag variables based on structure
  lag_vars <- switch(lag_structure,
                     "short" = c("dYield_L1"),
                     "medium" = c("dYield_L1", "shock_L1", "Market_debt_growth_L1"),
                     "long" = c("dYield_L1", "dYield_L2", "shock_L1", "shock_L2", "Market_debt_growth_L1", "Market_debt_growth_L2")
  )
  
  temp_data <- data %>%
    mutate(dYield_L1 = lag(dYield, 1), dYield_L2 = lag(dYield, 2),
           shock_L1 = lag(.data[[shock_var]], 1), shock_L2 = lag(.data[[shock_var]], 2),
           Market_debt_growth_L1 = lag(Market_debt_growth, 1),
           Market_debt_growth_L2 = lag(Market_debt_growth, 2))
  
  results <- list()
  for (h in 0:horizon) {
    temp_data$y_lead <- lead(temp_data[[outcome_var]], h)
    clean_data <- temp_data %>% slice(3:(n()-h)) %>% filter(complete.cases(.))
    if (nrow(clean_data) < 20) next
    
    base_formula <- paste("y_lead ~", shock_var, "+ Market_debt_growth")
    if (length(lag_vars) > 0) base_formula <- paste(base_formula, "+", paste(lag_vars, collapse = " + "))
    if (controls == "extended") base_formula <- paste(base_formula, "+ Month_factor")
    if (gold_standard) base_formula <- paste(base_formula, "+ Gold.Standard + I(", shock_var, " * Gold.Standard)")
    
    model <- lm(as.formula(base_formula), data = clean_data)
    coef_val <- coef(model)[shock_var]
    se_val <- summary(model)$coefficients[shock_var, "Std. Error"]
    
    results[[h+1]] <- data.frame(horizon = h, coefficient = coef_val, se = se_val,
                                 lower_ci = coef_val - 1.96 * se_val,
                                 upper_ci = coef_val + 1.96 * se_val, obs = nrow(clean_data))
  }
  return(do.call(rbind, results))
}

# Run complete robustness analysis
shock_vars <- c("Shock", "Partly.Expected.War")
enhanced_data <- analysis_data %>%
  mutate(Shock_abs = abs(Shock), War_abs = abs(`Partly.Expected.War`),
         Shock_std = scale(Shock)[,1], War_std = scale(`Partly.Expected.War`)[,1])

robustness_specs <- list(
  baseline = list(lag_structure = "medium", controls = "standard", sample_period = "full", gold_standard = FALSE),
  short_lags = list(lag_structure = "short", controls = "standard", sample_period = "full", gold_standard = FALSE),
  long_lags = list(lag_structure = "long", controls = "standard", sample_period = "full", gold_standard = FALSE),
  extended_controls = list(lag_structure = "medium", controls = "extended", sample_period = "full", gold_standard = FALSE),
  early_period = list(lag_structure = "medium", controls = "standard", sample_period = "early", gold_standard = FALSE),
  late_period = list(lag_structure = "medium", controls = "standard", sample_period = "late", gold_standard = FALSE),
  gold_interaction = list(lag_structure = "medium", controls = "standard", sample_period = "full", gold_standard = TRUE)
)

all_robustness <- list()
for (shock_var in shock_vars) {
  cat("Running robustness for", shock_var, "\n")
  robustness_results <- list()
  
  for (spec_name in names(robustness_specs)) {
    spec <- robustness_specs[[spec_name]]
    tryCatch({
      robustness_results[[spec_name]] <- do.call(run_lp_robust, 
                                                 c(list(data = enhanced_data, shock_var = shock_var), spec))
    }, error = function(e) NULL)
  }
  
  all_robustness[[shock_var]] <- robustness_results
}

