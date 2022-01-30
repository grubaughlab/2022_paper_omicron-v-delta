#############################################################################################################################
# About ----------------------------------------------
#############################################################################################################################

# Code for Figure 1
## 1A - Delta and Omicron total (SGTF + NGS) daily cases + % omicron over time
## 1B - Delta and Omicron logistic growth curves over emergence periods + doubling times in days

#############################################################################################################################
# Libraries & Data ----------------------------------------------
#############################################################################################################################

library(readr)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(ggpubr)

# path to data files
path_data <- paste0(getwd(), "/Data/")
dat_cases <- read_csv(paste0(path_data, "Delta-Omicron_frequencies.csv")) # saved original file as CSV to prevent date corruption
dat_freq_delta <- read_csv(paste0(path_data, "Growth-rates_first-detection_delta.csv"))
dat_freq_omicron <- read_csv(paste0(path_data, "Growth-rates_first-detection_omicron.csv"))

# path to save figures to
path_figures <- paste0(getwd(), "/Figures/")

# select colors for plotting
all_dark2_colors <- brewer.pal(8, "Dark2") 
customPalette <- c(all_dark2_colors[1], all_dark2_colors[3])

#############################################################################################################################
# Figure 1A - Variant Cases over Time  ----------------------------------------------
#############################################################################################################################

# Clarify column names and select only needed columns
dat_cases <- dat_cases[2:nrow(dat_cases), 1:9] # also drop 1st row - additional headers
colnames(dat_cases) <- c("Date", 
                         "num_delta_sgtf", "num_omicron_sgtf", 
                         "num_delta_ngs", "num_omicron_ngs",
                         "num_delta_total", "num_omicron_total",
                         "prop_delta_total", "prop_omicron_total")

dat_cases <- dat_cases %>% dplyr::select("Date", "num_delta_total", "num_omicron_total", "prop_omicron_total")

# Fix formats
dat_cases$Date <- as.Date(dat_cases$Date, format = "%m/%d/%y")
dat_cases$num_delta_total <- as.numeric(dat_cases$num_delta_total)
dat_cases$num_omicron_total <- as.numeric(dat_cases$num_omicron_total)
dat_cases$prop_omicron_total <- as.numeric(dat_cases$prop_omicron_total)

# Change format from wide to long to ready for formatting
dat_cases_long <- reshape2::melt(dat_cases, id.vars = "Date", measure.vars = c("num_delta_total", "num_omicron_total"))
colnames(dat_cases_long) <- c("Date", "Variant", "Num_Cases")
dat_cases_long$Variant <- as.character(dat_cases_long$Variant)
dat_cases_long[which(dat_cases_long$Variant == "num_delta_total"), "Variant"] <- "Delta"
dat_cases_long[which(dat_cases_long$Variant == "num_omicron_total"), "Variant"] <- "Omicron"
dat_cases_sub <- dat_cases[, c("Date", "prop_omicron_total")]
dat_cases_long <- dat_cases_long %>% dplyr::left_join(dat_cases_sub, by = "Date")

# Fig 1A - plot Delta vs Omicron cases (total) and percent Omicron
sec_axis_scale <- 1000 / 1

fig_1A <- ggplot(dat_cases_long, aes(x = Date, y = Num_Cases, fill = Variant)) +
  geom_bar(stat = "identity") +
  geom_line(aes(x = Date, y = prop_omicron_total*sec_axis_scale)) +
  scale_y_continuous(name = "Daily Cases (SGTF and NGS)", 
                     sec.axis = sec_axis(~./sec_axis_scale, name = "Percent Omicron", 
                                         labels = function(b) { paste0(round(b * 100, 0), "%")})) + 
  scale_fill_manual(values = customPalette) +
  theme_bw() +
  theme(legend.position = "bottom", legend.box = "horizontal")

#############################################################################################################################
# Figure 1B - Doubling Rates  ----------------------------------------------
#############################################################################################################################

# Clarify column names
dat_freq_delta <- dat_freq_delta %>% 
  dplyr::rename(Counter = `days since detection`,
                Value = Delta)

dat_freq_omicron <- dat_freq_omicron %>% 
  dplyr::rename(Counter = `days since detection`,
                Value = Omicron)

# Restrict data to same emergence period length
max_emergence_day <- min(max(dat_freq_delta$Counter), max(dat_freq_omicron$Counter), na.rm = TRUE)
dat_freq_delta_emerge <- dat_freq_delta %>% dplyr::filter(Counter <= max_emergence_day)
dat_freq_omicron_emerge <- dat_freq_omicron %>% dplyr::filter(Counter <= max_emergence_day)

# Calculate doubling time (days) - Delta
dat_freq_delta_emerge$Cumsum <- cumsum(dat_freq_delta_emerge$Value)
y_log <- log(dat_freq_delta_emerge$Cumsum)
x_days <- dat_freq_delta_emerge$Counter 
model <- lm(y_log ~ x_days)

## estimated doubling time
double_t_delta <- log(2) / model$coefficients[2]  
double_t_delta <- round(double_t_delta, digits = 2)

## estimated 95% confidence intervals for doubling time
model_CI <- confint(model)
double_t_delta_lo <- log(2) / model_CI[2, 1]
double_t_delta_lo <- round(double_t_delta_lo, digits = 2)
double_t_delta_hi <- log(2) / model_CI[2, 2]
double_t_delta_hi <- round(double_t_delta_hi, digits = 2)

# Calculate doubling time (days) - Omicron
dat_freq_omicron_emerge$Cumsum <- cumsum(dat_freq_omicron_emerge$Value)
y_log <- log(dat_freq_omicron_emerge$Cumsum)
x_days <- dat_freq_omicron_emerge$Counter 
model <- lm(y_log ~ x_days)

## estimated doubling time
double_t_omicron <- log(2) / model$coefficients[2]  
double_t_omicron <- round(double_t_omicron, digits = 2)

## estimated 95% confidence intervals for doubling time
model_CI <- confint(model)
double_t_omicron_lo <- log(2) / model_CI[2, 1]
double_t_omicron_lo <- round(double_t_omicron_lo, digits = 2)
double_t_omicron_hi <- log(2) / model_CI[2, 2]
double_t_omicron_hi <- round(double_t_omicron_hi, digits = 2)


#############################################################################################################################
# Figure 1B - Emergence Logistic Growth Curves ----------------------------------------------
#############################################################################################################################

# Change format so only 1 sample per line - either 0 or 1
expand_freq_data <- function(dat) {
  
  dat_expand <- NULL
  
  for (each_row in 1:nrow(dat)) {
    
    dat_row <- dat[each_row, ]
    
    # create line for each sample belonging to this variant
    dat_date <- rep(dat_row$Date, dat_row$Value)
    dat_counter <- rep(dat_row$Counter, dat_row$Value)
    dat_binary <- rep(1, dat_row$Value) # 1 = belongs
    dat_row_expand <- cbind.data.frame(dat_date, dat_counter, dat_binary)
    dat_expand <- rbind.data.frame(dat_expand, dat_row_expand)
    
    # create line for each sample NOT belonging to this variant
    dat_date <- rep(dat_row$Date, dat_row$Other)
    dat_counter <- rep(dat_row$Counter, dat_row$Other)
    dat_binary <- rep(0, dat_row$Other) # 0 = belongs to Other variant
    dat_row_expand <- cbind.data.frame(dat_date, dat_counter, dat_binary)
    dat_expand <- rbind.data.frame(dat_expand, dat_row_expand)
    
  }
  colnames(dat_expand) <- c("Date", "Counter", "Value")
  dat_expand
}

dat_freq_delta_emerge_exp <- expand_freq_data(dat_freq_delta_emerge)
dat_freq_omicron_emerge_exp <- expand_freq_data(dat_freq_omicron_emerge)

# Add variant category 
dat_freq_delta_emerge_exp$Variant <- rep("Delta", nrow(dat_freq_delta_emerge_exp))
dat_freq_omicron_emerge_exp$Variant <- rep("Omicron", nrow(dat_freq_omicron_emerge_exp))

# Order by counter
dat_freq_delta_emerge_exp <- dat_freq_delta_emerge_exp[order(dat_freq_delta_emerge_exp$Counter, decreasing = FALSE), ]
dat_freq_omicron_emerge_exp <- dat_freq_omicron_emerge_exp[order(dat_freq_omicron_emerge_exp$Counter, decreasing = FALSE), ]

# Combine Delta and Omicron together
dat_freq_both_emerge <- rbind.data.frame(dat_freq_delta_emerge_exp, dat_freq_omicron_emerge_exp)

# Plot logistic growth curves
plotLineage_cat_Xdays <- function(dat, num_days){
  alphapoint = 0.5
  alphaline = 0.8
  alphashade = 0.2
  ggplot(dat, aes(x = Counter, y = Value, color = Variant, fill = Variant)) + 
    theme_bw()+ 
    scale_fill_manual(values = customPalette) +
    scale_color_manual(values = customPalette) +
    geom_line(method = "glm", method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") +
    stat_smooth(method = "glm", method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0) +
    labs(color = "Variant", fill = "Variant") + 
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(0, num_days, by = 5)) +
    xlab("Days Since Initial Detection") +
    ylab("Probability of Given Sample Belonging to Variant") +
    theme(legend.position = "bottom", legend.box = "horizontal")
}

fig_1B <- plotLineage_cat_Xdays(dat_freq_both_emerge, num_days = max_emergence_day) +
  annotate("text", x = 4.9, y = 0.97, label = "Doubling time", fontface = 2) +
  annotate("text", x = 5.7, y = 0.92, label = paste0("Delta: ", double_t_delta, " days"), color = customPalette[1]) +
  annotate("text", x = 6.2, y = 0.87, label = paste0("Omicron: ", double_t_omicron, " days"), color = customPalette[2])

#############################################################################################################################
# Figure 1A/B - Combine Plots ----------------------------------------------
#############################################################################################################################

get_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- get_legend(fig_1B  + theme(legend.position='bottom'))

p <- grid.arrange(ggarrange(fig_1A + theme(legend.position="none"),
                            fig_1B + theme(legend.position="none"),
                            ncol = 2, labels = c("A", "B")),
                  mylegend,
                  nrow = 2,
                  heights = c(9.5, 0.5))

ggsave(paste(path_figures, "Fig_1.pdf", sep=""), p, height = 5, width =10)

#############################################################################################################################
# Figure 1A/B - Summary Stats ----------------------------------------------
#############################################################################################################################

# run binomial logistic regression separately so can pull out coefficients
dat <- dat_freq_both_emerge
model_list <- list()
CI_list <- list()
predict_list <- list()

for(i in 1:length(unique(dat$Variant))){ # for each variant category
  dat_var <- dat[which(dat$Variant== unique(dat$Variant)[[i]]), ]
  model_list[[i]] <- glm(Value ~ Counter, data=dat_var, family = "binomial") # for each variant, calculate binomial logistic regression
  CI_list[[i]] <- confint(model_list[[i]])[2, 1:2] # calculate confidence intervals around coefficients
  predict_logodds <- predict(model_list[[i]], newdata = data.frame(Counter = c(seq(0, max_emergence_day, by = 1)), type = "response")) # give fitted values for emergence period
  predict_odds <- exp(predict_logodds)
  predict_probs <- predict_odds / (1 + predict_odds) # transform log odds into probabilities
  predict_list[[i]] <- predict_probs
}

model_coefs <- sapply(model_list, function(x) coefficients(summary(x))[2,1]) # coefficient for X variable
coef_std_err <- sapply(model_list, function(x) coefficients(summary(x))[2,2]) # std error for X variable
model_p <- sapply(model_list, function(x) coefficients(summary(x))[2,4]) # p value for X variable
CI_vals <- sapply(CI_list, function(x) x[1:2]) # confidence intervals for X variable

model_df <- data.frame(lineage = unique(dat$Variant),
                       coefs = model_coefs,
                       lowci = CI_vals[1, ],
                       upci = CI_vals[2, ])

model_df$coefs <- round(model_df$coefs, digits = 2)
model_df$lowci <- round(model_df$lowci, digits = 2)
model_df$upci <- round(model_df$upci, digits = 2)


