###############################################################################################
### Replication data for: Attitudes toward Internal Migrants and Support for Redistribution ###
### Authors: Hsu Yumin Wang and Eddy S. F. Yeung                                            ###
### Date: November 8, 2022                                                                  ###
###############################################################################################

### Set-up ----
## Clean the R environment and set the working directory
rm(list = ls())
setwd("~/Desktop/mig-and-redist/analysis")

## Load the required packages
library(haven)
library(tidyverse)
library(estimatr)
library(texreg)

## Import the dataset
df <- read_dta("main.dta")

## Remove duplicated respondents
df <- df %>% filter(duplicate_location == 0 & duplicate_ip == 0)

## Remove respondents without Shanghai hukou
df$hukou <- ifelse(df$hukou == 1, 0, 1)
df <- df %>% filter(hukou == 1)

### Recode variables ----
## Female (= 1)
df$female <- df$gender - 1
table(df$female)

## Han (= 1)
df$han <- ifelse(df$han == 1, 1, 0)
table(df$han)

## Pretreatment attitude toward internal migration (-3 = decrease greatly; 3 = increase greatly)
df <- df %>% 
  mutate(
    mig_view = case_when(
      mig_view_dec == 1 ~ -3,
      mig_view_dec == 2 ~ -2,
      mig_view_neu == 2 ~ -1,
      mig_view_neu == 3 ~ 0,
      mig_view_neu == 1 ~ 1,
      mig_view_inc == 2 ~ 2,
      mig_view_inc == 1 ~ 3
    )
  )
table(df$mig_view)

## Treatment group (1 = control; 2 = fiscal burden; 3 = cultural difference;
## 4 = labor competition; 5 = less exclusionary)
df <- df %>% 
  mutate(
    treatment = case_when(
      group == "Very Exclusionary - No Prime"             ~ 1,
      group == "Very Exclusionary - Fiscal Burden"        ~ 2,
      group == "Very Exclusionary - Cultural Differences" ~ 3,
      group == "Very Exclusionary - Labor Market Threat"  ~ 4,
      group == "Less Exclusionary - No Prime"             ~ 5
    )
  )
df$treatment <- as.factor(df$treatment)
table(df$treatment)
df$any_treatment <- ifelse(df$treatment == 1, 0, 1)
df$group2 <- ifelse(df$treatment == 2, 1, ifelse(df$treatment == 1, 0, NA))
df$group3 <- ifelse(df$treatment == 3, 1, ifelse(df$treatment == 1, 0, NA))
df$group4 <- ifelse(df$treatment == 4, 1, ifelse(df$treatment == 1, 0, NA))
df$group5 <- ifelse(df$treatment == 5, 1, ifelse(df$treatment == 1, 0, NA))

## Redistribution measure 1 (-3 = strongly oppose; 3 = strongly support)
df$redist1 <- df$redist1 - 4
table(df$redist1)

## Redistribution measure 2 (-3 = strongly disagree; 3 = strongly agree)
df$redist2 <- df$redist2 - 4
table(df$redist2)

## Redistribution measure 3 (-3 = strongly disagree; 3 = strongly agree)
df$redist3 <- df$redist3 - 4
table(df$redist3)

## Attention level (0 = lowest; 2 = highest)
df$attention2 <- ifelse(df$attention2 == 5, 1, 0)
df$attention <- df$attention1 + df$attention2
table(df$attention)

## Support for redistribution (mean effects index with normalization)
redist_matrix <- cbind(df$redist1, df$redist2, df$redist3)
cal_MEI <- 
  function(Z, outcome_mat) {
    c_mean <- apply(X = outcome_mat[Z == 0, ], MARGIN = 2, FUN = mean, na.rm = T)
    c_sd <- apply(X = outcome_mat[Z == 0, ], MARGIN = 2, FUN = sd, na.rm = T)
    z_score <- t(t(sweep(outcome_mat, 2, c_mean)) / c_sd)
    index_numerator <- rowSums(z_score)
    n_outcomes <- ncol(outcome_mat)
    index <- index_numerator / n_outcomes
    index <- (index - mean(index[Z == 0], na.rm = T)) / sd(index[Z == 0], na.rm = T)
    return(index)
  }
df <- df %>% 
  mutate(redist = cal_MEI(Z = df$any_treatment, outcome_mat = redist_matrix))
pdf("DV_distribution.pdf", width = 6, height = 4)
hist(df$redist, main = "",
     xlab = "Support for Redistribution", ylab = "Number of Respondents",
     family = "Times")
dev.off()

### Analysis ----
## Mean support for redistribution by experimental group
summary <- df %>% 
  group_by(treatment) %>% 
  do(tidy(lm_robust(redist ~ 1, data = .))) %>% 
  mutate(support = estimate)
p <- ggplot(summary, aes(x = treatment, y = support)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), size = .9, width = 0) +
  theme_bw() +
  scale_x_discrete(labels = c("1" = "Group 1\n(n = 609)",
                              "2" = "Group 2\n(n = 584)",
                              "3" = "Group 3\n(n = 602)",
                              "4" = "Group 4\n(n = 587)",
                              "5" = "Group 5\n(n = 633)")) +
  xlab("") +
  ylab("Average Support for Redistribution") +
  coord_cartesian(ylim = c(-.5, .5)) +
  theme(text = element_text(family = "Times", size = 14),
        axis.text = element_text(family = "Times", size = 13))
ggsave(file = "fig_avg_support.pdf", p, width = 8, height = 5)

## Difference-in-means between treatment and control groups
# Subset the dataset by experimental group
df_group1 <- subset(df, treatment == 1)
df_group2 <- subset(df, treatment == 2)
df_group3 <- subset(df, treatment == 3)
df_group4 <- subset(df, treatment == 4)
df_group5 <- subset(df, treatment == 5)

# Group 1 vs. Group 2
t.test(df_group1$redist, df_group2$redist)

# Group 1 vs. Group 3
t.test(df_group1$redist, df_group3$redist)

# Group 1 vs. Group 4
t.test(df_group1$redist, df_group4$redist)

# Group 1 vs. Group 5
t.test(df_group1$redist, df_group5$redist)

### Why null? Testing the theoretical prior ----
## Are anti-migrant attitudes and redistribution support positively correlated in Shanghai?
df$anti_mig <- -1 * df$mig_view
reg1 <- lm_robust(redist ~ anti_mig,
                  data = df)
reg2 <- lm_robust(redist ~ anti_mig + factor(age) + gender + han,
                  data = df)
reg3 <- lm_robust(redist ~ anti_mig + factor(age) + gender + han + 
                    educ + income + hu_fluency,
                  data = df)
texreg(list(reg1, reg2, reg3),
       stars = c(0.01, 0.05, 0.10),
       include.ci = F,
       custom.header = list("Dependent Variable: Support for Redistribution" = 1:3),
       custom.coef.names = 
         c("Constant", "Anti-Migrant Attitudes", "Age: 30 to 39", "Age: 40 to 49",
           "Age: 50 to 59", "Age: 60 or above", "Female", "Han", "Education", 
           "Household Income", "Hu Fluency"),
       custom.note = "Entries are OLS estimates with robust standard errors in parentheses.
       All significance tests are two-tailed with the following notations:
       $^{***}p<0.01$; $^{**}p<0.05$; $^{*}p<0.1$.",
       fontsize = "small")
