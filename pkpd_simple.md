PKPD
================
Hans Brunner

## R Markdown

In this analysis, I explore pharmacokinetic (PK) and pharmacodynamic
(PD) data from a Multiple Ascending Dose (MAD) study using simulated
data from xgx’s PKPD Datasets. The focus is on understanding how drug
concentrations change over time (PK) and their effects on the body (PD),
without assessing safety aspects.

I first visualize the PK and PD data separately, followed by an analysis
of the dose-effect relationship to understand how different doses
influence the drug’s efficacy. This helps inform the optimal dosing
regimen while balancing effectiveness.

``` r
library(dplyr)
library(nlmixr2)
library(ggplot2)
library(patchwork)
library(knitr)
library(nlmixr2)
library(xgxr)

# Load the MAD dataset
pkpd_data <- read.csv("../Data/Multiple_Ascending_Dose_Dataset2.csv")


DOSE_CMT = 1
PK_CMT = 2
PD_CMT = 3
SS_PROFDAY = 6 # steady state prof day

# Prepare data
pkpd_data = pkpd_data %>%
  mutate(ID      = ID,     #ID   column
         TIME    = TIME,   #TIME column name 
         NOMTIME = NOMTIME,#NOMINAL TIME column name
         PROFDAY = 1 + floor(NOMTIME / 24), #PROFILE DAY day associated with profile, e.g. day of dose administration
         LIDV    = LIDV,   #DEPENDENT VARIABLE column name
         CENS    = CENS,   #CENSORING column name
         CMT     = CMT,    #COMPARTMENT column
         DOSE    = DOSE,   #DOSE column here (numeric value)
         TRTACT  = TRTACT, #DOSE REGIMEN column here (character, with units),
         LIDV_NORM = LIDV/DOSE,
         LIDV_UNIT    = EVENTU,
         DAY_label = ifelse(PROFDAY > 0, paste("Day", PROFDAY), "Baseline")
  )

#create a factor for the treatment variable for plotting
pkpd_data = pkpd_data %>%
  arrange(DOSE) %>%
  mutate(TRTACT_low2high = factor(TRTACT, levels = unique(TRTACT)),
         TRTACT_high2low = factor(TRTACT, levels = rev(unique(TRTACT))))

#create pk and pd datasets
pk_data <- pkpd_data %>%
  filter(CMT==PK_CMT)

pd_data <- pkpd_data %>%
  filter(CMT==PD_CMT)

pkpd_data_wide <- pd_data %>%
  select(ID, NOMTIME, PD = LIDV) %>%
  right_join(pk_data) %>%
  rename(CONC = LIDV)%>%
  filter(!is.na(PD))%>%
  filter(!is.na(CONC))

# labels and units
time_units_dataset = "hours"
time_units_plot    = "days"
conc_units         = unique(pk_data$LIDV_UNIT) %>% as.character()
conc_label         = paste0("Concentration (", conc_units, ")")
pd_units           = unique(pd_data$LIDV_UNIT) %>% as.character()
pd_label           = paste0("Continuous PD Marker (", pd_units, ")") 
```

## Including Plots

``` r
gg <- ggplot(data = pk_data, aes(x = NOMTIME, y = LIDV, color = TRTACT_high2low, fill = TRTACT_high2low)) +
  xgx_stat_ci(conf_level = 0.95) +
  xgx_scale_x_time_units(units_dataset = time_units_dataset, units_plot = time_units_plot) +
  xgx_scale_y_log10() +
  guides(color = guide_legend(""), fill = guide_legend("")) +
  labs(y = conc_label)
print(gg)
```

![](pkpd_simple_files/figure-gfm/plot%20concentration%20time%20curve-1.png)<!-- -->

``` r
gg %+% (data = pd_data) +
  scale_y_continuous() +
  labs(y = pd_label)
```

![](pkpd_simple_files/figure-gfm/plot%20PD%20marker-1.png)<!-- -->

``` r
max_pd <- pkpd_data %>%
  filter(CMT == 3) %>%
  group_by(ID, DOSE) %>%
  summarize(Max_Effect = max(LIDV, na.rm = TRUE), .groups = 'drop')
# Emax model to find minimum effective dose
emax_model <- nls(
  Max_Effect ~ (Emax * DOSE) / (EC50 + DOSE),
  data = max_pd,
  start = list(Emax = 4, EC50 = 50)
)

emax_params <- coef(emax_model)
Emax <- emax_params["Emax"]
EC50 <- emax_params["EC50"]

target_effect <- 0.8 * Emax

# Dose needed to achieve the target effect
med_dose <- (EC50 * target_effect) / (Emax - target_effect)
cat("Minimum Effective Dose (MED) for 80% of Emax:", round(med_dose, 2), "mg\n")
```

    ## Minimum Effective Dose (MED) for 80% of Emax: 213.05 mg

``` r
dose_seq <- seq(0, max(max_pd$DOSE), length.out = 100)
predicted_effect <- (Emax * dose_seq) / (EC50 + dose_seq)

plot_data <- data.frame(DOSE = dose_seq, Effect = predicted_effect)

# Plot the dose-response curve and MED
ggplot(plot_data, aes(x = DOSE, y = Effect)) +
  geom_line(color = "blue") +
  geom_point(data = max_pd, aes(x = DOSE, y = Max_Effect), color = "red") +
  geom_vline(xintercept = med_dose, linetype = "dashed", color = "green") +
  annotate("text", x = med_dose, y = 0.8 * Emax, 
           label = paste("MED =", round(med_dose, 2), "mg"), 
           vjust = -1, hjust = 1, color = "green") +
  labs(
    title = "Dose-Response Curve with Minimum Effective Dose (MED)",
    x = "Dose (mg)",
    y = "Effect"
  ) +
  theme_minimal()
```

![](pkpd_simple_files/figure-gfm/plot%20dose%20effect%20relationship-1.png)<!-- -->
