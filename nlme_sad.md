NLME SAD
================
Hans Brunner

## R Markdown

This report presents a pharmacokinetic (PK) analysis of single ascending
dose (SAD) data using a one-compartment model with first-order
absorption and elimination. The goal of this analysis is to estimate
clearance (CL), volume of distribution (Vd), and absorption rate
constant (Ka), using a nonlinear mixed-effects (NLME) approach.

Using the nlmixr2 package in R, I fit a one-compartment model to the
observed concentration-time data from the simulated SAD study. The
fitted model accounts for both fixed effects (population averages) and
random effects (inter-individual variability).

To assess the adequacy of the model, I generate various goodness-of-fit
plots, including observed versus predicted concentrations and residuals
analysis. These visual diagnostics help evaluate the model’s predictive
performance and its ability to accurately describe the PK data.

``` r
library(ggplot2)
library(dplyr)
library(tidyr)
library(xgxr)
library(readr)
library(caTools)
library(patchwork)
library(nlmixr2)
library(gt)

sad_data <- read_csv('../data/Single_Ascending_Dose_Dataset2.csv')

# Prepare data for one compartmental model fitting using nonlinear mixed effects model (nlmixr2)

one_comp_data <- sad_data %>%
  mutate(
    DV = ifelse(EVID == 1, 0, LIDV),  # Set DV to 0 when EVID is 1, otherwise use LIDV
    ID = as.integer(as.character(ID)) 
  ) %>%
  filter(!is.na(DV) & !is.nan(DV)) %>% # Discard rows where DV is NA or NaN
  select(TIME, ID, DV, AMT, EVID)     
```

``` r
one_comp_model<- function() {
  ini({
    tka <- log(0.3)  # Typical value of absorption rate constant (Ka)
    tcl <- log(0.02) # Typical value of clearance (CL)
    tv <- log(10)    # Typical value of volume of distribution (V)
    
    # Inter-individual variability
    eta.ka + eta.cl + eta.v ~ c(1, 
                                0.01, 1, 
                                0.01, 0.01, 1)
    add.err <- 0.1    # Residual variability
  })
  model({
    ka <- exp(tka + eta.ka)  # Individual value of absorption rate constant
    cl <- exp(tcl + eta.cl)  # Individual value of clearance
    v <- exp(tv + eta.v)     # Individual value of volume of distribution
    
    ke <- cl / v             # Elimination rate constant
    
    # Differential equations
    d/dt(A1) = -ka * A1                # Drug absorption from dosing compartment
    d/dt(A2) = ka * A1 - ke * A2       # Drug in central compartment (elimination)
    
    cp = A2 / v                        # Plasma concentration in central compartment
    cp ~ add(add.err)                  # Define error model
  })
}


fit <- nlmixr(one_comp_model, one_comp_data, "saem",
              control=list(print=0), 
              table=list(cwres=TRUE, npde=TRUE))
```

    ## [====|====|====|====|====|====|====|====|====|====] 0:00:00

    ## [====|====|====|====|====|====|====|====|====|====] 0:00:00

## Plots

Predicted vs Observed

``` r
ggplot(fit, aes(x = IPRED, y = DV)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    title = "One-Compartment Model: Predicted vs Observed",
    x = "Predicted Concentration",
    y = "Observed Concentration"
  )
```

![](nlme_sad_files/figure-gfm/plot%20predicted%20vs%20observed-1.png)<!-- -->

``` r
ggplot(fit, aes(x = TIME, y = RES)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  labs(
    title = "Residuals vs Time",
    x = "Time (hours)",
    y = "Residuals"
  ) +
  theme_minimal()
```

![](nlme_sad_files/figure-gfm/residuals%20vs%20time-1.png)<!-- -->

``` r
ggplot(fit, aes(x = IPRED, y = RES)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed") +
  labs(
    title = "Residuals vs Predicted Concentrations",
    x = "Predicted Concentration",
    y = "Residuals"
  ) +
  theme_minimal()
```

![](nlme_sad_files/figure-gfm/plot%20residuals%20vs%20predicted%20concentraion-1.png)<!-- -->

``` r
ggplot(fit, aes(x = TIME, y = DV, color = factor(ID))) +
  geom_point() +
  geom_line(aes(y = IPRED), linetype = "solid") +
  labs(
    title = "Concentration-Time Profiles (Observed vs Predicted)",
    x = "Time (hours)",
    y = "Concentration"
  ) +
  theme_minimal() +
  theme(legend.position = "none")
```

![](nlme_sad_files/figure-gfm/plot%20predicted%20concentration%20time%20curve-1.png)<!-- -->

``` r
set.seed(123) 
selected_ids <- sample(unique(fit$ID), 4)
subset_data <- fit %>% filter(ID %in% selected_ids)
ggplot(subset_data, aes(x = TIME, y = DV, color = factor(ID))) +
  geom_point() +
  geom_line(aes(y = IPRED), linetype = "solid") +
  labs(
    title = "Concentration-Time Profiles (Observed vs Predicted)",
    x = "Time (hours)",
    y = "Concentration"
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  facet_wrap(~ ID, ncol = 2)
```

![](nlme_sad_files/figure-gfm/plot%204%20psuedo-random%20examples-1.png)<!-- -->
