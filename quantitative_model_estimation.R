library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(gridExtra)


file <- file.path("Stan_models/chagas_efficacy_model_CT_values.stan")
mod <- cmdstan_model(file)
