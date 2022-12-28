library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)

balanced_model_coef <- read_csv("results/Balanced weight feature coefficients.csv") |> rename(feature_name = `...1` )
unbalanced_model_coef <- read_csv("results/Unbalanced weight feature coefficients.csv")|> rename(feature_name = `...1` )

single_assay_base_mod <- read_csv("results/unadjusted_fsp3_estimates.csv") |> mutate(model="Base model")
single_assay_drug_mod <- read_csv("results/drug_adjusted_fsp3_estimates.csv") |> mutate(model="Drug status adjusted")
single_assay_toxprint_pcr_mod_df <- read_csv("results/toxprint_adjusted_fsp3_estimates.csv") |> mutate(model="Drug status w/Toxprint PCR")
single_assay_top_01_tox_mod <- read_csv("results/toxprint_adjusted_fsp3_estimates_top01.csv") |> mutate(model="Drug status w/Top 82 Toxprints")
single_assay_top_05_tox_mod <- read_csv("results/toxprint_adjusted_fsp3_estimates_top05.csv") |> mutate(model="Drug status w/Top 200 Toxprints")

single_assay_results <- bind_rows(list(single_assay_base_mod,single_assay_drug_mod,single_assay_toxprint_pcr_mod_df,single_assay_top_01_tox_mod,single_assay_top_05_tox_mod))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(single_assay_results,aes(estimate,model,fill=model)) + 
  geom_boxplot() +
  scale_fill_manual(values=cbbPalette) +
  theme_classic() +guides(fill="none") + 
  xlab("Fsp3 regression coefficient estimate") +
  ylab("Model")

p
ggsave(plot = p,"fig/fsp3-single-assay-boxplot-summary.png",width = 6,height=4)


'tox21-ar-mda-kb2-luc-agonist-p3'
fsp3_coef_ranks_unbalanced <- unbalanced_model_coef  |> 
  mutate(across(-feature_name, ~rank(-abs(.x)))) |> 
  filter(feature_name == 'FractionCSP3') |>
  pivot_longer(-feature_name)

fsp3_coef_ranks_balanced <- balanced_model_coef |> 
  mutate(across(-feature_name, ~rank(-abs(.x)))) |> 
  filter(feature_name == 'FractionCSP3') |>
  pivot_longer(-feature_name)
