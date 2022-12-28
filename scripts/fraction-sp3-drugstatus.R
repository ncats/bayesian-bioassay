library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(lme4)

###Dependence of cell-tox & fsp3 adjusted for drug status

###TODO
# - only use drugs
# - explain the difference in results across assays, aka additional fixed effects 

#' Perform logistic regression for each assay
#' 
#' @param model_formula Specifies the model formula to use
single_assay_reg <- function(data,model_formula){
  
  model_results <- data |> 
    nest(data = -PROTOCOL_NAME) |>
    mutate(
      fit = map(data, ~ glm(model_formula,data=.x,family=binomial(link='logit'))),
      tidied = map(fit,tidy)
    ) |> 
    unnest(tidied) |> filter(term=='FractionCSP3') |> select(-fit,-data)
  
  return(model_results)  
}

bioassay_data <- read_csv("data/all_assays_with_drugstatus.csv")

## Model formulas, put here for ease of comparison
single_assay_base_form <- "ASSAY_OUTCOME ~ 1 + FractionCSP3"
single_assay_drug_form <- "ASSAY_OUTCOME ~ 1 + drug + FractionCSP3"
mixed_model_base_form <- "ASSAY_OUTCOME ~ 1 + FractionCSP3 + (1 + FractionCSP3 | PROTOCOL_NAME)"
mixed_model_drug_form <- "ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + (1 + FractionCSP3 | PROTOCOL_NAME)"

#Single assay models
single_assay_base_mod <- single_assay_reg(bioassay_data,single_assay_base_form)
single_assay_drug_mod <- single_assay_reg(bioassay_data,single_assay_drug_form)

#mixed effect models
mixed_effect_base_mod <- glmer(mixed_model_base_form,data=bioassay_data,family=binomial(link='logit'))
mixed_effect_drug_mod <- glmer(mixed_model_drug_form,data=bioassay_data,family=binomial(link='logit'))

#Plotting
p <- ggplot(single_assay_drug_mod,aes(reorder(PROTOCOL_NAME, estimate),estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymax=estimate+1.96*std.error,ymin=estimate-1.96*std.error,width=0.5)) +
  geom_hline(yintercept=fixef(mixed_effect_drug_mod)[3],color='red',linetype='dashed') +
  geom_hline(yintercept=0,color='black') +
  xlab("Protocol Name") +ylab("Effect of Fraction C-SP3") +
  coord_flip() + theme_classic() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
  )

ggsave(plot = p,"fig/fsp3-logistic-regression-with-drugstatus.png",width = 10,height=8)

### Mixed effect model tables
sjPlot::tab_model(mixed_effect_drug_mod,
                  pred.labels = c("Intercept","Drug status","Fraction C-SP3"),
                  title="Mixed effect logistic regression with drug status adjustment",
                  dv.labels="Cell viability counterscreen outcome",
                  file="fig/mixed_effect_drug_mod_tab.html")


write_csv(single_assay_base_mod,"results/unadjusted_fsp3_estimates.csv")
write_csv(single_assay_drug_mod,"results/drug_adjusted_fsp3_estimates.csv")



