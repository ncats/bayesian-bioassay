library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(corrr)
library(ggplot2)
library(lme4)

###Dependence of cell-tox & fsp3 adjusted for know toxic features (ToxPrints)

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


bioassay_data_toxprint <- read_csv("data/all_assays_with_toxprints.csv") |> 
  na.omit()

single_assay_toxprint_pcr_mod <- by(
  bioassay_data_toxprint,
  INDICES=bioassay_data_toxprint$PROTOCOL_NAME, FUN=function(df){
    
    toxprint_only <- df |> select(-SMILES,-FractionCSP3,-ASSAY_OUTCOME,-drug,-InChiKey,-DTXSID,-PREFERRED_NAME,-PROTOCOL_NAME) 
    
    # Use numeric fields for the PCA
    pca <- prcomp(toxprint_only[,unlist(lapply(toxprint_only, FUN=class))=="numeric"])
    
    var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
    pca_data <- as.data.frame(pca$x[,var_explained<0.90])
    
    pcr_model_formula <- paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + ",paste0(colnames(pca_data),collapse = " + "))
    
    pcr_mod <- glm(pcr_model_formula,
                   data=bind_cols(df[,c("ASSAY_OUTCOME","drug","FractionCSP3")],pca_data),
                   family=binomial(link='logit'))
    
    pcr_mod_summary <- tidy(pcr_mod) %>% filter(term=='FractionCSP3')
    
    return(pcr_mod_summary)
  })

single_assay_toxprint_pcr_mod_df <-  do.call("rbind",single_assay_toxprint_pcr_mod)
single_assay_toxprint_pcr_mod_df$PROTOCOL_NAME <- names(single_assay_toxprint_pcr_mod)

p <- ggplot(single_assay_toxprint_pcr_mod_df,aes(reorder(PROTOCOL_NAME, estimate),estimate)) + 
  geom_point() + 
  geom_errorbar(aes(ymax=estimate+1.96*std.error,ymin=estimate-1.96*std.error,width=0.5)) +
  #geom_hline(yintercept=fixef(mixed_effect_drug_mod)[3],color='red',linetype='dashed') +
  geom_hline(yintercept=0,color='black') +
  xlab("Protocol Name") +ylab("Effect of Fraction C-SP3") +
  coord_flip() + theme_classic() + 
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=16),
  )

ggsave(plot = p,"fig/fsp3-logistic-regression-with-drugstatus-toxprint-pcr.png",width = 10,height=8)
write_csv(single_assay_toxprint_pcr_mod_df,"results/toxprint_adjusted_fsp3_estimates.csv")


#Most correlated toxprints
sp3_toxprint <- bioassay_data_toxprint |>
  select(-PROTOCOL_NAME,-ASSAY_OUTCOME,-drug,-InChiKey,-DTXSID,-PREFERRED_NAME) |>
  distinct()

sp3_toxprint_cor <- correlate(sp3_toxprint[,-1])

N = nrow(sp3_toxprint)
r = sp3_toxprint_cor[1,3:625]
t_stat = r*sqrt((N-2)/(1-r^2))
p_values <- pt(as.numeric(abs(t_stat)),df=N-2,lower.tail=FALSE)

top_01_toxprint <- names(sp3_toxprint)[which(abs(r)>0.1)]
top_05_toxprint <- names(sp3_toxprint)[which(abs(r)>0.05)]

top_01_tox_form <- paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + `",paste0(top_01_toxprint,collapse = "` + `"),"`")
top_05_tox_form <- paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + `",paste0(top_05_toxprint,collapse = "` + `"),"`")

#Single assay models
single_assay_top_01_tox_mod <- single_assay_reg(bioassay_data_toxprint,top_01_tox_form)
single_assay_top_05_tox_mod <- single_assay_reg(bioassay_data_toxprint,top_05_tox_form)

write_csv(single_assay_top_01_tox_mod,"results/toxprint_adjusted_fsp3_estimates_top01.csv")
write_csv(single_assay_top_05_tox_mod,"results/toxprint_adjusted_fsp3_estimates_top05.csv")

top_01_tox_form_mm <- paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + `",paste0(top_01_toxprint,collapse = "` + `"),"` + (1 + FractionCSP3 | PROTOCOL_NAME)")
mixed_effect_top_01_tox_mod <- glmer(top_01_tox_form_mm,data=bioassay_data_toxprint,family=binomial(link='logit'))
### Mixed effect model tables
sjPlot::tab_model(top_01_tox_form_mm,
                  pred.labels = c("Intercept","Drug status","Fraction C-SP3"),
                  title="Mixed effect logistic regression with drug status adjustment",
                  dv.labels="Cell viability counterscreen outcome",
                  file="fig/mixed_effect_toxprint_mod_tab.html")
