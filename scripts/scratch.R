df <- bioassay_data_toxprint %>% filter(PROTOCOL_NAME=="tox21-hdac-p1")



toxprint_only <- df |> select(-SMILES,-FractionCSP3,-ASSAY_OUTCOME,-drug,-InChiKey,-DTXSID,-PREFERRED_NAME,-PROTOCOL_NAME) 

toxprint_only <- toxprint_only[,colSums(toxprint_only != 0)!=0]

# Use numeric fields for the PCA
pca <- prcomp(toxprint_only[,unlist(lapply(toxprint_only, FUN=class))=="numeric"])

var_explained <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
pca_data <- as.data.frame(pca$x[,var_explained<0.98])

pcr_model_formula <- paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + ",paste0(colnames(pca_data),collapse = " + "))

pcr_mod <- glm(pcr_model_formula,
               data=bind_cols(df[,c("ASSAY_OUTCOME","drug","FractionCSP3")],pca_data),
               family=binomial(link='logit'))

#185 perfect seperation
tidy(glm(paste0("ASSAY_OUTCOME ~ 1 + drug + FractionCSP3 + ",paste0(colnames(pca_data)[1:163],collapse = " + ")),
    data=bind_cols(df[,c("ASSAY_OUTCOME","drug","FractionCSP3")],pca_data),
    family=binomial(link='logit')))



pcr_mod_summary <- tidy(pcr_mod) %>% filter(term=='FractionCSP3')
