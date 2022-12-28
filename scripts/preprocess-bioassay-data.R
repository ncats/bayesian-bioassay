library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(purrr)
library(stringr)

#' Read all sheets from an excel file and return a dataframe
#'
#' @param xl_path The full path of the excel file
xlsheet_df <- function(xl_path){
  
  xldf <- xl_path |> 
    excel_sheets() |> 
    set_names() |> 
    map_dfr(read_excel,path=xl_path)
  
  return(xldf)
}

#' Convert bioassay data from wide format, with a column for each bioassay result
#' to a long format with a single column for with bioassay protocol names
#'
#' @param assay_df The dataframe containing wide format bioassay data
pivot_protocols <- function(assay_df){
  
  protocol_names <- colnames(assay_df)[str_detect(colnames(assay_df),'tox21')]
  
  assay_df_long <- assay_df |> 
    select(-`...1`) |> 
    select(SMILES,FractionCSP3,all_of(protocol_names)) |> 
    pivot_longer(cols=all_of(protocol_names),
                 names_to='PROTOCOL_NAME',
                 values_to='ASSAY_OUTCOME',
                 values_drop_na = TRUE)
  
  return(assay_df_long)
  
}

###Read in cell viability counter-screen data and annotations
# Annotation data from supplemental information of https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7887805/

xl_paths <-  list.files("data/50_assay/",full.names = TRUE)
tox21_compound_smiles <- read_csv("data/outcome_matrix_updated.csv") |> select(SAMPLE_ID,InChiKey,SMILES)
tox21_drugs <- read_xlsx("data/tx0c00264_si_003.xlsx",sheet="S3.TOX21SL list overlaps") |> select(DTXSID,DRUGS)
tox21_id_maps <- read_xlsx("data/tx0c00264_si_003.xlsx",sheet="S1.Tox21 ID Map")|> select(TOX21_ID,SAMPLE_ID,DTXSID)
tox21_toxprints <-  read_xlsx("data/tx0c00264_si_003.xlsx",sheet="S5.TOX21SL ToxPrints") 

### Reformat and join
#TODO:
# - wrap drug & toxprint to function for clarity

bioassay_data <- xl_paths |> 
  map(xlsheet_df) |> 
  map_dfr(pivot_protocols)

tox21_drug_status <- left_join(
  merge(tox21_id_maps,tox21_drugs),
  tox21_compound_smiles,on="SAMPLE_ID") |> 
  group_by(SMILES,InChiKey,DTXSID) |>
  summarize(drug = max(DRUGS))

bioassay_data <- left_join(
  bioassay_data,tox21_drug_status[,c("drug","SMILES","InChiKey")],
  on="SMILES")

tox21_dtx_smiles <- merge(tox21_id_maps,tox21_compound_smiles,on="SAMPLE_ID") |> 
  select(DTXSID,SMILES) |> unique()

tox21_toxprints <- left_join(tox21_toxprints,tox21_dtx_smiles,by="DTXSID") |> unique()

tox21_toxprints <- tox21_toxprints[,colSums(tox21_toxprints != 0)!=0]

colnames(tox21_toxprints) <- stringr::str_replace_all(colnames(tox21_toxprints),":","")

toxprint_names <- setdiff(colnames(tox21_toxprints),c("PREFERRED_NAME","SMILES","DTXSID"))

#Not every compound has a toxprint, so this dataset is slightly reduced
bioassay_data_toxprint <- right_join(bioassay_data,tox21_toxprints,on="SMILES")

write_csv(bioassay_data,"data/all_assays_with_drugstatus.csv")
write_csv(bioassay_data_toxprint,"data/all_assays_with_toxprints.csv")


