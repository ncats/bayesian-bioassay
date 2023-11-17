#!/usr/bin/env python
# coding: utf-8

import sys
import logging

import pandas as pd
import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,LabelBinarizer
import pickle

def main() -> int:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)

    xls = pd.ExcelFile('Tox21/assay_list.xls')

    # to read just one sheet to dataframe:
    df = []
    for i in range(len(xls.sheet_names)):
        df1 = xls.parse(i)
        df.append(df1)

    xls2 = pd.ExcelFile('Tox21/assay_list2.xls')

    for i in range(len(xls2.sheet_names)):
        df1 = xls2.parse(i)
        df.append(df1)

    logging.info("read in main data")

    t = pd.DataFrame({'Cell_Line': ['DT40','DT40', 'DT40'], 'ProtocolName':['tox21-dt40-p1_100', 'tox21-dt40-p1_653', 'tox21-dt40-p1_657']})

    assay_info_tox21 = pd.read_csv('Tox21/assay_meta_data.csv')
    assay_info_tox21 = assay_info_tox21.drop(assay_info_tox21.columns[0], axis = 1)
    assay_info_tox21['Cell_Line'] = assay_info_tox21['Cell_Line'].str.replace(r'*', '', regex=True)
    assay_info_tox21 = assay_info_tox21[['Cell_Line', 'ProtocolName']]
    assay_info_tox21 = assay_info_tox21.drop([49])
    assay_info_tox21 = pd.concat([assay_info_tox21, t]).reset_index(drop=True)

    assay_info = pd.read_excel('Tox21/SampleMeta_Data_update.xlsx')
    assay_info = assay_info.drop(assay_info.columns[6:], axis=1).iloc[:13,:].rename(columns = {assay_info.columns[5]: 'Tissue_Type2'})
    assay_info = assay_info.rename(columns = {assay_info.columns[0]: 'Cell_Line'})
    assay_info['Cell_Line'] = assay_info['Cell_Line'].str.replace(r'*', '', regex=True).str.replace(' ', '')
    assay_info['Cell_Type'] = assay_info['Cell_Type'].str.replace(r'*', '', regex=True)
    assay_info = assay_info.join(assay_info_tox21.set_index('Cell_Line'), on='Cell_Line')
   
    select_assays = ['tox21-shh-3t3-gli3-agonist-p1', 'tox21-shh-3t3-gli3-antagonist-p1', 'tox21-rar-antagonist-p2', 'tox21-dt40-p1_657', 'tox21-vdr-bla-agonist-p1','tox21-vdr-bla-antagonist-p1','tox21-hdac-p1', 'tox21-p53-bla-p1']

    df_list = []
    for i in range(len(df)):
        one_assay = df[i]
        one_assay = one_assay.drop(one_assay.columns[0], axis=1).drop_duplicates()
        if one_assay.columns[0] in select_assays:
            one_assay['ProtocolName'] = one_assay.columns[0]
            print(one_assay.columns[0])
            one_assay = one_assay.rename(columns = {one_assay.columns[0]: 'Outcome'})
            df_list.append(one_assay)

    df_stack_50 = pd.concat(df_list)
    df_stack_50 = df_stack_50.join(assay_info.set_index('ProtocolName'), on='ProtocolName')
    
    logging.info(f"Data shape is {df_stack_50.shape}")

    logging.info("completed pre-processing")
    # ### Train/Test Split & Bayesian model coordination

    df_stack_50.Outcome = df_stack_50.Outcome.astype(object)

    # Scale numeric chemical descriptors
    numeric_columns=list(df_stack_50.select_dtypes(['float64', 'int64']).columns)
    categorical_columns=list(df_stack_50.select_dtypes('object').columns)

    pipeline=ColumnTransformer([
        ('num',StandardScaler(),numeric_columns),
        ('cat', 'passthrough', categorical_columns)])

    chem_des_scale = pipeline.fit_transform(df_stack_50)
    chem_des_scale = pd.DataFrame(chem_des_scale)

    col_name = list(df_stack_50.columns)
    col_name.pop(0)
    col_name.pop(0)
    col_name.insert(204, 'Outcome')
    col_name.insert(205,'SMILES')
    chem_des_scale.columns = col_name

    chem_des_scale.to_pickle("data/processed_data.pkl")
    logging.info('data saved')
    return 0


if __name__=="__main__":
    sys.exit(main())
