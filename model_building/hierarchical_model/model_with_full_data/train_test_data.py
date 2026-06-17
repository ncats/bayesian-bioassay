#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import multiprocessing as mp

import pymc as pm
import pandas as pd
import numpy as np
import arviz as az
import pytensor.tensor as pt
from sklearn.model_selection import train_test_split
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder
from sklearn import metrics
from sklearn.preprocessing import StandardScaler,LabelBinarizer
import pickle

# ### Stack 50 assays together along with assay meta info
def main() -> int:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    logging.info(f"save training and testing for hierarchical model")
    # ### Stack 50 assays together along with assay meta info
    mp.set_start_method('forkserver', force=True)


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

    protocols = [df[i].columns[1] for i in range(len(df))]

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
    assay_info = assay_info[assay_info.ProtocolName.isin(protocols)].reset_index(drop=True)

   
    df_list = []
    for i in range(len(df)):
        one_assay = df[i]
        one_assay = one_assay.drop(one_assay.columns[0], axis=1).drop_duplicates(subset=['SMILES'])
        one_assay['ProtocolName'] = one_assay.columns[0]
        one_assay = one_assay.rename(columns = {one_assay.columns[0]: 'Outcome'})
        df_list.append(one_assay)
        

    df_stack_50 = pd.concat(df_list)
    df_stack_50 = df_stack_50.join(assay_info.set_index('ProtocolName'), on='ProtocolName')
    df_stack_50 = df_stack_50[df_stack_50.ProtocolName.isin(assay_info.ProtocolName)]

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


    factorized_protocols = pd.factorize(chem_des_scale.ProtocolName)

    feature_top_40 = ['Chi4v', 'SMR_VSA4', 'EState_VSA5', 'SlogP_VSA10', 'EState_VSA8',
       'Kappa3', 'fr_Al_COO', 'SlogP_VSA8', 'PEOE_VSA10', 'fr_amide',
       'PEOE_VSA14', 'VSA_EState10', 'FractionCSP3', 'PEOE_VSA4',
       'VSA_EState9', 'VSA_EState3', 'MolLogP', 'VSA_EState7', 'BalabanJ',
       'fr_ketone_Topliss', 'fr_methoxy', 'MaxAbsPartialCharge',
       'MinAbsPartialCharge', 'FpDensityMorgan3', 'BCUT2D_MWHI',
       'MinAbsEStateIndex', 'fr_ester', 'fr_bicyclic', 'SlogP_VSA3',
       'BCUT2D_LOGPLOW', 'BCUT2D_MRHI', 'fr_aniline', 'fr_para_hydroxylation',
       'EState_VSA2', 'VSA_EState4', 'EState_VSA9', 'PEOE_VSA9', 'SlogP_VSA1',
       'SlogP_VSA12', 'TPSA']


    X = chem_des_scale[chem_des_scale.columns.intersection(feature_top_40)]

    y = chem_des_scale['Outcome']
    assay_info_try = chem_des_scale.iloc[:,206:]

    ohe = OneHotEncoder(categories='auto')
    feature_arr = ohe.fit_transform(assay_info_try[['Tissue_Type4', 'Cell_Type', 'Gender', 'Organism', 'Tissue_Type2']]).toarray()
    feature_labels = ohe.get_feature_names_out(['Tissue_Type4', 'Cell_Type', 'Gender', 'Organism', 'Tissue_Type2'])

    meta_features = pd.DataFrame(feature_arr, columns=list(feature_labels))
    meta_features['ProtocolName'] = assay_info_try['ProtocolName']

    df_X = pd.concat([X,meta_features], axis=1)


    # In[83]:


#    X_train, X_test, y_train, y_test = train_test_split(df_X, y, test_size=0.2, random_state=42,
#                                                        stratify=df_X['ProtocolName'], shuffle=True)
    X_train, X_test, y_train, y_test, assay_info_train, assay_info_test = train_test_split(X, y, assay_info_try, test_size=0.2, random_state=42, stratify=assay_info_try['ProtocolName'], shuffle=True)

    print('X_train info')
    print(X_train.shape)
    print(X_train.columns)
    
    print(X_test)

    X_train.to_csv('X_train.csv', index=False)
    y_train.to_csv('y_train.csv', index=False)
    X_test.to_csv('X_test.csv', index=False)
    y_test.to_csv('y_test.csv', index=False)
    assay_info_train.to_csv('assay_info_train.csv', index=False)
    assay_info_test.to_csv('assay_info_test.csv', index=False)
    print('X_train, y_train, X_test, y_test saved')

    Y4_bayes = pd.factorize(y_train)
    print('Y4_bayes length', len(Y4_bayes[0]))

    return 0

if __name__=="__main__":
    sys.exit(main())





