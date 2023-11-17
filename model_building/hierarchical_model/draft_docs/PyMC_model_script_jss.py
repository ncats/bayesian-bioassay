#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import multiprocessing as mp

import pymc as pm
import pandas as pd
import numpy as np
import arviz as az
#import aesara as ae
#import aesara.tensor as T
import pytensor.tensor as pt
from sklearn.model_selection import train_test_split
from sklearn.compose import ColumnTransformer
from sklearn import metrics
from sklearn.preprocessing import StandardScaler,LabelBinarizer
import pickle

#TODO:
# convert the sections below into readable functions
# increase the tune and draw values in the model
# add the prediction components back in
def main() -> int:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    n_features=int(sys.argv[1])
    logging.info(f"starting main with {n_features} features")
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

    factorized_protocols = pd.factorize(chem_des_scale.ProtocolName)

    #X = chem_des_scale.iloc[:, :204]
    X = chem_des_scale.iloc[:,:n_features]
    y = chem_des_scale['Outcome']
    assay_info_try = chem_des_scale.iloc[:,206:]

    X_train, X_test, y_train, y_test, assay_info_train, assay_info_test = train_test_split(X, y, assay_info_try, test_size=0.2, random_state=42)
    X4_bayes = X_train

    X4_bayes.insert(0,'Intercept',1)
    X4_bayes = np.asarray(X4_bayes, dtype="float64")

    factorized_protocols = pd.factorize(assay_info_train.ProtocolName)

    Y4_bayes = pd.factorize(y_train)

    coords_simulated = {
        'obs_id': np.arange(X4_bayes.shape[0]),
        'chem_descrip': np.arange(X4_bayes.shape[1]),
        'protocol':list(factorized_protocols[1]),
        'params':['beta_{0}'.format(i) for i in range(X4_bayes.shape[1])]
    }

    organisms = pd.factorize(assay_info_train.drop_duplicates().Organism)
    tissues_4 = pd.factorize(assay_info_train.drop_duplicates().Tissue_Type4)
    tissues_2 = pd.factorize(assay_info_train.drop_duplicates().Tissue_Type2)
    gender = pd.factorize(assay_info_train.drop_duplicates().Gender)
    cell_type = pd.factorize(assay_info_train.drop_duplicates().Cell_Type)

    #lambdas = pd.read_csv('Tox21/log_reg_w_lambda_sig.csv')['0'].tolist()

    logging.info('start model setup')
    with pm.Model(coords=coords_simulated) as assay_level_model:
        x = pm.Data('x', X4_bayes, mutable = True)
        protocol_idx = pm.Data("protocol_idx", list(factorized_protocols[0]), mutable=True)
        organism_idx = pm.Data("organism_idx", organisms[0], mutable=True)
        #tissue_2_idx = pm.Data("tissue_2_idx", tissues_2[0], mutable=True)
        #tissue_4_idx = pm.Data("tissue_4_idx", tissues_4[0], mutable=True)
        gender_idx = pm.Data("gender_idx", gender[0], mutable=True)
        #cell_type_idx = pm.Data("cell_type_idx", cell_type[0], mutable=True)

        # prior stddev in intercepts & slopes (variation across protocol):
        sd_dist = pm.Exponential.dist(1.0)
        
        # get back standard deviations and rho:
        chol, corr, stds = pm.LKJCholeskyCov("chol", n=X4_bayes.shape[1], eta=2.0, sd_dist=sd_dist, compute_corr=True)
        
        #hyperpriors and priors for average betas:
        beta_list = []
        for i in range(X4_bayes.shape[1]):
            gbeta = pm.Normal("g_beta_{0}".format(i), mu=0.0, sigma=10.0, shape=3)

            mu_gbeta = gbeta[0] + gbeta[1] * organism_idx + gbeta[2] * gender_idx
            sigma_beta = pm.Exponential('sigma_beta_{0}'.format(i),10) ## pass the saved lambda from the logistic regression
            #sigma_beta = np.sqrt(lambdas)
            betas = pm.Normal('beta_{0}'.format(i), mu=mu_gbeta,sigma=sigma_beta,dims="protocol")
            beta_list.append(betas)

        #population of varying protocol effects:
        beta_protocol = pm.MvNormal("beta_protocol", mu=pt.stack(beta_list, axis=1), chol=chol,dims=('protocol', 'params'))
        
        #Expected value per protocol:
        theta = beta_protocol[protocol_idx,0]* x[:,0]

        for i in range(1,X4_bayes.shape[1]):
            theta += beta_protocol[protocol_idx,i] * x[:,i]

        p = 1.0 / (1.0 + pt.exp(-theta))
        likelihood = pm.Bernoulli('likelihood', p, observed=Y4_bayes[0], shape = p.shape)
        

    logging.info('Done model setup')
    logging.info('Running model')

    # Run the model
    with assay_level_model:
        tr_assay = pm.sample(draws=1000, tune=1000,cores=8,chains=8,init="adapt_diag")


    filename = f"results/tr_8_assay_model_noextra.pkl"
    pickle.dump(tr_assay, open(filename, 'wb'))
    logging.info('model saved')
    return 0


if __name__=="__main__":
    sys.exit(main())
