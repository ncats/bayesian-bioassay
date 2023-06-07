#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import multiprocessing as mp

import pymc as pm
import pymc.sampling_jax as pmjax
import pandas as pd
import numpy as np
import arviz as az
import pytensor.tensor as pt
from sklearn.model_selection import train_test_split
from sklearn.compose import ColumnTransformer
from sklearn import metrics
from sklearn.preprocessing import StandardScaler,LabelBinarizer
import pickle

# ### Stack 50 assays together along with assay meta info
def main() -> int:
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
#    n_features=int(sys.argv[1])
    logging.info(f"starting main with 204 features")
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

#    assay_idx = [6,7,8,9,28,29,30,31]
#    meta_idx = [3,4,6]
#    assay_info = assay_info.reset_index(drop=True).iloc[assay_idx,meta_idx]  
   
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

#    X = chem_des_scale.iloc[:, :204]
#    X = chem_des_scale.iloc[:, :20]
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



    # In[83]:


    X_train, X_test, y_train, y_test, assay_info_train, assay_info_test = train_test_split(X, y, assay_info_try, test_size=0.2, random_state=42)
    X4_bayes = X_train


    # In[84]:


    X4_bayes.insert(0,'Intercept',1)
    X4_bayes = np.asarray(X4_bayes, dtype="float64")
    print(X4_bayes.shape)

    # In[85]:


    factorized_protocols = pd.factorize(assay_info_train.ProtocolName)

    Y4_bayes = pd.factorize(y_train)
    print('Y4_bayes length', len(Y4_bayes[0]))
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



    # lambdas = pd.read_csv('Tox21/log_reg_w_lambda_sig.csv')['0'].tolist()


    print('start model setup')
    with pm.Model(coords=coords_simulated) as assay_level_model:
        x = pm.Data('x', X4_bayes, mutable = True)
        protocol_idx = pm.Data("protocol_idx", list(factorized_protocols[0]), mutable=True)
        organism_idx = pm.Data("organism_idx", organisms[0], mutable=True)
        tissue_2_idx = pm.Data("tissue_2_idx", tissues_2[0], mutable=True)
        tissue_4_idx = pm.Data("tissue_4_idx", tissues_4[0], mutable=True)
        gender_idx = pm.Data("gender_idx", gender[0], mutable=True)
        cell_type_idx = pm.Data("cell_type_idx", cell_type[0], mutable=True)

        # prior stddev in intercepts & slopes (variation across protocol):
        sd_dist = pm.Exponential.dist(1.0)
        
        # get back standard deviations and rho:
        chol, corr, stds = pm.LKJCholeskyCov("chol", n=X4_bayes.shape[1], eta=2.0, sd_dist=sd_dist, compute_corr=True)
        
        #hyperpriors and priors for average betas:
        beta_list = []
        for i in range(X4_bayes.shape[1]):
            gbeta = pm.Normal("g_beta_{0}".format(i), mu=0.0, sigma=10.0, shape=6)

            mu_gbeta = gbeta[0] + gbeta[1] * organism_idx + gbeta[2] * tissue_2_idx + gbeta[3] * tissue_4_idx + gbeta[4] * gender_idx + gbeta[5] * cell_type_idx
            sigma_beta = pm.Exponential('sigma_beta_{0}'.format(i),1.0) ## pass the saved lambda from the logistic regression
    #        sigma_beta = np.sqrt(lambdas)
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
        

    print('Done model setup')
    logging.info('Running model')
    # Run the model
    with assay_level_model:
        tr_assay = pmjax.sample_numpyro_nuts(draws=500, tune=500,chains=4)
    print('Finish training model')

#    RANDOM_SEED = 10
#    with assay_level_model:
#        tr_assay.extend(pm.sample_posterior_predictive(tr_assay))


    filename = 'tr_assay_full_model_40features.pkl'
    pickle.dump(tr_assay, open(filename, 'wb'))
    print('model saved')



    print('generate testing data')
    predictors_out_of_sample = X_test.copy()
    predictors_out_of_sample.insert(0,'Intercept',1)
    predictors_out_of_sample = np.asarray(predictors_out_of_sample, dtype='float64')
    print(predictors_out_of_sample.shape)
    outcomes_out_of_sample = y_test
    print(outcomes_out_of_sample.shape)

    filename_pred_sample = 'predictors_out_of_sample.pkl'
    filename_out_sample = 'outcome_out_of_sample.pkl'
    pickle.dump(predictors_out_of_sample, open(filename_pred_sample, 'wb'))
    pickle.dump(outcomes_out_of_sample, open(filename_out_sample, 'wb'))
    print('X_test, y_test saved')

    # In[7]:


    RANDOM_SEED = 10
    ### generate sample for prediction
    pd.DataFrame(assay_info_test).to_csv('assay_info_test.csv', index=False)
    print('assay_info_test saved')
    protocol_pred = pd.factorize(assay_info_test.ProtocolName)
    
    prediction_coords = {
        'obs_id': np.arange(predictors_out_of_sample.shape[0]),
        'chem_descrip': np.arange(predictors_out_of_sample.shape[1]),
        'protocol':list(protocol_pred[1]),
        'params':['beta_{0}'.format(i) for i in range(predictors_out_of_sample.shape[1])]
    }

    print('run prediction on current model')
    with assay_level_model:
        pm.set_data(
            new_data = {"x": predictors_out_of_sample,
                        "protocol_idx": list(protocol_pred[0]),
                        "organism_idx": organisms,
                        "tissue_2_idx": tissues_2,
                        "tissue_4_idx": tissues_4,
                        "gender_idx": gender,
                        "cell_type_idx": cell_type
                        }
        )

        posterior_predictive = pm.sample_posterior_predictive(
            tr_assay, var_names=["likelihood"], random_seed=RANDOM_SEED
        )


    # In[8]:


    # posterior_predictive.posterior_predictive['likelihood'].mean(('chain', 'draw'))
    test_pred = pd.DataFrame(posterior_predictive.posterior_predictive['likelihood'].mean(('chain', 'draw')))
    pd.DataFrame(test_pred).to_csv('test_pred_full_model.csv', index=False)
    pd.DataFrame([RANDOM_SEED]).to_csv('random_seed_full_model.csv', index=False)



    # In[9]:


    # y_test_r = y_test.replace('inactive', 0.0)
    # y_test_r = y_test_r.replace('active antagonist', 1.0)
    # y_test_r
    print('visualize y_test')
    print(y_test)

    # In[10]:

    print('classification report for testing set')
    print(metrics.classification_report(list(y_test), pd.DataFrame(test_pred).round()))
    print('Balanced Accuracy Test')
    print(metrics.balanced_accuracy_score(list(y_test), pd.DataFrame(test_pred).round()))


    # In[11]:


    print('AUCROC Test')
    print(metrics.roc_auc_score(list(y_test), test_pred))
    return 0

if __name__=="__main__":
    sys.exit(main())




