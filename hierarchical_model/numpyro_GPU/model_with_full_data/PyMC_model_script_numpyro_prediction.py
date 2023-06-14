#!/usr/bin/env python
# coding: utf-8

import sys
import logging
import multiprocessing as mp

import pymc as pm
#import pymc.sampling_jax as pmjax
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
    logging.info(f"try prediction with updated pymc package")
    # ### Stack 50 assays together along with assay meta info
    mp.set_start_method('forkserver', force=True)


    X_train = pd.read_csv('X_train.csv')
    X_test = pd.read_csv('X_test.csv')
    y_train = pd.read_csv('y_train.csv')
    y_test = pd.read_csv('y_test.csv')
    assay_info_train = pd.read_csv('assay_info_train.csv')
    assay_info_test = pd.read_csv('assay_info_test.csv')
    X4_bayes = X_train
    X_train_v, X_validation, y_train_v, y_validation, assay_info_train_v, assay_info_validation = train_test_split(X_train, y_train, assay_info_train, test_size = 0.2, random_state=42)


    # In[84]:
    X4_bayes.insert(0,'Intercept',1)
    X4_bayes = np.asarray(X4_bayes, dtype="float64")
    print(X4_bayes.shape)

    # In[85]:


    factorized_protocols = pd.factorize(assay_info_train.ProtocolName)

#    Y4_bayes = pd.factorize(y_train)
    print(y_train)
    Y4_bayes = list(y_train['Outcome'])
#    print(Y4_bayes[0].shape)
#    print('Y4_bayes length', len(Y4_bayes[0]))
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
            betas = pm.Normal('beta_{0}'.format(i), mu=mu_gbeta,sigma=sigma_beta,dims="protocol")
            beta_list.append(betas)

        #population of varying protocol effects:
        beta_protocol = pm.MvNormal("beta_protocol", mu=pt.stack(beta_list, axis=1), chol=chol,dims=('protocol', 'params'))
        
        #Expected value per protocol:
        theta = beta_protocol[protocol_idx,0]* x[:,0]

        for i in range(1,X4_bayes.shape[1]):
            theta += beta_protocol[protocol_idx,i] * x[:,i]

        p = 1.0 / (1.0 + pt.exp(-theta))
        likelihood = pm.Bernoulli('likelihood', p, observed=Y4_bayes, shape = p.shape)
        

    print('Done model setup')
    logging.info('Running model')
    # Run the model
    with assay_level_model:
#        tr_assay = pickle.load(open("tr_assay_full_model_40features.pkl", "rb"))
        tr_assay = pd.read_pickle("tr_assay_full_model_40features.pkl")  
    print('read saved trained model')

#    RANDOM_SEED = 10
#    with assay_level_model:
#        tr_assay.extend(pm.sample_posterior_predictive(tr_assay))



    print('generate testing data')
    predictors_out_of_sample = X_test.copy()
    print(predictors_out_of_sample)
    X_test.insert(0,'Intercept',1)
#    X_validation = np.asarray(X_validation, dtype='float')
    X_test = np.array(X_test)
    print(X_test.shape)
#    print(X_validation.dtype)
    print(y_test)
    outcomes_out_of_sample = list(y_test['Outcome'])
#    print(outcomes_out_of_sample)
#    print(outcomes_out_of_sample.shape)
    
    print(X4_bayes)
    print(X_test)

    # In[7]:


    RANDOM_SEED = 10
    ### generate sample for prediction
    protocol_pred = pd.factorize(assay_info_test.ProtocolName)
    
    prediction_coords = {
        'obs_id': np.arange(predictors_out_of_sample.shape[0]),
        'chem_descrip': np.arange(predictors_out_of_sample.shape[1]),
        'protocol':list(protocol_pred[1]),
        'params':['beta_{0}'.format(i) for i in range(predictors_out_of_sample.shape[1])]
    }


    organisms = pd.factorize(assay_info_test.drop_duplicates().Organism)
    tissues_4 = pd.factorize(assay_info_test.drop_duplicates().Tissue_Type4)
    tissues_2 = pd.factorize(assay_info_test.drop_duplicates().Tissue_Type2)
    gender = pd.factorize(assay_info_test.drop_duplicates().Gender)
    cell_type = pd.factorize(assay_info_test.drop_duplicates().Cell_Type)

    print('run prediction on current model')
    with assay_level_model:
        pm.set_data(
            new_data = {'x': X_test,
                        "protocol_idx": list(protocol_pred[0]),
                        "organism_idx": organisms[0],
                        "tissue_2_idx": tissues_2[0],
                        "tissue_4_idx": tissues_4[0],
                        "gender_idx": gender[0],
                        "cell_type_idx": cell_type[0]
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


    print('visualize y_test')
    print(y_test.Outcome.value_counts())
    print(test_pred)

    print('classification report for testing set')
    print(metrics.classification_report(list(y_test['Outcome']), pd.DataFrame(test_pred).round()))
    print('Balanced Accuracy Test')
    print(metrics.balanced_accuracy_score(list(y_test['Outcome']), pd.DataFrame(test_pred).round()))


    # In[11]:


    print('AUCROC Test')
    print(metrics.roc_auc_score(list(y_test['Outcome']), test_pred))
    return 0

if __name__=="__main__":
    sys.exit(main())





