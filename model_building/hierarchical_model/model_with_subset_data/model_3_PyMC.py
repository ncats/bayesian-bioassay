#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pymc as pm
import pandas as pd
import numpy as np
import arviz as az
import aesara as ae
import aesara.tensor as T
import pytensor.tensor as pt
import pickle
from sklearn import metrics
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler,LabelBinarizer
from sklearn.model_selection import train_test_split



# In[2]:


pymc_data = pd.read_csv('Tox21/pymc_data.csv')
pymc_data


# In[3]:


pymc_data.loc[pymc_data["outcome"] == 1, "outcome"] = 'active antagonist'
pymc_data.loc[pymc_data["outcome"] == 0, "outcome"] = 'inactive'


# In[5]:


assay_info = pymc_data.iloc[:, 53:]



# ### Setup for bayesian model

# In[7]:


assay_info['ProtocolName'] = pymc_data['ProtocolName']



# In[8]:


assay_info.describe()


# In[9]:


chem_des = pymc_data.iloc[:, 1:53]

# In[10]:


chem_des.outcome = chem_des.outcome.astype(object)


# In[12]:


numeric_columns=list(chem_des.select_dtypes(['float64', 'int64']).columns)
categorical_columns=list(chem_des.select_dtypes('object').columns)
print(len(numeric_columns))

pipeline=ColumnTransformer([
    ('num',StandardScaler(),numeric_columns),
    ('cat', 'passthrough', categorical_columns)])

chem_des_scale = pipeline.fit_transform(chem_des)
chem_des_scale = pd.DataFrame(chem_des_scale)


# In[13]:


col_name = list(chem_des.columns)
col_name.remove('ProtocolName')
col_name.remove('outcome')
col_name.append('ProtocolName')
col_name.append('outcome')
print(col_name)


# In[14]:


chem_des_scale.columns = col_name
chem_des_scale.head()


# In[15]:


factorized_protocols = pd.factorize(chem_des.ProtocolName)

X = chem_des_scale.iloc[:, :50]
y = chem_des_scale['outcome']


# In[16]:


X_train, X_test, y_train, y_test, assay_info_train, assay_info_test = train_test_split(X, y, assay_info, test_size=0.2, random_state=42)
# X4_bayes = X_train[['BCUT2D_CHGHI', 'BCUT2D_MRLOW', 'Chi0', 'Chi0n']] # descrip with higher importance rate
X4_bayes = X_train.iloc[:,:10]
# X4_bayes = X_train
print(X4_bayes.shape)
print(assay_info_train.shape)


# In[17]:


X4_bayes.insert(0,'Intercept',1)
X4_bayes = np.asarray(X4_bayes, dtype="float64")
print(X4_bayes.dtype)


# In[18]:


factorized_protocols = pd.factorize(assay_info_train.ProtocolName)
len(factorized_protocols[0])
# assay_info_


# In[19]:


Y4_bayes = pd.factorize(y_train)
# proto_name = 
#['tox21-ar-bla-antagonist-p1', 'tox21-gh3-tre-antagonist-p1', 'tox21-ahr-p1', 'tox21-erb-bla-p1']



# In[20]:


assay_info_train.iloc[0:200, :]


# In[21]:


print(X4_bayes.shape)


# In[22]:


organisms = [0,0,1,0]
tissues = [0,0,1,1]
gender = [0,0,1,0]
tissue_4 = [1,1,0,0]

# In[27]:


coords_simulated = {
    'obs_id': np.arange(X4_bayes.shape[0]),
    'chem_descrip': np.arange(X4_bayes.shape[1]),
    'protocol':list(factorized_protocols[1]),
    'params':['beta_{0}'.format(i) for i in range(X4_bayes.shape[1])]
}
coords_simulated


# In[158]:


#X4_bayes[:,0]
#X4_bayes.shape[0]


# In[188]:


lambdas = [1, 0.35938137,0.00599484,2.7825594,2.7825594]


# In[189]:

with pm.Model(coords=coords_simulated) as assay_level_model:
    x = pm.Data('x', X4_bayes, mutable = True)
    protocol_idx = pm.Data("protocol_idx", list(factorized_protocols[0]), mutable=True)
    organism_idx = pm.Data("organism_idx", organisms, mutable=True)
#    lambda_sig = pm.Data('lambda_sig', lambdas, mutable=True)
    tissue_idx = pm.Data('tissue_idx', tissues, mutable=True)
    gender_idx = pm.Data('gender_idx', gender, mutable=True)
    tissue_4_idx = pm.Data('tissue_4_idx', tissue_4, mutable=True)

    # prior stddev in intercepts & slopes (variation across protocol):
    sd_dist = pm.Exponential.dist(1.0)
    
    # get back standard deviations and rho:
    chol, corr, stds = pm.LKJCholeskyCov("chol", n=X4_bayes.shape[1], eta=2.0, sd_dist=sd_dist, compute_corr=True)
    
    #hyperpriors and priors for average betas:
    beta_list = []
    for i in range(X4_bayes.shape[1]):
        gbeta = pm.Normal("g_beta_{0}".format(i), mu=0.0, sigma=10.0, shape=5)

        mu_gbeta = gbeta[0] + gbeta[1] * organism_idx + gbeta[2] * tissue_idx + gbeta[3]*gender_idx+gbeta[4]*tissue_4_idx
        sigma_beta = pm.Exponential('sigma_beta_{0}'.format(i),1.0) ## pass the saved lambda from the logistic regression
#        sigma_beta = np.sqrt(lambda_sig[i])
        # print(sigma_beta)
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


# In[190]:


pm.model_to_graphviz(assay_level_model)


# In[191]:


with assay_level_model:
    tr_assay = pm.sample(1000, tune=2000, init="adapt_diag",chains=4,cores=8)


# In[192]:


with assay_level_model:
    az.plot_forest(
        tr_assay,
        combined=True,
        var_names=["beta_protocol"],
        figsize=(10, 10),
        textsize=14,
    )


# In[193]:


RANDOM_SEED = 10
with assay_level_model:
    tr_assay.extend(pm.sample_posterior_predictive(tr_assay))


# In[194]:

filename = 'tr_assay.pkl'
pickle.dump(tr_assay, open(filename, 'wb'))
#tr_assay


# In[195]:


pd.DataFrame(tr_assay.posterior_predictive['likelihood'].mean(('chain', 'draw'))).round().value_counts()


# In[134]:


# predictors_out_of_sample = X_test[['BCUT2D_CHGHI', 'BCUT2D_MRLOW', 'Chi0', 'Chi0n']]
predictors_out_of_sample = X_test.iloc[:,:20]
predictors_out_of_sample.insert(0,'Intercept',1)
predictors_out_of_sample = np.asarray(predictors_out_of_sample, dtype='float64')
print(predictors_out_of_sample.shape)
outcomes_out_of_sample = y_test
print(outcomes_out_of_sample.shape)


# In[197]:


RANDOM_SEED = 10
### generate sample for prediction
protocol_pred = pd.factorize(assay_info_test.ProtocolName)
    
with assay_level_model:
    pm.set_data(
        new_data = {"x": predictors_out_of_sample,
                    "protocol_idx": list(protocol_pred[0]),
                    "organism_idx": organisms,
	            "tissue_idx": tissues,
		    "gender_idx": gender,
		    "tissue_4_idx": tissue_4}
#                   "lambda_sig": lambdas}
    )

    posterior_predictive = pm.sample_posterior_predictive(
        tr_assay, var_names=["likelihood"], random_seed=RANDOM_SEED
    )


# In[198]:


posterior_predictive.posterior_predictive['likelihood'].mean(('chain', 'draw'))
test_pred = pd.DataFrame(posterior_predictive.posterior_predictive['likelihood'].mean(('chain', 'draw')))
pd.DataFrame(test_pred).to_csv('test_pred.csv', index=False)
pd.DataFrame([RANDOM_SEED]).to_csv('random_seed.csv', index=False)

# In[199]:

y_test_r = y_test.replace('inactive', 0.0)
y_test_r = y_test_r.replace('active antagonist', 1.0)
y_test_r


# In[200]:


print(metrics.classification_report(list(y_test_r), pd.DataFrame(test_pred).round()))
metrics.balanced_accuracy_score(list(y_test_r), pd.DataFrame(test_pred).round())


# In[201]:


print(metrics.roc_auc_score(list(y_test_r), test_pred))


train_pred = [0.165 , 0.1805, 0.37  , 0.389 , 0.678 , 0.6065, 0.442 , 0.1355,
              0.258 , 0.5585, 0.4955, 0.352 , 0.2505, 0.473 , 0.681 , 0.4945,
              0.8075, 0.5455, 0.409 , 0.514 , 0.3485, 0.289 , 0.4575, 0.6545,
       0.806 , 0.078 , 0.863 , 0.1995, 0.907 , 0.2135, 0.3875, 0.638 ,
       0.36  , 0.7415, 0.742 , 0.442 , 0.8345, 0.631 , 0.3655, 0.604 ,
       0.517 , 0.5395, 0.441 , 0.3765, 0.431 , 0.563 , 0.105 , 0.2885,
       0.742 , 0.522 , 0.4095, 0.256 , 0.367 , 0.6775, 0.337 , 0.699 ,
       0.318 , 0.4805, 0.349 , 0.193 , 0.2055, 0.367 , 0.2155, 0.433 ,
       0.647 , 0.746 , 0.262 , 0.5455, 0.5955, 0.5215, 0.28  , 0.2775,
       0.386 , 0.5325, 0.267 , 0.2865, 0.3295, 0.295 , 0.301 , 0.739 ,
       0.187 , 0.3945, 0.5745, 0.921 , 0.7165, 0.9885, 0.285 , 0.6825,
       0.188 , 0.773 , 0.418 , 0.6865, 0.265 , 0.743 , 0.331 , 0.109 ,
       0.618 , 0.4985, 0.645 , 0.677 , 0.236 , 0.582 , 0.1835, 0.671 ,
       0.7835, 0.646 , 0.2835, 0.3345, 0.3795, 0.5165, 0.806 , 0.4525,
       0.9195, 0.7785, 0.4085, 0.384 , 0.3605, 0.3435, 0.249 , 0.178 ,
       0.925 , 0.4705, 0.3135, 0.2695, 0.615 , 0.6945, 0.5535, 0.2835,
       0.199 , 0.236 , 0.8535, 0.6455, 0.5655, 0.661 , 0.553 , 0.609 ,
       0.646 , 0.529 , 0.3875, 0.481 , 0.4925, 0.906 , 0.5125, 0.163 ,
       0.642 , 0.7435, 0.3325, 0.999 , 0.438 , 0.405 , 0.165 , 0.3075,
       0.664 , 0.4265, 0.1825, 0.5935, 0.325 , 0.825 , 0.633 , 0.759 ,
       0.3415, 0.158 , 0.557 , 0.2465, 0.1775, 0.8805, 0.548 , 0.318 ,
       0.3285, 0.8395, 0.761 , 0.258 , 0.515 , 0.947 , 0.237 , 0.5455,
       0.446 , 0.8505, 0.423 , 0.4985, 0.558 , 0.642 , 0.429 , 0.229 ,
       0.415 , 0.263 , 0.4465, 0.5275, 0.3485, 0.846 , 0.4435, 0.445 ,
       0.627 , 0.1495, 0.712 , 0.284 , 0.1495, 0.2895, 0.9605, 0.1035,
       0.348 , 0.5885, 0.805 , 0.331 , 0.3205, 0.5365, 0.302 , 0.679 ,
       0.773 , 0.337 , 0.493 , 0.165 , 0.402 , 0.407 , 0.2545, 0.3745,
       0.5275, 0.1785, 0.3615, 0.461 , 0.383 , 0.6405, 0.493 , 0.387 ,
       0.6495, 0.4685, 0.4015, 0.4135, 0.7635, 0.4435, 0.942 , 0.637 ,
       0.4045, 0.8525, 0.3045, 0.3085, 0.304 , 0.374 , 0.397 , 0.5915,
       0.5065, 0.463 , 0.2135, 0.378 , 0.689 , 0.511 , 0.811 , 0.52  ,
       0.31  , 0.5815, 0.0855, 0.3225, 0.5005, 0.3235, 0.9555, 0.174 ,
       0.4935, 0.2015, 0.545 , 0.5945, 0.374 , 0.2285, 0.483 , 0.8605,
       0.7645, 0.6475, 0.1585, 0.6485, 0.6175, 0.294 , 0.4575, 0.983 ,
       0.478 , 0.597 , 0.381 , 0.4515, 0.767 , 0.151 , 0.6455, 0.4985,
       0.7905, 0.633 , 0.637 , 0.2465, 0.7745, 0.797 , 0.929 , 0.86  ,
       0.283 , 0.4   , 0.976 , 0.8885, 0.295 , 0.455 , 0.66  , 0.686 ,
       0.5425, 0.387 , 0.8405, 0.2625, 0.632 , 0.2655, 0.146 , 0.203 ,
       0.5445, 0.4035, 0.509 , 0.6905, 0.888 , 0.1045, 0.3045, 0.4525,
       0.998 , 0.743 , 0.352 , 0.3645, 0.4795, 0.9115, 0.1785, 0.1715]


# In[69]:


pd.DataFrame(train_pred).round().value_counts()


# In[78]:


y_train = y_train.replace('inactive', 0.0)
y_train = y_train.replace('active antagonist', 1.0)
y_train


# In[86]:



print(metrics.classification_report(list(y_train), pd.DataFrame(train_pred).round()))





