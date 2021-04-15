# /restricted/projectnb/fhs-methylation/Projects_H/DNAm_Mortality/DeepSurv/results_allcause/raw

import os
cwd = os.getcwd()

import sys
sys.path.append('/restricted/projectnb/fhs-methylation/Projects_H/DNAm_Mortality/DeepSurv/deepsurv')
import deep_surv

from deepsurv_logger import DeepSurvLogger, TensorboardLogger
import utils
import viz

import numpy as np
import pandas as pd

import lasagne
import matplotlib
import matplotlib.pyplot as plt

import numpy
import sklearn

# methy sig

def dataframe_to_deepsurv_ds(df, event_col = 'Event', time_col = 'Time', data_col='CpGs'):
    # Extract the event and time columns as numpy arrays
    e = df[event_col].values.astype(np.int32)
    t = df[time_col].values.astype(np.float32)
    # Extract the patient's covariates as a numpy array
    #x_df = df.drop([event_col, time_col], axis = 1)	
    x = df[data_col].values.astype(np.float32)
    # Return the deep surv dataframe
    return {
        'x' : x,
        'e' : e,
        't' : t
    }

meta_file="/restricted/projectnb/fhs-methylation/Projects_H/copy_fhsge/methy_huan/Mortality/Pred_2020_v1/meta_mortality_EA_results_Jan2020_repFHSARIC_disP05.csv"
meta= pd.read_csv(meta_file)

Pvals=[1e-7,1e-6,1e-5,1e-4,1e-3,0.05]
#Pvals=[1e-7,1e-6]
#p=Pvals[0]

i=10
pwd_file="/restricted/projectnb/fhs-methylation/Projects_H/copy_fhsge/methy_huan/Mortality/Pred_2020_v1/split_data_resid/"
Train_name=pwd_file + "Train_Pheno_methy_1304sample_invResid_p05_file" +str(i)+".csv"
Test_name=pwd_file+"Test_Pheno_methy_877sample_invResid_p05_file"+str(i)+".csv"

train_df= pd.read_csv(Train_name)
#train_df.head()
test_df= pd.read_csv(Test_name)
#test_df.head()
test_pheno=test_df.iloc[:,0:20]

## age, sex
#col_pheno=["age","sex","gramsperday","smoke","bmi","diabetes", "prevCHD","hrtnow","cancer","g2ex8phy", "EduYears","prevCHF","prevABI"]
col_pheno=["age","sex"]

train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'death', time_col= 'days_to_event_num',data_col=col_pheno)
test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'death', time_col= 'days_to_event_num',data_col=col_pheno)

hyperparams = {
        'L2_reg': 10.0,
        'batch_norm': True,
        'dropout': 0.4,
        'hidden_layers_sizes': [50, 50],
        'learning_rate': 0.02,
        'lr_decay': 0.001,
        'momentum': 0.9,
        'n_in': train_data['x'].shape[1],
        'standardize': False,
        'activation': 'selu'
    }

model = deep_surv.DeepSurv(**hyperparams)
logger = None
update_fn=lasagne.updates.nesterov_momentum
n_epochs = 500
metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
predRisk_test = model.predict_risk(test_data['x'])	
c_col='AgeSex_file' + str(i)
sLength = len(test_pheno['age'])
test_pheno[c_col]=predRisk_test	


## age, sex, clini
col_pheno=["age","sex","gramsperday","smoke","bmi","diabetes", "prevCHD","hrtnow","cancer","g2ex8phy", "EduYears","prevCHF","prevABI"]
#col_pheno=["age","sex"]

train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'death', time_col= 'days_to_event_num',data_col=col_pheno)
test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'death', time_col= 'days_to_event_num',data_col=col_pheno)

hyperparams = {
        'L2_reg': 10.0,
        'batch_norm': True,
        'dropout': 0.4,
        'hidden_layers_sizes': [50, 50],
        'learning_rate': 0.02,
        'lr_decay': 0.001,
        'momentum': 0.9,
        'n_in': train_data['x'].shape[1],
        'standardize': False,
        'activation': 'selu'
    }

model = deep_surv.DeepSurv(**hyperparams)
logger = None
update_fn=lasagne.updates.nesterov_momentum
n_epochs = 500
metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
predRisk_test = model.predict_risk(test_data['x'])	
c_col='AgeSexClin_file' + str(i)
sLength = len(test_pheno['age'])
test_pheno[c_col]=predRisk_test	


for p in Pvals:
    meta1=meta.loc[meta['dis.pval.random'] < p]
    
    CpGs=meta1["CpG"].values.astype('string')
    ##train_df[CpGs]
    # read data
    train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    #
    hyperparams = {
        'L2_reg': 10.0,
        'batch_norm': True,
        'dropout': 0.4,
        'hidden_layers_sizes': [200, 100, 50],
        'learning_rate': 0.02,
        'lr_decay': 0.001,
        'momentum': 0.9,
        'n_in': train_data['x'].shape[1],
        'standardize': False,
        'activation': 'selu'
    }
    # Create an instance of DeepSurv using the hyperparams defined above
    model = deep_surv.DeepSurv(**hyperparams)
    logger = None
    # Now we train the model
    update_fn=lasagne.updates.nesterov_momentum # The type of optimizer to use. \
                                                # Check out http://lasagne.readthedocs.io/en/latest/modules/updates.html \
                                                # for other optimizers to use
    n_epochs = 500
    # If you have validation data, you can add it as the second parameter to the function
    metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
    #print('Train C-Index:', metrics['c-index'][-1])
    predRisk_test = model.predict_risk(test_data['x'])
		
    c_col='CpGs_p'+ str(p)+'_file' + str(i)
    sLength = len(test_pheno['age'])
    test_pheno[c_col]=predRisk_test	
    ## age, sex, cpgs
    valuesInsert=["age","sex"]
    CpGs =np.append(CpGs,valuesInsert,axis=0)
    train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    hyperparams = {
        'L2_reg': 10.0,
        'batch_norm': True,
        'dropout': 0.4,
        'hidden_layers_sizes': [200, 100, 50],
        'learning_rate': 0.02,
        'lr_decay': 0.001,
        'momentum': 0.9,
        'n_in': train_data['x'].shape[1],
        'standardize': False,
        'activation': 'selu'
    }
	# Create an instance of DeepSurv using the hyperparams defined above
    model = deep_surv.DeepSurv(**hyperparams)
    logger = None
    # Now we train the model
    update_fn=lasagne.updates.nesterov_momentum # The type of optimizer to use. 
    n_epochs = 500
    # If you have validation data, you can add it as the second parameter to the function
    metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
	#print('Train C-Index:', metrics['c-index'][-1])
    predRisk_test = model.predict_risk(test_data['x'])
    #	
    c_col='AgeSex_CpGs_p'+ str(p)+'_file' + str(i)
    sLength = len(test_pheno['age'])
    test_pheno[c_col]=predRisk_test
	#
	## age, sex, clini, cpgs
    valuesInsert=["age","sex","gramsperday","smoke","bmi","diabetes", "prevCHD","hrtnow","cancer","g2ex8phy", "EduYears","prevCHF","prevABI"]
    CpGs =np.append(CpGs,valuesInsert,axis=0)
    train_data = dataframe_to_deepsurv_ds(train_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    test_data = dataframe_to_deepsurv_ds(test_df, event_col = 'death', time_col= 'days_to_event_num',data_col=CpGs)
    hyperparams = {
        'L2_reg': 10.0,
        'batch_norm': True,
        'dropout': 0.4,
        'hidden_layers_sizes': [200, 100, 50],
        'learning_rate': 0.02,
        'lr_decay': 0.001,
        'momentum': 0.9,
        'n_in': train_data['x'].shape[1],
        'standardize': False,
        'activation': 'selu'
    }
	# Create an instance of DeepSurv using the hyperparams defined above
    model = deep_surv.DeepSurv(**hyperparams)
    logger = None
	# Now we train the model
    update_fn=lasagne.updates.nesterov_momentum # The type of optimizer to use. \
    n_epochs = 500
    # If you have validation data, you can add it as the second parameter to the function
    metrics = model.train(train_data, n_epochs=n_epochs, logger=logger, update_fn=update_fn)
	#print('Train C-Index:', metrics['c-index'][-1])
    predRisk_test = model.predict_risk(test_data['x'])
		
    c_col='AgeSexClini_CpGs_p'+ str(p)+'_file' + str(i)
    sLength = len(test_pheno['age'])
    test_pheno[c_col]=predRisk_test

	

pwd_out='/restricted/projectnb/fhs-methylation/Projects_H/copy_fhsge/methy_huan/Mortality/Pred_2020_v1/dl/results_allsample/'	
out_name=pwd_out + 'testData_resid_predvalues_file' + str(i) +'.csv'
test_pheno.to_csv(out_name,index=False)
