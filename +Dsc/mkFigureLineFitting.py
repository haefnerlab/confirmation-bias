#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  4 09:41:32 2020

@author: liushizhao
"""

import numpy as np
import pickle as pk
import os
import matplotlib.pyplot as plt
from scipy import stats

def calculate_discrepancy(samples1,samples2,option = 'parametric'):
    if option == 'parametric':
        nParams = samples1.shape[1]
        nFeatures = 5
        
        discrepancy = np.zeros((nFeatures,1))
        discrepancy[0] = abs(np.mean(samples1[:,0]) - np.mean(samples2[:,0]))
        discrepancy[1] = abs(np.mean(samples1[:,1]) - np.mean(samples2[:,1]))
        samples1_var = np.cov(samples1.T)
        samples2_var = np.cov(samples2.T)
        discrepancy[2] = abs(samples1_var[0,0] - samples2_var[0,0])
        discrepancy[3] = abs(samples1_var[1,1] - samples2_var[1,1])
        discrepancy[4] = abs(samples1_var[0,1] - samples2_var[0,1]) 
        
    elif option == 'nonparametric':
        nParams = samples1.shape[1]  
        discrepancy = np.zeros((nParams,1))
        
        for i in range(nParams):
            discrepancy[i],_ = stats.ks_2samp(samples1[:,i],samples2[:,i])
    return discrepancy.squeeze()
# load analytical results
name = 'LineFitting/analytical.pk'
with open(name, 'rb') as handle:
    file  = pk.load(handle)
analytical_posterior = file[0]
alpha_space = file[1]
beta_space = file[2]
# normalize the posterior, from log space to probabilty space summing up to 1
analytical_posterior = np.exp(analytical_posterior)
analytical_posterior = analytical_posterior /np.sum(analytical_posterior[:])
# draw samples from analytical posterior
# create a population of alph&beta
alphav, betav = np.meshgrid(alpha_space , beta_space)

k = 10000
analytical_posterior_flat = analytical_posterior.flatten()
index = np.random.choice(np.arange(len(analytical_posterior_flat)),size = k,replace = True,p = analytical_posterior_flat)
analytical_samples = np.zeros((k,2))
analytical_samples[:,0] = alphav.flatten()[index]
analytical_samples[:,1] = betav.flatten()[index]
#%%
N_trial = 1000
samples_list = [2000,3000,5000,10000]

resultsFolder = '../dscData/Line Fitting Results/Results-a/selfConsistencyOverParamDraws'

nonpara_discrepancy_apt = np.zeros((len(samples_list),2,9))
nonpara_discrepancy_mh = np.zeros((len(samples_list),2,9))

para_discrepancy_apt = np.zeros((len(samples_list),5,9))
para_discrepancy_mh = np.zeros((len(samples_list),5,9))
for i,N in enumerate(samples_list):
    # load apt results
    subfolder = os.path.join(resultsFolder,'%d_%d_2/'%(N_trial,N))
    n1 = 0
    n2 = 0
    for file in os.listdir(subfolder):
    
        if file.endswith('.pk') & file.startswith('posterior_apt'):
            with open(os.path.join(subfolder,file),'rb') as handle:
                res = pk.load(handle)

    
                n_rounds = len(file)
                posterior = res[-1]
                apt_samples = posterior.gen(k)
                
                nonpara_discrepancy_apt[i,:,n1] = calculate_discrepancy(apt_samples,analytical_samples,option = 'nonparametric')
                para_discrepancy_apt[i,:,n1] = calculate_discrepancy(apt_samples,analytical_samples,option = 'parametric')
                n1 += 1                     
        elif file.endswith('.pk') & file.startswith('mh_samples'):
            with open(os.path.join(subfolder,file),'rb') as handle:
                res = pk.load(handle)
                mh_samples = res
                
                nonpara_discrepancy_mh[i,:,n2] = calculate_discrepancy(mh_samples,analytical_samples,option = 'nonparametric')
                para_discrepancy_mh[i,:,n2] = calculate_discrepancy(mh_samples,analytical_samples,option = 'parametric')
                n2 += 1
#%%  
plt.figure()
feature_labels = ['alpha','beta']
for i in range(2):
    plt.subplot(2,1,i+1)
    plt.errorbar(np.asarray(samples_list)*2,y = np.mean(nonpara_discrepancy_apt[:,i,:],axis = 1),\
                 yerr = np.std(nonpara_discrepancy_apt[:,i,:],axis = 1))
    plt.errorbar(np.asarray(samples_list)*2,y = np.mean(nonpara_discrepancy_mh[:,i,:],axis = 1),\
                 yerr = np.std(nonpara_discrepancy_mh[:,i,:],axis = 1))
    plt.title(feature_labels[i])
    plt.xticks(np.asarray(samples_list)*2,['4000','6000','10000','20000'])
    if i == 0:
       plt.legend(['APT','mh_sampling'])
plt.figure()
feature_labels = ['alpha_mean','beta_mean','alpha_variance','beta_variance','covariance']
for i in range(5):
    plt.subplot(5,1,i+1)
    plt.errorbar(np.asarray(samples_list)*2,y = np.mean(para_discrepancy_apt[:,i,:],axis = 1),\
                 yerr = np.std(para_discrepancy_apt[:,i,:],axis = 1))
    plt.errorbar(np.asarray(samples_list)*2,y = np.mean(para_discrepancy_mh[:,i,:],axis = 1),\
                 yerr = np.std(para_discrepancy_mh[:,i,:],axis = 1))
    plt.title(feature_labels[i])
    plt.xticks(np.asarray(samples_list)*2,['4000','6000','10000','20000'])
    if i == 0:
       plt.legend(['APT','mh_sampling'])