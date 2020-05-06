#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 12:10:05 2020

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
        corr_samples1 = np.corrcoef(samples1[:,0],samples1[:,1])[0,1]
        corr_samples2 = np.corrcoef(samples2[:,0],samples2[:,1])[0,1]
        discrepancy[2] = abs(samples1_var[0,0] - samples2_var[0,0])
        discrepancy[3] = abs(samples1_var[1,1] - samples2_var[1,1])
        discrepancy[4] = abs(corr_samples1  - corr_samples2) 
        
    elif option == 'nonparametric':
        nParams = samples1.shape[1]  
        discrepancy = np.zeros((nParams,1))
        
        for i in range(nParams):
            discrepancy[i],_ = stats.ks_2samp(samples1[:,i],samples2[:,i])
    return discrepancy.squeeze()

#%%
N_trial_list = [20,50]
samples_list = [2000,3000,5000,10000]


    
# nonpara_discrepancy_apt = np.zeros((len(samples_list),2,9))
# nonpara_discrepancy_mh = np.zeros((len(samples_list),2,9))

para_discrepancy_apt = np.zeros((2,len(samples_list),9,5))
para_discrepancy_mh = np.zeros((2,len(samples_list),9,5))

for i,N_trial in enumerate(N_trial_list):
    resultsFolder = '../dscData/Psycho Results/CP Statistics/selfConsistencyOverParamDraws%d' % N_trial
    # load analytical results
    name = os.path.join(resultsFolder,'analytical.pk')
    with open(name, 'rb') as handle:
        file  = pk.load(handle)
    analytical_posterior = file['posterior_1']
    alpha_space = file['alpha']
    beta_space = file['beta']
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
    for j,N in enumerate(samples_list):
        

        
        
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
                    
                    #nonpara_discrepancy_apt[i,j,:,n1] = calculate_discrepancy(apt_samples,analytical_samples,option = 'nonparametric')
                    para_discrepancy_apt[i,j,n1,:] = calculate_discrepancy(apt_samples,analytical_samples,option = 'parametric')
                    n1 += 1                     
            elif file.endswith('.pk') & file.startswith('mh_samples'):
                with open(os.path.join(subfolder,file),'rb') as handle:
                    res = pk.load(handle)
                    mh_samples = res
                    
                    #nonpara_discrepancy_mh[i,j,:,n2] = calculate_discrepancy(mh_samples,analytical_samples,option = 'nonparametric')
                    para_discrepancy_mh[i,j,n2,:] = calculate_discrepancy(mh_samples,analytical_samples,option = 'parametric')
                    n2 += 1
#%%  
  
logpara_discrepancy_apt = np.log(para_discrepancy_apt)
logpara_discrepancy_mh = np.log(para_discrepancy_mh)
# from scipy.io import savemat
# savemat('discrepancy',{})         
# plt.figure()
# feature_labels = ['alpha','beta']
# for i in range(2):
#     plt.subplot(2,1,i+1)
#     plt.errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(nonpara_discrepancy_apt[:,i,:],axis = 1),\
#                  yerr = np.std(nonpara_discrepancy_apt[:,i,:],axis = 1))
#     plt.errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(nonpara_discrepancy_mh[:,i,:],axis = 1),\
#                  yerr = np.std(nonpara_discrepancy_mh[:,i,:],axis = 1))
#     plt.title(feature_labels[i])
#     plt.xticks(np.log(np.asarray(samples_list)*2),['4000','6000','10000','20000'])
#     if i == 0:
#        plt.legend(['APT','mh_sampling'])
fig,axs = plt.subplots(5,sharex = True,sharey = False)
fig.set_size_inches(6, 10)
fig.subplots_adjust(hspace=0.5, wspace=0.4)
feature_labels = ['sensitivity_mean','bias_mean','sensitivity_variance','bias_variance','covariance']
for i in range(5):

    axs[i].errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(logpara_discrepancy_apt[0,:,:,i],axis = 1),\
                  yerr = np.std(logpara_discrepancy_apt[0,:,:,i],axis = 1),label ='apt_20',color = 'b',ls= '--')
        
    axs[i].errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(logpara_discrepancy_mh[0,:,:,i],axis = 1),\
                  yerr = np.std(logpara_discrepancy_mh[0,:,:,i],axis = 1),label = 'mh_sampling_20',color = 'g',ls = '--')
    
    axs[i].errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(logpara_discrepancy_apt[1,:,:,i],axis = 1),\
                  yerr = np.std(logpara_discrepancy_apt[1,:,:,i],axis = 1),label ='apt_50',color = 'b')
        
    axs[i].errorbar(np.log(np.asarray(samples_list)*2),y = np.mean(logpara_discrepancy_mh[1,:,:,i],axis = 1),\
                  yerr = np.std(logpara_discrepancy_mh[1,:,:,i],axis = 1),label = 'mh_sampling_50',color='g')    
    axs[i].set_title(feature_labels[i],size = 12)
    axs[i].set_ylabel('Discrepancy')
    plt.xticks(np.log(np.asarray(samples_list)*2),['4000','6000','10000','20000'])
    if i == 4:
        plt.legend(loc='l', bbox_to_anchor=(1,1))

    plt.xlabel('Number of simulations',size = 16)