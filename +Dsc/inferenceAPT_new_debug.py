#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  3 09:49:11 2020

@author: liushizhao
"""






#%% 
import pickle
import matlab.engine
from delfi.simulator.BaseSimulator import BaseSimulator
import os

import delfi.distribution as dd
import delfi.generator as dg
import numpy as np
import delfi.inference as infer
from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
import copy
from mat4py import loadmat
import logging
import datetime
from scipy.io import savemat
from APTCB_utils import ConfirmationBiasSimulator, ConfirmationBias,ConfirmationBiasStats
from APTCB_utils import makeFigures
#%%
# load the MATLAB engine
if __name__ == "__main__":
    
    """"
    variables to change for different inference tasks:
    filename
    fieldstofit and corresponding range
    """
    # create a log file
    time_stamp = datetime.datetime.now()
    currentTime  = time_stamp.strftime('%Y.%m.%d-%H:%M')
    logname = 'inferenceAPT' + currentTime
    logging.basicConfig(filename = logname,level = logging.DEBUG)
    # create a savefolder
    savefolder = '../dscData/resultsSyntheticDataCB_debug/'
    if not os.path.exists(savefolder):  
        os.mkdir(savefolder)
    
    # set some hyperparameters
    # training schedule
    n_train = 10000
    n_rounds = 4
    # fitting setup

    val_frac = 0.05
    # network setup
    
    n_hiddens = [50,50]

    # MAF parameters
    density = 'mog'
    n_mades = 5         # number of MADES
    
   
    
    seed_inf = 1
    pilot_samples = 2000


    
    # network setup
    
    minibatch = 500
    epochs = 100
    
    prior_norm = True
    

    hyperparameters = dict()
    hyperparameters['n_train'] = n_train
    hyperparameters['n_rounds'] = n_rounds
    hyperparameters['val_frac'] = val_frac
    hyperparameters['n_hiddens'] = n_hiddens
    hyperparameters['density'] = density
    hyperparameters['n_mades'] = n_mades
    
    # start matlab enigne and set path
    engine = matlab.engine.start_matlab()
    try:
        ps = engine.genpath('/scratch/sliu88/APT')
        engine.addpath(ps)
    except:
        pass
    # save inference method settings (hyperparameters)
    savemat(os.path.join(savefolder,'hyperparameters.mat'),{'hyperparameters':hyperparameters})
    C = 2
    g = 2
    l = 2
    s = 2
    
                   
    # load data
    datafolder = '../dscData/syntheticDataCB'
    filename = 'Trial1priorC%dgamma%dlapse%dsamples%d' %(C,g,l,s)
    dataname = os.path.join(datafolder,filename+'.mat')
    savename = os.path.join(savefolder,'lshc_res'+ filename+'.mat')
    savenamepkl = os.path.join(savefolder,'lshc_res'+ filename+'.pkl')
    if not os.path.isfile(savename):
        logging.info('Running %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))
        # load mat
        data = loadmat(dataname)
            
        lshc_choices = np.asarray(data['lshc_sim_results']['choices'])
        hslc_choices = np.asarray(data['hslc_sim_results']['choices'])
        lshc_signals = np.asarray(data['lshc_signals'])
        hslc_signals = np.asarray(data['hslc_signals'])
 
        lshc_params = data['lshc_params']
        hslc_params = data['hslc_params']
     
        lshc_infoType = np.ones(lshc_choices.shape)
        hslc_infoType = np.ones(hslc_choices.shape) * -1
        #%%
        # summary statistics of lshc condition
        lshc_s = ConfirmationBiasStats(lshc_signals,lshc_infoType,condition = 'choices',useStd = False)
        lshc_obs = {'choices':lshc_choices} 
        lshc_obs_stats = lshc_s.calc([lshc_obs])
        # summary statistics of hslc condition
        hslc_s = ConfirmationBiasStats(hslc_signals,hslc_infoType,condition = 'choices',useStd = False)
        hslc_obs = {'choices':hslc_choices}
        hslc_obs_stats = hslc_s.calc([hslc_obs])
        
    
        #%%
        
        
        fieldstofit = ['prior_C','gamma','lapse']
                           
        prior_min = np.array([0,0,0])
        prior_max = np.array([1,0.5,0.5])
       
        seed_p = 2
        prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
        #%%
   
     
        lshc_m = ConfirmationBias(lshc_signals,fieldstofit,lshc_params,engine)
        
        
        gt = dg.Default(model=lshc_m, prior=prior, summary=lshc_s)


   
   

        # inference object
        res = infer.SNPEC(gt,
                        obs=lshc_obs_stats,
                        n_hiddens=n_hiddens,
                        seed=seed_inf,
                        pilot_samples=pilot_samples,
                        n_mades=n_mades,
                        prior_norm=prior_norm,
                        density=density)
        
        # train
        log, _, posterior = res.run(
                            n_train=n_train,
                            n_rounds=n_rounds,
                            minibatch=minibatch,
                        epochs=epochs,
                        silent_fail=False,
                        proposal='gaussian',
                        val_frac=val_frac,
                        verbose=True,)
    
        # save posterior 
        with open(savenamepkl,'wb') as f:
            pickle.dump([log,posterior],f)
        # generate 2000 samples from posterior
        posteriorSamples = []
        for n in range(n_rounds):
            posteriorSamples.append(posterior[n].gen(2000))
        # save these samples to a .matfile
    
   
    
    
        savemat(savename,{'posterSamples':posteriorSamples,'log':log})
        #logging.info('Finished %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))
        
        
        
        paraInfo = loadmat(os.path.join('../dscData/SyntheticDataCB','genDataInfo.mat'))
        prior_C_list = paraInfo['genDataInfo']['prior_C_list'].tolist()[0][0][0]
        gamma_list = paraInfo['genDataInfo']['gamma_list'].tolist()[0][0][0]
        lapse_list = paraInfo['genDataInfo']['lapse_list'].tolist()[0][0][0]
        samples_list = paraInfo['genDataInfo']['samples_list'].tolist()[0][0][0]
        true_params = [prior_C_list[C-1],gamma_list[g-1],lapse_list[l-1],samples_list[s-1]]
        makeFigures(log,posterior,true_params,fieldstofit,prior_min,prior_max)
