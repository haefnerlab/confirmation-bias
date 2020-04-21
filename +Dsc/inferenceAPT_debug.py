#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 09:39:03 2020

@author: liushizhao
"""

#%% 

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
#%%
def ConfirmationBiasSimulator(xval,fieldstofit,params, dataPath, engine):
    # run matlab code run.vectorized and read results
    # print("Hello I am in Simulator")

    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        if type(params[f]) is float:
            paramsNew[f] = float(xval[i])
        elif type(params[f]) is int:
            paramsNew[f] = int(xval[i])


    sim_results = engine.Model.runVectorizedAPT(paramsNew, dataPath)
    #print("Hello I am out of simulator")
    
    sim_choice = np.squeeze(np.asarray(sim_results['choices']))



    return sim_choice
# define a simulator class linking to the real simulator

class ConfirmationBias(BaseSimulator):
    
    def __init__(self, dataPath,signals,fieldstofit,params,engine, seed = None):
        
        dim_param = 4
        super().__init__(dim_param = dim_param,seed = seed)
        self.dataPath = dataPath
        self.signals = signals
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        choice = self.simulate(xval,self.fieldstofit,self.params,self.dataPath, self.engine)
        return {'choice':choice,'signals':self.signals}
#%%

class ConfirmationBiasStats(BaseSummaryStats):
    def __init__(self):
        self.n_summary = 2
    def calc(self,repetition_list):
        stats = []
        for r in range(len(repetition_list)):
            x = repetition_list[r]
            choice = np.asarray(x['choice'])
            #signals = np.asarray(x['signals']).flatten('C')
            
            #sum_stats_vec = np.concatenate((np.asarray(choice),np.asarray(signals)))
            
            stats.append(np.asarray(choice))
        return stats
#%%
def runAPTinference(dataname,fieldstofit,prior,hyperparameters):
    # load the data

    data = loadmat(dataname)
    choice = np.squeeze(np.asarray(data['sim_results']['choices']))
    params = data['params']
    signals = data['signals']
    
    # inference parameters
    seed_inf = 1
    pilot_samples = 2000
    # training schedule
    n_train = hyperparameters['n_train']
    n_rounds =  hyperparameters['n_rounds']
    # fitting setup
    minibatch = 500
    epochs = 100
    val_frac =  hyperparameters['val_frac']
    # network setup
    
    n_hiddens =  hyperparameters['n_hiddens']
    # convenience
    prior_norm = True
    # MAF parameters
    density =  hyperparameters['density']
    n_mades =  hyperparameters['n_mades']         # number of MADES
    
    
    # define generator class
    s = ConfirmationBiasStats()
   
    n_processes = 2
    seeds_m = np.arange(1,n_processes+1,1)
    
    #
    engine = matlab.engine.start_matlab()
    
    ps = engine.genpath('/gpfs/fs2/scratch/sliu88/APT')
    engine.addpath(ps)
    m = ConfirmationBias(dataname, signals,fieldstofit,params,engine, seeds_m[0])
    g = dg.Default(model=m, prior=prior, summary=s)
    
    obs = {'choice':choice,'signals':signals}
    obs_stats = s.calc([obs])
    # do inference
    # inference object
    res = infer.SNPEC(g,
                    obs=obs_stats,
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
                    proposal='mog',
                    val_frac=val_frac,
                    verbose=True,)
    return log,posterior
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
    n_train = 8000
    n_rounds = 3
    # fitting setup

    val_frac = 0.05
    # network setup
    
    n_hiddens = [200,50,50]

    # MAF parameters
    density = 'mog'
    n_mades = 5         # number of MADES
    
    hyperparameters = dict()
    hyperparameters['n_train'] = n_train
    hyperparameters['n_rounds'] = n_rounds
    hyperparameters['val_frac'] = val_frac
    hyperparameters['n_hiddens'] = n_hiddens
    hyperparameters['density'] = density
    hyperparameters['n_mades'] = n_mades
                    
    # save inference method settings (hyperparameters)
    savemat(os.path.join(savefolder,'hyperparameters.mat'),{'hyperparameters':hyperparameters})
    C = 1
    g = 2
    l = 2
    s = 2

   
    # load data
    datafolder = '../dscData/syntheticDataCB'
    filename = 'Trial1priorC%dgamma%dlapse%dsamples%d.mat' %(C,g,l,s)
    dataname = os.path.join(datafolder,filename)
    savename =os.path.join(savefolder,'res'+ filename)
    if not os.path.isfile(savename):
        # load mat
        data = loadmat(dataname)
        choice = np.squeeze(np.asarray(data['sim_results']['choices']))
      
        # define prior over model parameters
        fieldstofit = ['prior_C']
       
        prior_min = np.array([0])
        prior_max = np.array([1])
       
        seed_p = 2
        prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    #%%
       # print('Running %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))
        logging.info('Running %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))
        log,posterior = runAPTinference(dataname,fieldstofit,prior,hyperparameters)
        
        
        # generate 1000 samples from posterior
        posteriorSamples = posterior[0].gen(2000)
        # save these samples to a .matfile
        
       
        
        
        savemat(savename,{'posterSamples':posteriorSamples,'log':log})
        logging.info('Finished %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))

