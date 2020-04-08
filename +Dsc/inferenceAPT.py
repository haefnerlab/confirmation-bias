#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 15:26:46 2020

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
    
    def __init__(self, dataPath, fieldstofit,params,engine, seed = None):
        
        dim_param = 5
        super().__init__(dim_param = dim_param,seed = seed)
        self.dataPath = dataPath
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        sim_results = self.simulate(xval,self.fieldstofit,self.params,self.dataPath, self.engine)
        return sim_results
#%%

class ConfirmationBiasStats(BaseSummaryStats):
    def __init__(self):
        self.n_summary = 1
    def calc(self,choice):
        stats = choice

     
        return stats
#%%
def runAPTinference(dataname,fieldstofit,prior):
    # load the data

    data = loadmat(dataname)
    choice = np.squeeze(np.asarray(data['sim_results']['choices']))
    params = data['params']

    
    # inference parameters
    seed_inf = 1
    pilot_samples = 2000
    # training schedule
    n_train = 2000
    n_rounds = 3
    # fitting setup
    minibatch = 100
    epochs = 100
    val_frac = 0.05
    # network setup
    
    n_hiddens = [50,50]
    # convenience
    prior_norm = True
    # MAF parameters
    density = 'maf'
    n_mades = 5         # number of MADES
    
    
    # define generator class
    s = ConfirmationBiasStats()
   
    n_processes = 2
    seeds_m = np.arange(1,n_processes+1,1)
    
    #
    engine = matlab.engine.start_matlab()
    
    ps = engine.genpath('/gpfs/fs2/scratch/sliu88/APT')
    engine.addpath(ps)
    m = ConfirmationBias(dataname, fieldstofit,params,engine, seeds_m[0])
    g = dg.Default(model=m, prior=prior, summary=s)
    

    # do inference
    # inference object
    res = infer.SNPEC(g,
                    obs=choice,
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
                    proposal='prior',
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
    savefolder = '../dscData/resultsSyntheticDataCB'+ currentTime
    os.mkdir(savefolder)
    
    priorC_list = [1,2,3]
    gamma_list = [1,2,3]
    lapse_list = [1,2,3]
    samples_list = [1,2,3]
    for C in priorC_list:
        for g in gamma_list:
            for l in lapse_list:
                for s in samples_list:
                   
                    # load data
                    datafolder = '../dscData/syntheticDataCB'
                    filename = 'Trial1priorC%dgamma%dlapse%dsamples%d.mat' %(C,g,l,s)
                    dataname = os.path.join(datafolder,filename)
                    
                    
                    # define prior over model parameters
                    fieldstofit = ['prior_C','gamma','samples','lapse']
                   
                    prior_min = np.array([0,0,1,0])
                    prior_max = np.array([1,1,100,1])
                   
                    seed_p = 2
                    prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
                
                
                    logging.info('Running %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))
                    log,posterior = runAPTinference(dataname,fieldstofit,prior)
                    
                    #%%
                    # generate 1000 samples from posterior
                    posteriorSamples = posterior[0].gen(1000)
                    # save these samples to a .matfile
                    from scipy.io import savemat
                   
                    
                    savename ='res'+ filename
                    savemat(os.path.join(savefolder,savename),{'posterSamples':posteriorSamples,'log':log})
                    logging.info('Finished %dpriorC %dgamma %dlapse %dsamples' %(C,g,l,s))