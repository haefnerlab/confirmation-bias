#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  2 13:55:48 2020

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
#%%
def ConfirmationBiasSimulator(xval,fieldstofit,params, signals, engine):
    # run matlab code run.vectorized and read results
   

    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        if f == 'samples':
            paramsNew[f] = int(xval[i])
        elif type(params[f]) is int:
            paramsNew[f] = float(xval[i])


    sim_results = engine.Model.runVectorizedAPT(paramsNew, signals.tolist())
    sim_choice = np.squeeze(np.asarray(sim_results['choices']))



    return sim_choice

# define a simulator class linking to the real simulator
class ConfirmationBias(BaseSimulator):
    
    def __init__(self, signals,fieldstofit,params,engine, seed = None):
        
        dim_param = 4
        super().__init__(dim_param = dim_param,seed = seed)
        self.signals = signals
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        choices = self.simulate(xval,self.fieldstofit,self.params,self.signals, self.engine)
        return {'choices':choices}
    
class ConfirmationBiasStats(BaseSummaryStats):
    
    def __init__(self,signals,infoType,condition = 'choices_infoType',useStd = True):
        self.signals = signals
        self.infoType = infoType
        self.condition = condition
        self.useStd = useStd
    def calc(self,repetition_list):

        stats = []
    
        for r in range(len(repetition_list)):

            data = repetition_list[r]
            choices = np.asarray(data['choices'])
            signals = np.asarray(self.signals)
            infoType = np.asarray(self.infoType)
            
            if self.condition == 'choices_infoType':
                # group trials according to both choices and catergory information (4 groups)
                ind1 = np.where( (choices == 1) & (infoType == 1))[0]
                ind2 = np.where( (choices == 1) & (infoType == -1))[0]
                ind3 = np.where( (choices == -1) & (infoType == 1))[0]
                ind4 = np.where( (choices == -1) & (infoType == -1))[0]
                if self.useStd: 
                    # include standard error mean in summary statistics
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.mean(signals[ind3,:],axis = 0),np.mean(signals[ind4,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)),\
                                               np.std(signals[ind3,:],axis = 0)/np.sqrt(len(ind3)),\
                                               np.std(signals[ind4,:],axis = 0)/np.sqrt(len(ind4)))) 
                else:
                    # only use average of signals as summary statistics 
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                   np.mean(signals[ind3,:],axis = 0),np.mean(signals[ind4,:],axis = 0)))    
            elif self.condition == 'choices':
                # group trials according to just choices (2 groups)
                ind1 = np.where(choices == 1)[0]
                ind2 = np.where(choices == -1)[0]
                if self.useStd:
                    # include standard error mean in summary statistics
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    # only use average of signals as summary statistics 
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                              np.mean(np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1))),\
                                                  np.mean(np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))))
            elif self.condition == 'infoType':
                # group trials according to just category information (2 groups)
                ind1 = np.where(infoType == 1)[0]
                ind2 = np.where(infoType == -1)[0]
                if self.useStd:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                              np.mean(np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1))),\
                                                  np.mean(np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))))
                    
                    
            stats.append(sum_stats_vec)
        return stats
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
    savefolder = '../dscData/resultsSyntheticDataCB_version2/'
    if not os.path.exists(savefolder):  
        os.mkdir(savefolder)
    
    # set some hyperparameters
    # training schedule
    n_train = 5000
    n_rounds = 3
    # fitting setup

    val_frac = 0.05
    # network setup
    
    n_hiddens = [50,50]

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
    
    # start matlab enigne and set path
    engine = matlab.engine.start_matlab()
    try:
        ps = engine.genpath('/scratch/sliu88/APT')
        engine.addpath(ps)
    except:
        pass
    # save inference method settings (hyperparameters)
    savemat(os.path.join(savefolder,'hyperparameters.mat'),{'hyperparameters':hyperparameters})
    priorC_list = [1,2,3]
    gamma_list = [1,2,3]
    lapse_list = [1,2,3]
    samples_list = [1,2,3]
    for s in samples_list:
        for l in lapse_list:
            for g in gamma_list:
                for C in priorC_list:
                   
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
                        
                        
                        fieldstofit = ['prior_C','gamma','lapse','samples']
                                           
                        prior_min = np.array([0,0,0,1])
                        prior_max = np.array([1,0.5,0.5,50])
                       
                        seed_p = 2
                        prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
                        #%%
                   
                     
                        lshc_m = ConfirmationBias(lshc_signals,fieldstofit,lshc_params,engine)
                        
                        
                        gt = dg.Default(model=lshc_m, prior=prior, summary=lshc_s)
                        
                        # set hyparameters 
                        # training schedule
                        n_train = 5000
                        n_rounds = 2
                        seed_inf = 1
                        pilot_samples = 2000
                    
                    
                        val_frac = 0.05
                        # network setup
                        n_hiddens = [50,50]
                        minibatch = 500
                        epochs = 100
                        
                        prior_norm = False
                        
                        # MAF parameters
                        density = 'mog'
                        n_mades = 5         # number of MADES
    
   
   
            
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
