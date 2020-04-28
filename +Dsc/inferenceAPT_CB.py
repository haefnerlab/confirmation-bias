#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 20:09:21 2020

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
import matplotlib.pyplot as plt
from delfi.utils.viz import samples_nd
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
    
    def __init__(self, dataPath,fieldstofit,params,engine, seed = None):
        
        dim_param = 4
        super().__init__(dim_param = dim_param,seed = seed)
        self.dataPath = dataPath
        self.simulate = ConfirmationBiasSimulator
        self.engine = engine
        self.fieldstofit = fieldstofit
        self.params = params
        
    def gen_single(self,xval):
        
        choices = self.simulate(xval,self.fieldstofit,self.params,self.dataPath, self.engine)
        return {'choices':choices}
    
class ConfirmationBiasStats(BaseSummaryStats):
    
    def __init__(self,signals,categories,condition = 'choices_categories',useStd = True):
        self.signals = signals
        self.categories = categories
        self.condition = condition
        self.useStd = useStd
    def calc(self,repetition_list):

        stats = []
    
        for r in range(len(repetition_list)):

            data = repetition_list[r]
            choices = np.asarray(data['choices'])
            signals = np.asarray(self.signals)
            categories = np.asarray(self.categories)
            
            if self.condition == 'choices_categories':
                # group trials according to both choices and catergory information (4 groups)
                ind1 = np.where( (choices == 1) & (categories == 1))[0]
                ind2 = np.where( (choices == 1) & (categories == -1))[0]
                ind3 = np.where( (choices == -1) & (categories == 1))[0]
                ind4 = np.where( (choices == -1) & (categories == -1))[0]
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
                ind1 = np.where(choices == 1)
                ind2 = np.where(choices == -1)
                if self.useStd:
                    # include standard error mean in summary statistics
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    # only use average of signals as summary statistics 
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0)))   
            elif self.condition == 'categories':
                # group trials according to just category information (2 groups)
                ind1 = np.where(categories == 1)
                ind2 = np.where(categories == -1)
                if self.useStd:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0),\
                                               np.std(signals[ind1,:],axis = 0)/np.sqrt(len(ind1)),\
                                               np.std(signals[ind2,:],axis = 0)/np.sqrt(len(ind2)))) 
                else:
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0))) 
                    
                    
            stats.append(sum_stats_vec)
        return stats
#%%
if __name__ == "__main__":
    C = 2 # prior_C =  0.5
    g = 1 # gamma = 0
    l = 2 # lapse = 0.1 
    s = 2 # sample = 5
    datafolder = '../dscData/syntheticDataCB'
    filename = 'Trial1priorC%dgamma%dlapse%dsamples%d.mat' %(C,g,l,s)
    dataname = os.path.join(datafolder,filename)
    data = loadmat(dataname)
    
    choices = np.asarray(data['sim_results']['choices'])
    signals = np.asarray(data['signals'])
    categories = np.asarray(data['categories'])
    params = data['params']
    
    s = ConfirmationBiasStats(signals,categories,condition = 'choices_categories',useStd = True)
    obs = {'choices':choices}
    
    obs_stats = s.calc([obs])
    
    plt.figure()
    plt.errorbar(np.arange(10),obs_stats[0][0:10],obs_stats[0][40:50])
    plt.errorbar(np.arange(10),obs_stats[0][10:20],obs_stats[0][50:60])
    plt.errorbar(np.arange(10),obs_stats[0][20:30],obs_stats[0][60:70])
    plt.errorbar(np.arange(10),obs_stats[0][30:40],obs_stats[0][70:80])
    
    
    # #%%
    # engine = matlab.engine.start_matlab()
    
    # fieldstofit = ['prior_C','gamma','lapse','samples']
                       
    # prior_min = np.array([0,0,0,1])
    # prior_max = np.array([1,1,1,100])
   
    # seed_p = 2
    # prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    # m = ConfirmationBias(dataname,fieldstofit,params,engine)
    
    # g = dg.Default(model=m, prior=prior, summary=s)
    # # set hyparameters 
    # # training schedule
    # n_train = 5000
    # n_rounds = 2
    # seed_inf = 1
    # pilot_samples = 2000


    # val_frac = 0.05
    # # network setup
    # n_hiddens = [50,50]
    # minibatch = 500
    # epochs = 100
    
    # prior_norm = True
    
    # # MAF parameters
    # density = 'mog'
    # n_mades = 5         # number of MADES
    



    # # inference object
    # res = infer.SNPEC(g,
    #                 obs=obs_stats,
    #                 n_hiddens=n_hiddens,
    #                 seed=seed_inf,
    #                 pilot_samples=pilot_samples,
    #                 n_mades=n_mades,
    #                 prior_norm=prior_norm,
    #                 density=density)
    
    # # train
    # log, _, posterior = res.run(
    #                     n_train=n_train,
    #                     n_rounds=n_rounds,
    #                     minibatch=minibatch,
    #                 epochs=epochs,
    #                 silent_fail=False,
    #                 proposal='gaussian',
    #                 val_frac=val_frac,
    #                 verbose=True,)
    
    
    # #%%
    # fig = plt.figure(figsize=(15,5))

    # plt.plot(log[0]['loss'],lw=2)
    # plt.xlabel('iteration')
    # plt.ylabel('loss');
    
    
    
    


    # prior_min = g.prior.lower
    # prior_max = g.prior.upper
    # prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)
    
    # posterior_samples = posterior[0].gen(10000)
    
    # ###################
    # # colors
    # hex2rgb = lambda h: tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
    
    # # RGB colors in [0, 255]
    # col = {}
    # col['GT']      = hex2rgb('30C05D')
    # col['SNPE']    = hex2rgb('2E7FE8')
    # col['SAMPLE1'] = hex2rgb('8D62BC')
    # col['SAMPLE2'] = hex2rgb('AF99EF')
    
    # # convert to RGB colors in [0, 1]
    # for k, v in col.items():
    #     col[k] = tuple([i/255 for i in v])
    
    # ###################
    # # posterior
    # fig, axes = samples_nd(posterior_samples,
    #                        limits=prior_lims,
    #                        ticks=prior_lims,
    #                        labels=labels_params,
    #                        fig_size=(5,5),
    #                        diag='kde',
    #                        upper='kde',
    #                        hist_diag={'bins': 50},
    #                        hist_offdiag={'bins': 50},
    #                        kde_diag={'bins': 50, 'color': col['SNPE']},
    #                        kde_offdiag={'bins': 50},
    #                        points=[true_params],
    #                        points_offdiag={'markersize': 5},
    #                        points_colors=[col['GT']],
    #                        title='');
        
    
    
    
    
    
    
    
    
    
    
    
    
    