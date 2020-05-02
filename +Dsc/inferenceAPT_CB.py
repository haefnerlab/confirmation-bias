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
plt.rcParams.update({'font.size': 16})
from delfi.utils.viz import samples_nd

#%%
def ConfirmationBiasSimulator(xval,fieldstofit,params, signals, engine):
    # run matlab code run.vectorized and read results
   

    paramsNew = copy.deepcopy(params)
    for i,f in enumerate(fieldstofit):
        if type(params[f]) is float:
            paramsNew[f] = float(xval[i])
        elif type(params[f]) is int:
            paramsNew[f] = int(xval[i])


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
                    sum_stats_vec = np.hstack((np.mean(signals[ind1,:],axis = 0),np.mean(signals[ind2,:],axis = 0)))   
            elif self.condition == 'infoType':
                # group trials according to just category information (2 groups)
                ind1 = np.where(infoType == 1)[0]
                ind2 = np.where(infoType == -1)[0]
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
    C = 1 # prior_C =  0.5
    g = 2 # gamma = 0
    l = 2 # lapse = 0.1 
    s = 1 # sample = 5
    datafolder = '../dscData/syntheticDataCB'
    filename = 'Trial1priorC%dgamma%dlapse%dsamples%d.mat' %(C,g,l,s)
    dataname = os.path.join(datafolder,filename)
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
    lshc_s = ConfirmationBiasStats(lshc_signals,lshc_infoType,condition = 'choices',useStd = True)
    lshc_obs = {'choices':lshc_choices}
    
    lshc_obs_stats = lshc_s.calc([lshc_obs])
    
    # plt.subplot(121)
    # plt.errorbar(np.arange(10),lshc_obs_stats[0][0:10],lshc_obs_stats[0][20:30],color = 'red',label = 'choice 1')
    # plt.errorbar(np.arange(10),lshc_obs_stats[0][10:20],lshc_obs_stats[0][30:40],ls='--',color = 'red',label = 'choice -1')
    # plt.ylim([-2,2])
    # plt.xlabel('Frame')
    # plt.ylabel('Signal')
    # plt.title('LSHC')
    # plt.legend()
    # plt.tight_layout()
    # summary statistics of hslc condition
    hslc_s = ConfirmationBiasStats(hslc_signals,hslc_infoType,condition = 'choices',useStd = True)
    hslc_obs = {'choices':hslc_choices}
    
    hslc_obs_stats = hslc_s.calc([hslc_obs])
    
    # plt.subplot(122)
    # plt.errorbar(np.arange(10),hslc_obs_stats[0][0:10],hslc_obs_stats[0][20:30],color = 'blue',label = 'choice 1')
    # plt.errorbar(np.arange(10),hslc_obs_stats[0][10:20],hslc_obs_stats[0][30:40],ls='--',color = 'blue',label = 'choice -1')
    # plt.ylim([-2,2])
    # plt.xlabel('Frame')
    # plt.ylabel('Signal')
    # plt.title('HSLC')
    # plt.legend()
    # plt.tight_layout()
    # plt.suptitle('priorC = %.1f, gamma = %.1f, lapse = %.1f, samples = %d' %(hslc_params['prior_C'],\
    #                                                                    hslc_params['gamma'],\
    #                                                                 hslc_params['lapse'],\
    #                                                                 hslc_params['samples']),y = 0.99)
        
    #%%
    engine = matlab.engine.start_matlab()
    
    fieldstofit = ['prior_C','gamma','lapse','samples']
                       
    prior_min = np.array([0,0,0,1])
    prior_max = np.array([1,1,1,100])
   
    seed_p = 2
    prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    lshc_m = ConfirmationBias(lshc_signals,fieldstofit,lshc_params,engine)
    
    g = dg.Default(model=lshc_m, prior=prior, summary=lshc_s)
    
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
    
    prior_norm = True
    
    # MAF parameters
    density = 'mog'
    n_mades = 5         # number of MADES
    



    # inference object
    res = infer.SNPEC(g,
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
        
    
    
    
    
    
    
    
    
    
    
    
    
    