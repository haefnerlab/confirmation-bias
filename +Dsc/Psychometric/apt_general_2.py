#%%
from delfi.simulator.BaseSimulator import BaseSimulator
import os
from delfi.summarystats.BaseSummaryStats import BaseSummaryStats
import matplotlib.pyplot as plt
import numpy as np
import delfi.distribution as dd
import delfi.generator as dg
import delfi.inference as infer
from delfi.utils.viz import samples_nd

"""

This file contains the Noise Model class for 2 parameters, 
such so, since with the `x`, I also provide the err term to
the linearNoiseModel

"""

class psychometricModel(BaseSimulator):
    
    def __init__(self, simulator, labels, true_params, dim_param, seed = None):

        super().__init__(dim_param=dim_param, seed=seed)
        
        self.simulate = simulator
        self.labels = labels
        self.true_params = true_params

    def gen_single(self, params):

        params = np.asarray(params)

        assert params.ndim == 1, 'params.ndim must be 1'

        y = self.simulate(params)

        # return {'y': y.reshape(-1)}
        # print(f'y is {y}')
        # print(f'Along 0 axis: {np.sum(y, 1)}')
        # print(f'Along 1 axis: {np.sum(y, 0)}')
        # print(f'reshaped: {y.reshape(-1)}')
        # print(f'last: {np.sum(y.reshape(-1).reshape(-1, 50), 1)}')
        # print(f'last2: {np.sum(y.reshape(-1).reshape(-1, 50), 0)}')
        # 21*50 by 1
        return {'y': y.reshape(-1)}  
    
class psychometricStats(BaseSummaryStats):
    
    def __init__(self, N): # Number of trial per stimulus condition):
        super().__init__()
        self.N = N
    def calc(self,repetition_list):
        
        stats = []
        for r in range(len(repetition_list)):
            data = repetition_list[r]
            # make sure is one dimensional 
            # reshape back 21 by 50
            # sum over 50
            observed_choices_avg = np.sum(data['y'].reshape(-1, self.N), 1)
            
            # print(observed_choices_avg)
            stats.append(observed_choices_avg)
        return stats
    
class psychometricStats2(BaseSummaryStats):
    
    def __init__(self,signals,N):
        self.signals = signals
        self.N = N
    def calc(self,repetition_list):
        stats = []
        for r in range(len(repetition_list)):
            data = repetition_list[r]
            choices = np.asarray(data['y'].reshape(-1,N))
            signals = np.asarray(self.signals)

            # group trials according to just choices (2 groups)
            ind1 = np.where(choices == 1)[0]
            ind2 = np.where(choices == 0)[0]

            # include standard error mean in summary statistics
            sum_stats_vec = np.hstack((np.mean(signals[ind1],axis = 0),np.mean(signals[ind2],axis = 0),\
                                       np.std(signals[ind1],axis = 0)/np.sqrt(len(ind1)),\
                                       np.std(signals[ind2],axis = 0)/np.sqrt(len(ind2)))) 
            stats.append(sum_stats_vec)
        return stats


def plot_APT(posterior, g, labels, true_params, fignames):
    
    posterior_samples = posterior[0].gen(1000)

    prior_min = g.prior.lower
    prior_max = g.prior.upper
    
    prior_lims = np.concatenate((prior_min.reshape(-1,1),prior_max.reshape(-1,1)),axis=1)

    ###################
    # colors
    hex2rgb = lambda h: tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

    # RGB colors in [0, 255]
    col = {}
    col['GT']      = hex2rgb('30C05D')
    col['SNPE']    = hex2rgb('2E7FE8')
    col['SAMPLE1'] = hex2rgb('8D62BC')
    col['SAMPLE2'] = hex2rgb('AF99EF')

    # convert to RGB colors in [0, 1]
    for k, v in col.items():
        col[k] = tuple([i/255 for i in v])

    ###################
    # posterior
    fig, axes = samples_nd(posterior_samples,
                        limits=prior_lims,
                        ticks=prior_lims,
                        labels = labels,
                        fig_size=(5,5),
                        diag='kde',
                        upper='kde',
                        hist_diag={'bins': 50},
                        hist_offdiag={'bins': 50},
                        kde_diag={'bins': 50, 'color': col['SNPE']},
                        kde_offdiag={'bins': 50},
                        points=[true_params],
                        points_offdiag={'markersize': 5},
                        points_colors=[col['GT']],
                        title='');
    
    plt.savefig(fignames + '_plot.png')
    # plt.show()
    
def runAPT2psychometric(obs0, hyps, labels, true_params, fignames, plot = True):
    """
    obs0:
        The data observation we made
        
    This takes in a dict with the following parameters:
    
    prior       : prior estimate and bounds
    m           : The model we are going to run apt on
    s           : The summary statistics we care about
    n_train     : Number of param-samples drawn during a training round
    n_rounds    : Number of rounds you want to run APT for
    n_hidden    : A Number of Layers x 1 array containing the layer sizes
    pilot_samples : Number of samples drawn during the pilot run
    val_frac    : Fraction of validation samples 
    minibatch   : Size of a batch
    epochs      : Number of training iterations
    density     : 'mog' or 'maf' supported
    n_mades     : number of Mades, an MAF parameter
         
    """
    
    g = dg.Default(model=hyps['m'], prior=hyps['prior'], summary=hyps['s'])

    # define statsitics summary of observation
    obs_stats = hyps['s'].calc([obs0])

    seed_inf = hyps['seed_inf']

    prior_norm = True

    # MAF parameters
    density = 'maf'
    n_mades = 5         # number of MADES

    # inference object
    res = infer.SNPEC(g,
                    obs = obs_stats,
                    n_hiddens = hyps['n_hiddens'],
                    seed = seed_inf,
                    pilot_samples = hyps['pilot_samples'],
                    n_mades= hyps['n_mades'],
                    prior_norm= prior_norm,
                    density = hyps['density'])
    # train
    log, _, posterior = res.run(
                        n_train = hyps['n_train'],
                        n_rounds = hyps['n_rounds'],
                        minibatch = hyps['minibatch'],
                        epochs = hyps['epochs'],
                        silent_fail = False,
                        proposal = 'prior',
                        val_frac = hyps['val_frac'],
                        verbose = True,)
    
    if plot == True:
        plot_APT(posterior, g, labels, true_params, fignames)
    
    return posterior

# %%
