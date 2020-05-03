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

class linearNoiseModel(BaseSimulator):
    
    def __init__(self, simulator, labels, true_params, dim_param, seed = None):

        super().__init__(dim_param=dim_param, seed=seed)
        
        self.simulate = simulator
        self.labels = labels
        self.true_params = true_params

    def gen_single(self, params):

        params = np.asarray(params)

        assert params.ndim == 1, 'params.ndim must be 1'

        y = self.simulate(params)

        return {'y': y.reshape(-1)}  
    
class linearNoiseStats(BaseSummaryStats):
    
    def __init__(self, x):
        self.x = x
        
    def calc(self,repetition_list):
        seg = 20
        stats = []

        for r in range(len(repetition_list)):

            data = repetition_list[r]
            x = self.x
            y = np.asarray(data['y'])
            N = x.shape[0]
            x_sorted = np.asarray([x for x,y in sorted(zip(x,y))])
            y_sorted = np.asarray([y for x,y in sorted(zip(x,y))])
            sum_stats_vec = []
            for i in np.arange(seg):
                ind = np.arange(i*int(N/seg),(i+1)*int(N/seg)).astype(int)
                y_min = np.min(y_sorted[ind])
                y_max = np.max(y_sorted[ind])
                sum_stats_vec.extend([y_min,y_max])
            stats.append(sum_stats_vec)
        return stats   
    
class linearNoiseStatsAlternate(BaseSummaryStats):
    """
    NOTE: This code has not been verified yet
    """
    def __init__(self, x):
        self.x = x
        
    def calc(self,repetition_list):
        stats = []
    
        for r in range(len(repetition_list)):
            
            data = repetition_list[r]
            y = np.asarray(data['y'])
            y_sorted = [yi for xi,yi in sorted(zip(self.x,y))]
            y_sorted = np.array(y_sorted)
            sum_stats_vec = []
            bins = np.linspace(y_sorted[0], y_sorted[-1], 50)
            digitized = np.digitize(y_sorted, bins, 'left')
            bin_means = [y_sorted[digitized == i].mean() for i in range(len(bins))]
            bin_vars = [y_sorted[digitized == i].var() for i in range(len(bins))]
            var = np.mean(bin_vars)
            stats.extend(bin_means)
            stats.extend(bin_vars)
            stats.append(var)
        return [stats]
    
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
    
def runAPT2LinearNoise(obs0, hyps, labels, true_params, fignames, plot = True):
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