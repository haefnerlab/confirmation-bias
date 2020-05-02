## Fitting method Comparison
## A billion imports
#%%
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
import corner

import pickle as pk
import os.path as path
#%%
def linearNoiseSimulator(x, err, params, seed=None):
    
    alpha = params[0]
    beta = params[1]
    N = x.shape[0]
    #np.random.seed(123)
    y = x * alpha + beta + err * np.random.randn(N)
    return y

class linearNoiseModel(BaseSimulator):
    def __init__(self, x, err, dim_param, seed = None):

        super().__init__(dim_param=dim_param, seed=seed)
        self.x = x
        self.simulate = linearNoiseSimulator
        self.err = err

    def gen_single(self, params):

        params = np.asarray(params)

        assert params.ndim == 1, 'params.ndim must be 1'

        hh_seed = np.random.seed(1234)

        y = self.simulate(self.x, self.err, params,seed=hh_seed)

        return {'y': y.reshape(-1)}    
class linearNoiseStats(BaseSummaryStats):
    def __init__(self, x):
        self.x = x
        
    def calc(self,repetition_list):
        seg = 100
        stats = []

        for r in range(len(repetition_list)):
            # data = repetition_list[r]
            # x = np.asarray(data['x'])
            # y = np.asarray(data['y'])

            # sum_stats_vec = np.concatenate((x,y))
            # stats.append(sum_stats_vec)
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
                sum_stats_vec.extend([x_sorted[ind[0]], x_sorted[ind[-1]], y_min,y_max])
            stats.append(sum_stats_vec)
        return stats   
class linearNoiseStatsAlternate(BaseSummaryStats):
    
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
#%%
def generateLinearToyData(params, err, N):
    alpha = params[0]
    beta = params[1]
    
    np.random.seed(1234)
    
    x0 = np.random.uniform(0, 10, N)
    
    m = linearNoiseModel(x0, err, dim_param = 2)
    obs0 = m.gen_single(true_params)
    
    prior_min = np.array([-10,-10])
    prior_max = np.array([10,10])
    seed_p = 2
    
    prior =  dd.Uniform(lower = prior_min , upper = prior_max,seed = seed_p)
    
    return x0, obs0, m, prior

def runAPT(m, prior, obs0, x0, statsFunction):
    
    s = statsFunction(x0)
    g = dg.Default(model=m, prior=prior, summary=s)

    # define statsitics summary of observation
    obs_stats = s.calc([obs0])
    # training schedule
    n_train = 3000
    n_rounds = 2
    seed_inf = 1
    
    pilot_samples = 2000

    val_frac = 0.05
    # network setup
    n_hiddens = [200, 50,50]
    minibatch = 500
    epochs = 100

    prior_norm = True

    # MAF parameters
    density = 'maf'
    n_mades = 5         # number of MADES

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
                    proposal='prior',
                    val_frac=val_frac,
                    verbose=True,)
    
    return g, s, posterior
#%%
## Analytical
def posterior_calc_lineFitting(a, b, x, y):
    
    prior = {'mu_alpha': 10, 'var_alpha': 2, 'mu_beta': 2.5, 'var_beta': 10} 

    dy = y - a * x - b
    
    ivar = 1 / err **2
    lnlikeli =  -0.5 * (N* np.log(2*np.pi) + np.sum(2*np.log(err)) + np.sum(dy**2 * ivar))
    lnprior_a = -.5 * (a - prior['mu_alpha'])**2 / prior['var_alpha']**2
    lnprior_b = -.5 * (b - prior['mu_beta'])**2  / prior['var_beta']**2
    # lnprior_c = -.5 * (c - prior['mu_sigma']**2) / prior['var_sigma']**2
    lnposterior = lnlikeli + lnprior_a + lnprior_b
    if np.isnan(lnposterior):
        lnposterior = -np.inf
    if np.isnan(lnlikeli):
        lnlikeli = -np.inf
    return lnposterior

def runAnalytical(alpha_true_mean, beta_true_mean, x, y):
    
    alpha = np.linspace(-.1+alpha_true_mean, .1+alpha_true_mean, 400)
    beta = np.linspace(-.2+beta_true_mean, .2+beta_true_mean, 400)
    posterior = np.empty([np.size(beta), np.size(alpha)])
    
    for ib, b in enumerate(beta):
        for ia, a in enumerate(alpha):
                posterior[ib][ia] = posterior_calc_lineFitting(a, b, x, y)
                    
    posterior_1 = posterior - posterior.max()
    posterior_1[posterior_1 < -1000] = -1000
    analytic_dict = {'alpha_range': alpha, 'beta_range': beta, 'posterior': posterior_1}
    return analytic_dict, posterior_1

def sample_proposal(*sigmas):
    return np.random.normal(0., sigmas)
def run_metropolis_hastings(p0, x, y, n_steps, proposal_sigmas):
    """
    Run a Metropolis-Hastings MCMC sampler to generate samples from the input
    log-posterior function, starting from some initial parameter vector.
    
    Parameters
    ----------
    p0 : Initial parameter vector.
    n_steps : Number of steps to run the sampler for.
    model : function to calculate the posterior
    proposal_sigmas : list, array
        A list of standard-deviations passed to the sample_proposal 
        function. These are like step sizes in each of the parameters.
    """
    # p0 = np.array(p0)
    if len(proposal_sigmas) != len(p0):
        raise ValueError("Proposal distribution should have same shape as parameter vector.")
    
    # the objects we'll fill and return:
    chain = np.zeros((n_steps, len(p0))) # parameter values at each step
    ln_probs = np.zeros(n_steps) # log-probability values at each step
    
    # we'll keep track of how many steps we accept to compute the acceptance fraction                        
    n_accept = 0 
    
    # evaluate the log-posterior at the initial position and store starting position in chain
    ln_probs[0] = posterior_calc_lineFitting(p0[0], p0[1], x, y)
    chain[0] = p0
    
    # loop through the number of steps requested and run MCMC
    for i in range(1,n_steps):
        # proposed new parameters
        step = sample_proposal(*proposal_sigmas)
        new_p = chain[i-1] + step
        
        # compute log-posterior at new parameter values
        new_ln_prob = posterior_calc_lineFitting(new_p[0], new_p[1], x, y)
        
        # log of the ratio of the new log-posterior to the previous log-posterior value
        ln_prob_ratio = new_ln_prob - ln_probs[i-1]
        
        if (ln_prob_ratio > 0) or (ln_prob_ratio > np.log(np.random.uniform())):
            chain[i] = new_p
            ln_probs[i] = new_ln_prob
            n_accept += 1
            
        else:
            chain[i] = chain[i-1]
            ln_probs[i] = ln_probs[i-1]
    
    acc_frac = n_accept / n_steps
    return chain, ln_probs, acc_frac

def plot_analytical_line_fitting(posterior_1, good_samples, alpha, beta, alpha_true_mean, beta_true_mean):
    fig,ax = plt.subplots(1, 1, figsize=(5,5))
    a_i = 0
    a_e = 400
    b_i =  0
    b_e = 400

    assert a_i < a_e
    assert b_i < b_e 

    posterior_2 = posterior_1[np.ix_(np.arange(b_i, b_e), np.arange(a_i, a_e) )]
    ax.contourf(alpha[a_i:a_e], beta[b_i:b_e], posterior_2, cmap='Blues', levels=100, vmin=posterior_2.max()-128, vmax=posterior_2.max())
    
    ax.plot(alpha_true_mean, beta_true_mean, marker = 'o', zorder=999, color='r')

    return ax
#%%
def open_file(name):
    if path.exists(name):
        with open(name, 'rb') as handle:
            file = pk.load(handle)
    else:
        return -1
    return file
def save_file(name, file):
    with open(name, 'wb') as handle:
        pk.dump(file, handle, protocol=pk.HIGHEST_PROTOCOL)
        
def plot_apt(g, posterior_apt):
    posterior_samples = posterior_apt[0].gen(1000)

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
                        labels=labels_params,
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
#%%
    
N = 1000
alpha = .5
beta = -2
err = 0.4
true_params = [alpha,beta]
labels_params = ['alpha','beta']

file_paths = {
    'init_data': 'LineFitting/init_data.pk',
    'analytical': 'LineFitting/analytical.pk',
    'mh_samples': 'LineFitting/mh_samples.pk',
    'posterior_apt_dict': 'LineFitting/posterior_apt_dict.pk',
}

#%%
## Generate data
init_data = open_file(file_paths['init_data'])
if init_data == -1:
    x0, obs0, m, prior = generateLinearToyData(true_params, err, N)
    init_data = {'x0': x0, 'obs0': obs0, 'm': m, 'prior': prior}
    save_file(file_paths['init_data'], init_data)
else:
    x0 = init_data['x0']
    obs0 = init_data['obs0']
    m = init_data['m']
    prior = init_data['prior']
#%%
analytical = open_file(file_paths['analytical'])
if analytical == -1:
    analytic_dict, posterior_analytical = runAnalytical(true_params[0], true_params[1], x0, obs0['y'])
    save_file(file_paths['analytical'], analytic_dict)
else:
    posterior_analytical = analytical['posterior']
    analytic_dict = analytical
#%%
mh_samples = open_file(file_paths['mh_samples'])
if mh_samples == -1:
    p0 = [.47, -1.81]
    chain,_,acc_frac = run_metropolis_hastings(p0, x0, obs0['y'], n_steps=180000, proposal_sigmas=[0.05,0.05])
    print("Acceptance fraction: {:.1%}".format(acc_frac))
    mh_samples = chain[2000::8]
    save_file(file_paths['mh_samples'], mh_samples)
#%%
ax = plot_analytical_line_fitting(posterior_analytical,mh_samples,
    analytic_dict['alpha_range'], analytic_dict['beta_range'], true_params[0], true_params[1])
corner.hist2d(mh_samples[:, 0], mh_samples[:, 1], ax = ax,zorder=5, fill_contours=False)
#%%

posterior_apt_dict = open_file(file_paths['posterior_apt_dict'])

if posterior_apt_dict == -1:
    
    g, s, posterior_apt = runAPT(m, prior, obs0, x0, linearNoiseStats)
    posterior_apt_dict = {'g': g, 's': s, 'posterior_apt': posterior_apt}
    save_file(file_paths['posterior_apt_dict'], posterior_apt_dict)
else:
    g = posterior_apt_dict['g']
    s = posterior_apt_dict['s']
    posterior_apt = posterior_apt_dict['posterior_apt']
  
plot_apt(g, posterior_apt)
# %%
ranger = [[analytic_dict['alpha_range'].min(), analytic_dict['alpha_range'].max()], 
[analytic_dict['beta_range'].min(),  analytic_dict['beta_range'].max()]]
#%%
# apt_samples_r = np.reshape(apt_samples, apt_samples.size, order='F').reshape(apt_samples.shape[1], apt_samples.shape[0])
ax = plot_analytical_line_fitting(posterior_analytical,mh_samples,
analytic_dict['alpha_range'], analytic_dict['beta_range'], true_params[0], true_params[1])
corner.hist2d(apt_samples_r[0], apt_samples_r[1], ax = ax,zorder=5, fill_contours=False, 
              range = ranger)

#%%
ax = plot_analytical_line_fitting(posterior_analytical,mh_samples,
analytic_dict['alpha_range'], analytic_dict['beta_range'], true_params[0], true_params[1])
corner.hist2d(apt_samples_r[0], apt_samples_r[1], ax = ax,zorder=5, fill_contours=False, 
              range = None)



# %%
