#%%
from utilities import open_file, save_file
import apt_general_2
import mh_sampler
import analyticLineFittingToy
import numpy as np
import matplotlib.pyplot as plt
import corner
import delfi.distribution as dd

#%%
if __name__ == '__main__':
    
    file_paths = {
        'init_data': 'init_data.pk',
        'analytical': 'analytical.pk',
        'mh_samples': 'mh_samples.pk',
        'posterior_apt': 'posterior_apt.pk',
    }
    
    what_to_do = {
        'init_data': False,
        'analytical': False,
        'mh_samples': False,
        'apt': True
    }
    
    N = 1000
    alpha_true = .5
    beta_true = -2
    err = 0.4
    
    a_i, a_e, b_i, b_e = 0, 400, 0, 400

    true_params = [alpha_true, beta_true]
    labels = ['slope', 'intercept']
    
    alpha = np.linspace(-.1+alpha_true, .1+alpha_true, 400)
    beta = np.linspace(-.2+beta_true, .2+beta_true, 400)
    
    ranger = [[alpha[0], alpha[-1]], 
              [beta[0],  beta[-1]]]
    
    np.random.seed(1234)
    
    x0 = np.random.uniform(0, 10, N)
    
    prior_analytical = {'mu_alpha': 10, 'var_alpha': 2, 'mu_beta': 2.5, 'var_beta': 10}
    
    analyticModel = analyticLineFittingToy.analyticalLinearNoise(x0, err, true_params, prior_analytical)
    
    obs0 = {'y': analyticModel.y.reshape(-1)} 
    
    s = apt_general_2.linearNoiseStats(x0)
    
    apt_model = apt_general_2.linearNoiseModel(analyticModel.linearNoiseSimulator, labels, true_params, dim_param = 2)
    
    seed_p = 2
    prior_apt =  dd.Uniform(lower = np.asarray([-1, -3]), upper = np.asarray([1, 1]), seed = seed_p)
    
    hyps = {'prior': prior_apt,
            'm': apt_model,
            's': s,
            'n_train': 4000,
            'n_rounds': 2,
            'n_hiddens': [200, 50, 50],
            'pilot_samples': 2000,
            'val_frac': 0.05,
            'minibatch': 500,
            'epochs': 100,
            'density': 'maf',
            'n_mades': 5 }
    
    # apt_general_2.runAPT2LinearNoise(obs0, hyps)

    posterior_1 = open_file(file_paths['analytical'])
    
    if posterior_1 is -1 or what_to_do['analytical']:
        
        posterior = np.empty([np.size(beta), np.size(alpha)])
        for ib, b in enumerate(beta):
            for ia, a in enumerate(alpha):
                    t = [a, b]
                    posterior[ib][ia] = analyticModel.posteriorLinearNoise(t)
                
        posterior_1 = posterior - posterior.max()
        posterior_1[posterior_1 < -1000] = -1000
        
        save_file(file_paths['analytical'], posterior_1)
        
    mh_samples = open_file(file_paths['mh_samples'])
    
    if mh_samples is -1 or what_to_do['mh_samples']:
        
        p0 = [.47, -1.81]
        
        chain, _ ,acc_frac = mh_sampler.run_metropolis_hastings(p0, analyticModel, n_steps=180000, proposal_sigmas=[0.05,0.05])
        print("Acceptance fraction: {:.1%}".format(acc_frac))
        
        mh_samples = chain[2000::8]

        save_file(file_paths['mh_samples'], mh_samples)
        
    apt_posterior = open_file(file_paths['posterior_apt'])
    if apt_posterior is -1 or what_to_do['apt']:
        apt_posterior = apt_general_2.runAPT2LinearNoise(obs0, hyps, plot = True)
        save_file(file_paths['posterior_apt'], apt_posterior)
    
    ## Plotting the Metropolis on top of the analytical posterior
    fig, ax = plt.subplots(2, 1, figsize=(5,5))
    
    posterior_2 = posterior_1[np.ix_(np.arange(b_i, b_e), np.arange(a_i, a_e) )]
    
    ax[0].contourf(alpha[a_i:a_e], beta[b_i:b_e], posterior_2, 
                   cmap='Blues', levels=100, vmin=posterior_2.max()-128, vmax=posterior_2.max())    
    ax[0].plot(true_params[0], true_params[1], marker = 'o', zorder=999, color='r')
    
    corner.hist2d(mh_samples[:, 0], mh_samples[:, 1], ax = ax[0],zorder=5, fill_contours=False, range = ranger)
    
    ax[1].contourf(alpha[a_i:a_e], beta[b_i:b_e], posterior_2, 
                   cmap='Blues', levels=100, vmin=posterior_2.max()-128, vmax=posterior_2.max())
    
    ax[1].plot(true_params[0], true_params[1], marker = 'o', zorder=999, color='r')
    

    apt_samples = apt_posterior[0].gen(1000)
    apt_samples_r = np.reshape(apt_samples, apt_samples.size, order='F').reshape(apt_samples.shape[1], apt_samples.shape[0])
    corner.hist2d(apt_samples_r[0], apt_samples_r[1], ax = ax[1], zorder=5, fill_contours=False, 
              range = ranger)
    
    plt.show()

# %%
