# -*- coding: utf-8 -*-
"""Zenodo HDX+TrIQ sampling

# Prepare environment
"""


## Setup Colab environment
## 
## !pip install -q condacolab
## import condacolab
## condacolab.install()
## 
## from google.colab import drive
## drive.mount('/content/drive')
## 
## !pip install mdtraj
## !pip install openpyxl --upgrade --pre
## 
## condacolab.check()
## !conda install -y -q -c conda-forge openmm python=3.7 pdbfixer openmmtools pyemma 2>&1 1>/dev/null
## 

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.linalg import expm
from numpy.linalg.linalg import matrix_power

import mdtraj as md
import pyemma


## Sampling various thermodynamic models with ME forward model

### Prepare function for sampling
"""

def randomSM( n, probs, sigma=None, deltat = 1):
  '''
  create a sochastic transition matrix with given rates and time step
  '''

  k = np.zeros((n,n))
  for i,j,a in probs:
    if isinstance(a, tuple):
      if len(a)==2:
        k[i,j] = np.random.lognormal(mean=np.log(a[0]), sigma=a[1])
      else:
        k[i,j] = np.random.lognormal(mean=np.log(a[0]), sigma=sigma)
    else:
      k[i,j] = a
  
  rowsums = np.sum(k, axis=1)
  k = k - np.diag(rowsums)
  t = expm(k*deltat)
  if np.any(t<0):
    print('negatives')

  return t


def traj(nsteps, p0, t):
  #traj = [p0.reshape(1,len(p0))]
  traj = [p0]

  for i in range(nsteps):
    q = np.matmul( traj[i] , t)
    traj.append(q)

  return traj

def forward( obs_states, obs_sigma, traj):
  '''
  returns the observed values given probabilities and state_obs
  if state_obs is a single set, uses it for all frames,
  otherwise, state_obs has to be a list of the same lenght as the probs
  '''
  o = []
  if not len(obs_states)==len(traj):
    for pr in traj:
      obs_frame = np.concatenate([np.inner(obs,pr)+np.random.normal(0,obs_sigma,size=1) for obs in obs_states])
      o.append(obs_frame)
  else:
    print((len(obs_states),len(traj)))
    for ifr in range(len(traj)):
      print(ifr)
      pr = traj[ifr]
      obs_st = obs_states[ifr]
      obs_frame = np.concatenate([np.inner(obs,pr)+np.random.normal(0,obs_sigma,size=1) for obs in obs_st])
      o.append(obs_frame)
  return o

def eqprobs(mat):
  n = mat.shape[0]
  w,v = np.linalg.eig(np.transpose(mat))
  iw = np.argmin((w-1)**2)
  #print((w, iw))
  if np.abs(w[iw]-1)>1e-6:
    print(f' perr eigenval is {w[iw]}')
  perr_v = v[:,iw]/np.sum(v[:,iw])
  return perr_v.reshape((1,n))

"""### Define likelihoods"""

def llike_obs_static (piin, obs_list, obsexpS, sigmaexp=.1):
  llike = []
  for sample in obsexpS:
    obscalc = forward(obs_list, 0.0, piin)[0]
    llike.append(np.sum(-(obscalc - sample)**2/(2*sigmaexp**2)))
  return np.sum(np.array(llike))


def llike_pi_static2(piin, sigmafes=2, fes_priors=None):
  '''
  logprob for static probabilities with free-energy priors
  '''
  fesin = np.log(piin) - np.mean( np.log(piin) )
  if fes_priors is not None:
    if isinstance(fes_priors, tuple):
      feccenter, fesstd = fes_priors
      #feccenter = np.log(prior_centers) - np.mean( np.log(prior_centers) )
      fesin = fesin - feccenter
      return np.sum( -(fesin)**2/(2*fesstd**2))
    else:
      fesin = fesin - fes_priors
      return np.sum(-(fesin)**2/(2*sigmafes**2))
  else: 
    return np.sum(-(fesin)**2/(2*sigmafes**2))

def propose_pi(piin,fes_step=1):
  fesin = np.log(piin)
  irand = np.random.choice(range(len(piin[0])),1)
  fes_out = np.copy(fesin)
  fes_out[0][irand] = fes_out[0][irand] + np.sign(np.random.normal())*np.random.normal(fes_step,0.1)

  try:
    ef = np.exp(fes_out)
  except Exception as e:
    print(e)
    print(fes_out)
  return ef/np.sum(ef)


class metropolis(object):
  def __init__(self):
    self.trial_accepted=[]
  
  def accept_proposal(self,likeout,likein):
    if likeout > likein:
      self.trial_accepted.append(1)
      return True
    else:
      if np.exp(likeout-likein) > np.random.uniform():
        self.trial_accepted.append(1)
        return True
      else:
        self.trial_accepted.append(0)
        return False

  def get_fraction(self):
    return np.sum(self.trial_accepted)/len(self.trial_accepted)

##pi 0 with informative priors

print(len(set1))
print(len(set2))

eps = 1e-6
piin = (np.array([[int(i>=len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors = np.log(piin)

for i in range(20):
  piout = propose_pi(piin)
  pbs = ( sum(piout[0][0:len(set1)]), sum(piout[0][len(set1):nstates]))

  lin = llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors)  
  lout = llike_pi_static2(piout, sigmafes=2, fes_priors=fespriors)

  print(( pbs, lout-lin))

### pi inf with informative

eps = 1e-6
piin = (np.array([[int(i<len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors = np.log(piin)
## llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors)  

for i in range(20):
  piout = propose_pi(piin)
  pbs = ( sum(piout[0][0:len(set1)]), sum(piout[0][len(set1):nstates]))

  lin = llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors)  
  lout = llike_pi_static2(piout, sigmafes=2, fes_priors=fespriors)

  print(( pbs, lout-lin))

### pi0 with informative priors
eps = 1e-6
piin = (np.array([[int(i>=len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_t0 = np.log(piin)


### pi inf with informative
eps = 1e-6
piin = (np.array([[int(i<len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_inf = np.log(piin)

## llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors_inf)

"""### Sample pi at equilibrium with ME forward parameters"""

sigmaexp = 2
sigmafes = 2

#beta0 = [0,0,0,0,0,-2,-2,2,-3,-3,-2]
beta0=[0 for i in range(len(cl_to_name))]
theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
obs_list = obs_list_function_base( theta_in, beta0 )

# was July 12, 2022
#piin = (np.zeros((1,nstates))+1)/nstates
#likein = llike_pi_static2(piin,sigmafes=sigmafes) + llike_obs_static(piin, obs_list, exp_inf,sigmaexp=sigmaexp)

eps = 1e-6
piin = (np.array([[int(i<len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_inf = np.log(piin)
likein = llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors_inf) + \
  llike_obs_static(piin, obs_list, exp_inf,sigmaexp=sigmaexp) + \
  llike_TrQI_static(piin, QI_sw, trqi_exp_inf ,sigmaexp=sigmaexp_trqi)

pis_inf = []
acceptance = metropolis()
for i in range(1000):
  piout = propose_pi(piin, fes_step=1)
  likeout = llike_pi_static2(piout,sigmafes=2, fes_priors=fespriors_inf) + \
    llike_obs_static(piout, obs_list, exp_inf,sigmaexp=sigmaexp) + \
    llike_TrQI_static(piout, QI_sw, trqi_exp_inf , sigmaexp=sigmaexp_trqi)

  #acc = accept(likeout,likein)
  #acceptance.append(acc) 
  if acceptance.accept_proposal(likeout,likein):
    piin = piout
    likein = likeout
  pis_inf.append(piin)

probs_inf = pis_inf[400::50]

acceptance.get_fraction()

np.round(probs_inf[0],2)

## check TrQI vals
forward([QI_sw], 0, probs_inf)

#predicted_obsinf

predicted_obsinf = forward(obs_list, 0, probs_inf)
#predicted_obsinf

def compare_vals(val1, val2, labels=None):
  nbars = val1.shape[1]
  pos = np.arange(nbars)
  dodge=.3
  fig, ax = plt.subplots(figsize=(10,4))
  bp1=plt.boxplot(val1,positions=pos, patch_artist=True, boxprops=dict(facecolor='red', color='red'))
  bp2=plt.boxplot(val2,positions=pos+dodge,patch_artist=True,boxprops=dict(facecolor='darkblue', color='darkblue'))
  if labels==None:
    labels=['first','second']
  
  ax.legend([bp1["boxes"][0], bp2["boxes"][0]], labels, loc='upper right')

predicted_obsinf = forward(obs_list, 0, probs_inf)
val1=np.vstack(predicted_obsinf)
val2=np.vstack(exp_inf)
compare_vals(val1,val2, ['predicted','exp'])
pl1=plt.title('Equilibrium observables')

#probs_inf

"""### Sample static model at t=0"""

#theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
#obs_list = obs_list_function( theta_in )
#sigmaexp = 1

#piin = (np.zeros((1,nstates))+1)/nstates
####### set to 0 probs in set1
eps = 1e-6
piin = (np.array([[int(i>=len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_t0 = np.log(piin)


likein = llike_pi_static2(piin,sigmafes=sigmafes,fes_priors=fespriors_t0) + \
  llike_obs_static(piin, obs_list, exp_t0, sigmaexp=sigmaexp) + \
  llike_TrQI_static(piin, QI_sw, trqi_exp_t0 ,sigmaexp=.1)


pis_t0 = []
acceptance = metropolis()
for i in range(1000):
  #print(i)
  piout = propose_pi(piin,fes_step=1)
  likeout = llike_pi_static2(piout,sigmafes=sigmafes,fes_priors=fespriors_t0) + \
    llike_obs_static(piout, obs_list, exp_t0, sigmaexp=sigmaexp) + \
    llike_TrQI_static(piout, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

  #acc = accept(likeout,likein)
  #acceptance.append(acc) 
  if acceptance.accept_proposal(likeout,likein):
    piin = piout
    likein = likeout
  pis_t0.append(piin)

probs_t0 = pis_t0[400::50]
#np.sum(acceptance)/len(acceptance)
acceptance.get_fraction()

np.round(pis_t0[400::50][0],2)

## check TrQI vals
np.round(pis_t0[400::50][0],0)
forward([QI_sw], 0, pis_t0[400::50])

predicted_obst0 = forward(obs_list, 0, pis_t0[400::50])
val1=np.vstack(predicted_obst0)
val2=np.vstack(exp_t0)
compare_vals(val1,val2,['predicted','exp'])
pl1=plt.title('Initial observables')



"""### Sample static model at one lagtime (no TrIQ)"""

#theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
#obs_list = obs_list_function( (b_me, xc_me, xh_me, bc_me, bh_me, b0_me-.8) )
#sigmaexp = 1

exp_t1 = exp_time[0]
piin = (np.zeros((1,nstates))+1)/nstates
likein = llike_pi_static2(piin,sigmafes=sigmafes) + llike_obs_static(piin, obs_list, exp_t1, sigmaexp=sigmaexp)

pis_t1 = []
acceptance = metropolis()
for i in range(1000):
  piout = propose_pi(piin,fes_step=1)
  likeout = llike_pi_static2(piout,sigmafes=sigmafes) + llike_obs_static(piout, obs_list, exp_t1, sigmaexp=sigmaexp)
  if acceptance.accept_proposal(likeout,likein):
    piin = piout
    likein = likeout
  pis_t1.append(piin)

probs_t1 = pis_t1[400::50]
acceptance.get_fraction()

#theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
#obs_list2 = obs_list_function( (b_me, xc_me, xh_me, bc_me, bh_me, b0_me-.8) )

predicted_obst1 = forward(obs_list, 0, probs_t1)
val1=np.vstack(predicted_obst1)
val2=np.vstack(exp_t1)
compare_vals(val1,val2,['pred','exp'])
pl1=plt.title('t=1 observables')

val1=np.vstack(probs_t0)
val2=np.vstack(probs_t1)
compare_vals(val1,val2,['t0','t1'])
pl1=plt.title('probs at t=0 and at t=1')

"""### Setup initial transition matrix

"""

def instantneous_tmat(probs):
  return np.repeat(probs,3,axis=0)

def solveDD(A,pi0_flat, tol=1e-12, niter=20000):
  ndim = A.shape[0]
  def DD(vec):
    return pi0_flat/np.matmul(A, vec)
  vec = (np.zeros(ndim)+1)/ndim
  xi = 0.4
  for i in range(niter):
    dvec = DD(vec)
    vec = (1 - xi) * vec + xi * dvec
    if np.all(np.abs(vec-dvec)<tol):
      break
  #print(DD(vec) - vec)
  return vec

def riequil_tmat(m0_old, piout):
  phiS = np.sqrt(m0_old*np.transpose(m0_old))
  piout_flat = piout.reshape((piout.shape[1],))
  vec = solveDD(phiS, piout_flat, tol = 1e-18, niter = 300)
  matout = np.outer(vec/piout_flat, vec) * phiS
  return matout

from logging import warning
def lagDis2(alpha, Dis, Adj, pi0, tol=1e-18, niter=300, return_values=False):
  pi0_flat = pi0.reshape((pi0.shape[1],))
  #print(pi0_flat)
  if Adj is None:
    n = pi0.shape[1]
    Adj1 = np.zeros((n,n))*0+1
  else:
    Adj1=np.copy(Adj)
    #print(Adj)
  if len(alpha)!=len(Dis):
    warning("o")  
  alphadis = [Dis[i]*alpha[i] for i in range(len(alpha))]
  alphadis = np.sum(np.array(alphadis),axis=0)
  #print(alphadis)
  vec = solveDD(Adj1 * np.exp(alphadis), pi0_flat, tol = tol, niter = niter)
  #print(vec)
  T = np.outer(vec/pi0_flat, vec) * Adj1 * np.exp(alphadis)

  if not return_values:
    return T
  else:
      vals1 = [np.sum(np.matmul(pi0,T*dis)) for dis in Dis]
      return (vals1, T)


def symm(mat):
  return (np.transpose(mat)+mat)*1/2

## extract fes to use for priors
def get_fes_prior(probs_traj):
  probs = np.vstack(probs_traj)
  nstates = probs.shape[1]
  ff = np.log(probs)
  ffmean=np.repeat(np.mean(ff,axis=1)[:,np.newaxis],nstates,axis=1)
  fcent = ff-ffmean
  ffaverage = np.mean(fcent,axis=0)
  ffsted = np.std(fcent,axis=0)
  
  ## calc pis
  piaverage = np.exp(ffaverage)/np.sum(np.exp(ffaverage))
  #piupper = np.exp(ffaverage+ffsted)/np.sum(np.exp(ffaverage+ffsted))
  #pilower = np.exp(ffaverage-ffsted)/np.sum(np.exp(ffaverage-ffsted))
  #pi_std  = (piupperpilower)/2


  f, ax = plt.subplots(1, 2)

  val1 = np.vstack(fcent)
  nbars = val1.shape[1]
  pos = np.arange(nbars)
  dodge=.3
  pt1=ax[0].boxplot(val1,positions=pos, patch_artist=True, boxprops=dict(facecolor='red', color='red'))
  pl1=ax[0].errorbar(pos+dodge, ffaverage, yerr = ffsted, c='green', marker='x')
  
  # sample pis
  ff_sample=np.random.multivariate_normal(ffaverage, np.diag(ffsted**2),10)
  pi_sample = np.exp(ff_sample)
  pi_norm = np.sum(pi_sample,axis=1)
  pi_sample_norm = pi_sample/np.repeat(pi_norm[:,np.newaxis],nstates,axis=1)

  val1 = probs
  val2 = pi_sample_norm
  nbars = val1.shape[1]
  pos = np.arange(nbars)
  dodge=.3
  pt1=ax[1].boxplot(val1,positions=pos, patch_artist=True, boxprops=dict(facecolor='red', color='red'))
  pt1=ax[1].boxplot(val2,positions=pos+dodge,patch_artist=True,boxprops=dict(facecolor='darkblue', color='darkblue'))
  pt1=ax[1].plot(pos-dodge,   piaverage, marker='o', color='pink')


  return (ffaverage, ffsted, piaverage.reshape((1,nstates)))

ff_mean_inf, ff_std_inf, pi_mean_inf = get_fes_prior(probs_inf)
ff_mean_t0, ff_std_t0, pi_mean_t0 = get_fes_prior(probs_t0)
ff_mean_t1, ff_std_t1, pi_mean_t1 = get_fes_prior(probs_t1)

(pi_mean_inf>1e-2) | (pi_mean_t0>1e-2) | (pi_mean_t1>1e-2)

pi_mean_t0

#regul_prob(pi_mean_t0,.01)

from scipy.optimize import minimize
def regul_prob(vec, prior):
  return (vec+prior)/np.sum(vec+prior)

## use the probs from averages
if False:
  prob0_sample = regul_prob(probs_t0[0],.5) # pi_mean_t0
  probinf_sample = regul_prob(probs_inf[0],.5) #pi_mean_inf
else:
  prob0_sample = regul_prob(pi_mean_t0,.005)
  probinf_sample = regul_prob(pi_mean_inf,.005)

## vals = forward(obs_list, 0, [probinf_sample])[0]
## use average observables from the first lag
## vals = predicted_obst1[1] #exp_time[0][1]

subset=slice(0,1)
vals= np.mean(np.vstack(exp_time[0]),axis=0)[subset]
Diss = [symm(np.outer(obsv/probinf_sample,prob0_sample)) for obsv in  obs_list[subset]]

def optf(alphas): 
  calc_vals, mat = lagDis2(alphas, Diss, None, probinf_sample,return_values=True)
  return np.sum((np.array(calc_vals) - vals)**2)

alphas = [0.001 for v in vals]
res = minimize(optf, alphas, tol=1e-15)
calc_vals, mat = lagDis2(res.x, Diss, None, probinf_sample,return_values=True,niter=1000,tol=1e-18)

if res.success:
  print(res)
  print(f'difference in vals: {vals-calc_vals}')
  print(f'difference in probs: {eqprobs(mat)-probinf_sample}')
else:
  print('WARNING, WARNING: optimization failed')
  print(res)

  print(f'difference in vals: {vals-calc_vals}')
  print(f'difference in probs: {eqprobs(mat)-probinf_sample}')

np.round(mat,3)

prob10trjs = [pi_mean_t0]
pi_mean_t0.sum()

times = np.arange(0,60)

for i in times[1:]:
  prt = np.matmul(prob10trjs[i-1], mat)
  prob10trjs.append(prt)
  if np.abs(prt.sum()-1)>1e-5:
    print(prt.sum())

exp_time[0]

from IPython.core.pylabtools import figsize
#exp_time[i] is i-th intermediate time,


pred_time0 = forward(obs_list,0, prob10trjs)
pred_time = np.vstack(pred_time0)
# each col is one observation out of the nobs
# nobs = len(keep_peptides)

f, ax = plt.subplots(nobs,1,figsize=(10,6))
for iobs,obs in  enumerate(keep_peptides):
  ax[iobs].plot(times, pred_time[:,iobs], label=obs)
  exppts0 = [dp[iobs] for dp in exp_t0]
  exppts_inf = [dp[iobs] for dp in exp_inf]
  
  ax[iobs].plot(np.zeros((len(exppts0),)), exppts0,marker='o')
  ax[iobs].plot(0, np.mean(exppts0),marker='x')
  ax[iobs].plot(np.zeros((len(exppts0),))+np.max(times), exppts_inf ,marker='o')
  ax[iobs].plot(np.max(times), np.mean(exppts_inf) ,marker='x')

  for ix, timex in enumerate(time_powers):
    exppt_time = [dp[iobs] for dp in exp_time[ix]]
    ax[iobs].plot(np.zeros((len(exppt_time),))+timex, exppt_time,marker='o')
    ax[iobs].plot(timex, np.mean(exppt_time),marker='x')

  ax[iobs].set_ylim([0, 40])
  ax[iobs].legend()

plt.plot(forward(obs_list,0, prob10trjs),c='red')

#cols = ['red','orange','blue']
#for ot in otrjs:
#  ost = np.vstack(ot[::10]).transpose()
#  for i,a in enumerate(ost):
#    plt.plot(a, c=plt.get_cmap('Paired')(i),label=f'trj {i}')

"""## Sample thermodynamic models with Bayesian forward **model**

### Sample pi at equilibirum and forward model parameters
"""

import warnings
#warnings.filterwarnings("error")
warnings.filterwarnings("default")

#prob0_sample = regul_prob(pi_mean_t0,.01)
#probinf_sample = regul_prob(pi_mean_inf,.01)
#print(probinf_sample)

eps = 1e-6
piin = (np.array([[int(i<len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_inf = np.log(piin)


theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
beta0_in = beta0
#obs_list_in = obs_list_function( theta_in )
obs_list_in = obs_list_function_base( theta_in ,beta0_in)

pis_inf = []
theta_s2 = []
beta0_s2 = []
obs_s2 = []

acceptance_inf = metropolis()
acceptance_theta = metropolis()

nepochs=100

for i in range(nepochs):

  for j in range(nstates):
    likein = llike_pi_static2(piin,sigmafes=2, fes_priors=fespriors_inf) + \
      llike_obs_static(piin, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piin, QI_sw, trqi_exp_inf ,sigmaexp=sigmaexp_trqi)

    piout = propose_pi(piin, fes_step=.5)
    likeout = llike_pi_static2(piout,sigmafes=2, fes_priors=fespriors_inf) + \
      llike_obs_static(piout, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piout, QI_sw, trqi_exp_inf ,sigmaexp=sigmaexp_trqi)
      
    if acceptance_inf.accept_proposal(likeout,likein):
        piin = piout
  pis_inf.append(piin)

  for j in range(10):
    try:
      likein = llike_obs_static(piin, obs_list_in, exp_inf, sigmaexp=sigmaexp)
      theta_out, llog_theta = vpriors.propose_theta(theta_in)  
      beta0_out, llog_beta0 = vpriors.propose_beta0(beta0_in)   
      obs_list_out = obs_list_function_base( theta_out , beta0_out)
      likeout = llike_obs_static(piin, obs_list_out, exp_inf, sigmaexp=sigmaexp) + llog_theta + llog_beta0
    except Exception as e:
      print(f'{e}')

    if acceptance_theta.accept_proposal(likeout,likein):
      theta_in = theta_out
      obs_list_in = obs_list_out
      beta0_in = beta0_out
  theta_s2.append(theta_in)
  beta0_s2.append(beta0_in)
  obs_s2.append(obs_list_in)

probs_inf = pis_inf[40::5]
theta_s2str = theta_s2[40::5]
beta0_s2str = beta0_s2[40::5]
obs_s2str = obs_s2[40::5]

## check TrQI vals
np.round(probs_inf,0)
forward([QI_sw], 0, probs_inf)

theta_s2str

print((acceptance_theta.get_fraction(),acceptance_inf.get_fraction()))

#probs_t0

predicted_obs_inf = forward(obs_s2str, 0, probs_inf)
val1=np.vstack(predicted_obs_inf)
val2=np.vstack(exp_inf)
compare_vals(val1,val2,['pred','exp'])
pl1=plt.title('t=inf observables')

ff_mean_inf, ff_std_inf, pi_mean_inf = get_fes_prior(probs_inf)

"""### Sample t=0 with $\theta$'s"""

#piin = (np.zeros((1,nstates))+1)/nstates

eps = 1e-6
piin = (np.array([[int(i>=len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_t0 = np.log(piin)


theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
beta0_in = beta0
obs_list_in = obs_list_function_base( theta_in ,beta0_in)

pis_t0 = []
theta_s1 = []
beta0_s1 = []
obs_s1 = []

acceptance_t0 = metropolis()
acceptance_theta = metropolis()

nepochs=100

for i in range(nepochs):

  for j in range(nstates):
    likein = llike_pi_static2(piin,sigmafes=sigmafes, fes_priors=fespriors_t0) + \
      llike_obs_static(piin, obs_list_in, exp_t0, sigmaexp=sigmaexp) +\
      llike_TrQI_static(piin, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

    piout = propose_pi(piin, fes_step=.5)
    likeout = llike_pi_static2(piout,sigmafes=sigmafes, fes_priors=fespriors_t0) + \
      llike_obs_static(piout, obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piout, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

    if acceptance_t0.accept_proposal(likeout,likein):
        piin = piout
  pis_t0.append(piin)

  for j in range(10):
    try:
      likein = llike_obs_static(piin, obs_list_in, exp_t0, sigmaexp=sigmaexp)
      theta_out, llog_theta = vpriors.propose_theta(theta_in) 
      beta0_out, llog_beta0 = vpriors.propose_beta0(beta0_in)   
      obs_list_out = obs_list_function_base( theta_out, beta0_out )
      likeout = llike_obs_static(piin, obs_list_out, exp_t0, sigmaexp=sigmaexp) + llog_theta + llog_beta0
    except Exception as e:
      print(f'{e}')

    if acceptance_theta.accept_proposal(likeout,likein):
      theta_in = theta_out
      beta0_in = beta0_out
      obs_list_in = obs_list_out

  theta_s1.append(theta_in)
  beta0_s1.append(beta0_in)
  obs_s1.append(obs_list_in)

probs_t0 = pis_t0[40::5]
theta_s1str = theta_s1[40::5]
beta0_s1str = beta0_s1[40::5]

obs_s1str = obs_s1[40::5]

acceptance_theta.get_fraction()

theta_s1str

## check TrIQ vals
np.round(probs_t0,0)
forward([trQI_sw], 0, probs_t0)

for pii in probs_t0:
  aa = (sum(pii[0][0:len(set1)]), 
        sum(pii[0]))
  print(np.round(aa,2))

print((acceptance_theta.get_fraction(),acceptance_t0.get_fraction()))

predicted_obs_t0 = forward(obs_s1str, 0, probs_t0)
val1=np.vstack(predicted_obs_t0)
val2=np.vstack(exp_t0)
compare_vals(val1,val2,['pred','exp'])
pl1=plt.title('t=0 observables')



"""### Sample $\theta$'s with $\pi_0$ and $\pi_{\rm inf}$

We sample the static likelihoods and the thetas.

$$p(\pi_{inf}, \pi_0, \theta \vert D) = 
 p(D \vert \pi_{inf}, \pi_0, \theta ) p(\pi_{inf}) p(\pi_0) p(\theta) $$
 $$ 
 = p(D_0 \vert \pi_0 \theta ) 
 p(D_{inf} \vert \pi_{inf} \theta ) 
 p(\pi_{inf}) p(\pi_0) p(\theta) $$

 so:

* sample pi0 from   p(D0 \vert pi_0 theta ) p(pi_0) 
* sample piinf from  p(Dinf \vert pi_inf theta ) p(pi_inf) 
* sample theta from  p(D0 \vert pi_0 theta ) 
   p(Dinf \vert pi_inf theta ) p(theta)
"""

sigmaexp

#piin_0 = (np.zeros((1,nstates))+1)/nstates

eps = 1e-6
piin_0 = (np.array([[int(i>=len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_t0 = np.log(piin_0)

## was July 7, 2022
## piin_inf = (np.zeros((1,nstates))+1)/nstates
### pi inf with informative
eps = 1e-6
piin_inf = (np.array([[int(i<len(set1)) for i in range(nstates)]]) + eps)/(nstates*eps + len(set2))
fespriors_inf = np.log(piin_inf)
## llike_pi_static2(piin, sigmafes=2, fes_priors=fespriors_inf)  

# theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
# obs_list_in = obs_list_function_base( theta_in, beta0_in )


theta_in = (b_me, xc_me, xh_me, bc_me, bh_me, b0_me)
beta0_in = beta0
#obs_list_in = obs_list_function( theta_in )
obs_list_in = obs_list_function_base( theta_in, beta0_in)

pis_t0 = []
pis_inf = []
theta_s3 = []
beta0_s3=[]
obs_s3 = []

acceptance_t0 = metropolis()
acceptance_inf = metropolis()
acceptance_theta = metropolis()

nepochs=100

for i in range(nepochs):

  for j in range(nstates):
    likein = llike_pi_static2(piin_0,sigmafes=sigmafes, fes_priors=fespriors_t0) + \
      llike_obs_static(piin_0, obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piin_0, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

    piout_0 = propose_pi(piin_0, fes_step=1.)
    likeout = llike_pi_static2(piout_0,sigmafes=sigmafes,fes_priors=fespriors_t0) + \
      llike_obs_static(piout_0, obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piout_0, QI_sw, trqi_exp_t0 ,sigmaexp=.1)
    if acceptance_t0.accept_proposal(likeout,likein):
        piin_0 = piout_0
  pis_t0.append(piin_0)

  for j in range(nstates):
    likein = llike_pi_static2(piin_inf,sigmafes=sigmafes, fes_priors=fespriors_inf) + \
      llike_obs_static(piin_inf, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piin_inf, QI_sw, trqi_exp_inf ,sigmaexp=0.02)

    piout_inf = propose_pi(piin_inf, fes_step=1.)
    likeout = llike_pi_static2(piout_inf,sigmafes=sigmafes, fes_priors=fespriors_inf) + \
      llike_obs_static(piout_inf, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
      llike_TrQI_static(piout_inf, QI_sw, trqi_exp_inf ,sigmaexp=0.02)

    if acceptance_inf.accept_proposal(likeout,likein):
        piin_inf = piout_inf
  pis_inf.append(piin_inf)


  for j in range(10):
    try:
      likein = llike_obs_static(piin_0,   obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
               llike_obs_static(piin_inf, obs_list_in, exp_inf, sigmaexp=sigmaexp)
      theta_out, llog_theta = vpriors.propose_theta(theta_in)
      beta0_out, llog_beta0 = vpriors.propose_beta0(beta0_in)   
      obs_list_out = obs_list_function_base( theta_out, beta0_out )
      likeout = llike_obs_static(piin_0, obs_list_out,   exp_t0,  sigmaexp=sigmaexp) + \
                llike_obs_static(piin_inf, obs_list_out, exp_inf, sigmaexp=sigmaexp) + \
                llog_theta + llog_beta0
    except Exception as e:
      print(f'{e}')

    if acceptance_theta.accept_proposal(likeout,likein):
      theta_in = theta_out
      beta0_in = beta0_out
      obs_list_in = obs_list_out    
  theta_s3.append(theta_in)
  beta0_s3.append(beta0_in)
  obs_s3.append(obs_list_in)

probs_t0 = pis_t0[40::5]
probs_inf = pis_inf[40::5]

theta_s3str = theta_s3[40::5]
beta0_s3str = beta0_s3[40::5]
obs_s3str = obs_s3[40::5]

print((acceptance_theta.get_fraction(),acceptance_inf.get_fraction(),acceptance_t0.get_fraction()))

## check TrQI vals at inf
print(np.round(probs_inf[0] ,2))
forward([QI_sw], 0, probs_inf)

## check TrQI vals at t0
print(np.round(probs_t0 [0],2))
forward([QI_sw], 0, probs_t0)

predicted_obs_inf = forward(obs_s3str, 0, probs_inf)
val1=np.vstack(predicted_obs_inf)
val2=np.vstack(exp_inf)
compare_vals(val1,val2,['pred','exp'])
pl1=plt.title('t=inf observables')

predicted_obs_t0 = forward(obs_s3str, 0, probs_t0)
val1=np.vstack(predicted_obs_t0)
val2=np.vstack(exp_t0)
compare_vals(val1,val2,['pred','exp'])
pl1=plt.title('t=0 observables')

ff_mean_inf, ff_std_inf, pi_mean_inf = get_fes_prior(probs_inf)
ff_mean_t0, ff_std_t0, pi_mean_t0 = get_fes_prior(probs_t0)

"""### Guess tmat from observations at t=1"""

pi_mean_inf

from scipy.optimize import minimize
def regul_prob(vec, prior):
  return (vec+prior)/np.sum(vec+prior)

## use the probs from averages
if False:
  prob0_sample = regul_prob(probs_t0[0],.0005) # pi_mean_t0
  probinf_sample = regul_prob(probs_inf[0],.0005) #pi_mean_inf
else:
  prob0_sample = regul_prob(pi_mean_t0,.01)
  probinf_sample = regul_prob(pi_mean_inf,.01)

## vals = forward(obs_list, 0, [probinf_sample])[0]
## use average observables from the first lag
## vals = predicted_obst1[1] #exp_time[0][1]

subset=slice(0,1)
vals= np.mean(np.vstack(exp_time[0]),axis=0)[subset]
Diss = [symm(np.outer(obsv/probinf_sample, prob0_sample)) for obsv in  obs_list[subset]]

def optf(alphas): 
  calc_vals, mat = lagDis2(alphas, Diss, None, probinf_sample,return_values=True)
  return np.sum((np.array(calc_vals) - vals)**2)

alphas = [0.01 for v in vals]
res = minimize(optf, alphas, tol=1e-18)
calc_vals, mat = lagDis2(res.x, Diss, None, probinf_sample,return_values=True, niter=3000, tol=1e-19)

prob0_sample

print(np.sum(pi_mean_t0))
print(np.sum(probinf_sample))

np.sum(mat,axis=1)

if res.success:
  print(res)
  print(f'difference in vals: {vals-calc_vals}')
  print(f'difference in probs: {eqprobs(mat)-probinf_sample}')
else:
  print('WARNING, WARNING: optimization failed')
  print(res)

  print(f'difference in vals: {vals-calc_vals}')
  print(f'difference in probs: {eqprobs(mat)-probinf_sample}')

print(eqprobs(mat)[0][set1])
print(eqprobs(mat)[0][set2])

prob10trjs = [pi_mean_t0]
print(pi_mean_t0.sum())

times = np.arange(0,60)

for i in times[1:]:
  prt = np.matmul(prob10trjs[i-1], mat)
  prob10trjs.append(prt)
  if np.abs(prt.sum()-1)>1e-15:
    print(prt.sum())
  else:
    print(f'OK {prt.sum()}')

pred_time0 = forward(obs_s3str[5],0, prob10trjs)
pred_time = np.vstack(pred_time0)
# each col is one observation out of the nobs
# nobs = len(keep_peptides)

f, ax = plt.subplots(nobs,1,figsize=(10,6))
for iobs,obs in  enumerate(keep_peptides):
  ax[iobs].plot(times, pred_time[:,iobs], label=obs)
  exppts0 = [dp[iobs] for dp in exp_t0]
  exppts_inf = [dp[iobs] for dp in exp_inf]
  
  ax[iobs].plot(np.zeros((len(exppts0),)), exppts0,marker='o')
  ax[iobs].plot(0, np.mean(exppts0),marker='x')
  ax[iobs].plot(np.zeros((len(exppts0),))+np.max(times), exppts_inf ,marker='o')
  ax[iobs].plot(np.max(times), np.mean(exppts_inf) ,marker='x')

  for ix, timex in enumerate(time_powers):
    exppt_time = [dp[iobs] for dp in exp_time[ix]]
    ax[iobs].plot(np.zeros((len(exppt_time),))+timex, exppt_time,marker='o')
    ax[iobs].plot(timex, np.mean(exppt_time),marker='x')

  ax[iobs].set_ylim([0, 40])
  ax[iobs].legend()

print(np.round(mat,3))
plt.plot(forward(obs_list,0, prob10trjs),c='red')

theta_samp=theta_s3str[0]
beta0_samp = beta0_s3str[0]
obs_byobs2 = obs_list_function_base(theta_samp,beta0_samp)

ob_clusters2 = pd.DataFrame(
    np.vstack(obs_byobs2), 
    columns=[f'cl{i}' for i in range(len(dists_all_clusters))],
    index = keep_peptides)

ob_clusters2['mean_t0'] = mexp_t0
ob_clusters2['mean_inf'] = mexp_inf

ob_clusters2



"""### Read back saved data"""

### read all
if True:
  #saved_data = [f'save_00{i}.npy' for i in [1,2,3,4,5,6,7]]
  saved_data = [f'{basedir}/sampling_f05/save_{10:03}.npy']
  for filename in saved_data:
    with open(filename, 'rb') as f:
      r_probs_t0=np.load(f)
      rprobs_inf=np.load(f)
      rtheta_s4str=np.load(f)
      rbeta0_s4str=np.load(f)
      robs_s4str=np.load(f)
      rmat_s4str=np.load(f)
      print(f'{filename} -> {r_probs_t0.shape}')

probs_t0 = r_probs_t0
probs_inf=rprobs_inf  
theta_s4str=rtheta_s4str  
beta0_s4str=rbeta0_s4str  
obs_s4str  =robs_s4str  
mat_s4str= rmat_s4str

"""

### Initial guess of tmat before sampling the full model"""

# rescale tmat
# mat0=np.copy(mat)

mat0 = rmat_s4str[10]
print(sum(eqprobs(mat0)[0][set1]))
print(sum(eqprobs(mat0)[0][set2]))

mat = np.linalg.matrix_power(mat0,2)

