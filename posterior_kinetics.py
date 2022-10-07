# -*- coding: utf-8 -*-
"""Zenodo HDX+TrIQ sampling

# Prepare environment
"""

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


"""## Sample kinetic model with forward model parameters

We sample the full model:

$$p(K, \pi_{inf}, \pi_0, \theta \vert D) = $$

$$
 = p(D \vert K, \pi_{inf}, \pi_0, \theta ) 
 p(K\vert \pi_{inf}) p(\pi_{inf}) p(\pi_0) p(\theta) $$
 $$ 
= p(D_0 \vert \pi_0 \theta ) 
 p(D_{inf} \vert \pi_{inf} \theta ) 
 p(D_t \vert K \pi_{inf} \pi_0 \theta ) 
 p(\pi_{inf}) p(\pi_0) p(\theta) p(K \vert \pi_{inf} )$$

 so:



* sample pi0 from $p(D_0 \vert \pi_0 \theta ) 
p(D_t \vert K \pi_{inf} \pi_0 \theta )  p(\pi_0)$
* sample piinf,K from  $ p(D_{inf} \vert \pi_{inf} \theta ) 
 p(D_t \vert K \pi_{inf} \pi_0 \theta ) 
 p(\pi_{inf}) p(K \vert \pi_{inf} )$:
 first sample $\pi_{inf}$ updating $K$ to the max-cal matrix with the new equilibrium distribution; then - with fixed $\pi_{inf}$ sample $K$
* sample theta from  $p(D_0 \vert \pi_0 \theta ) 
 p(D_{inf} \vert \pi_{inf} \theta ) 
 p(D_t \vert K \pi_{inf} \pi_0 \theta ) 
 p(\theta)$
"""

def propose_Tmat(mat_in0, pout=None, Adj=None, verbose=False):
  '''
  propose a new tmat starting from an incoming tmat,
  changes of eq. probability are allowed, in case p_infout is not none, 

  '''
  nstates = mat_in0.shape[0]
  if not pout is None:
    mat_in = riequil_tmat(mat_in0, pout)
    probs = pout.reshape((nstates,))
  else:
    mat_in = mat_in0
    probs = eqprobs(mat_in0).reshape((nstates,))
  if Adj is None:
    Adj = mat_in0*0+1
  Adj = Adj - np.diag(np.diag(Adj))
  transitions=np.where(Adj==1)
  itrans = np.random.choice(range(len(transitions[0])))
  i = transitions[0][itrans]
  j = transitions[1][itrans]

  delta_min = max(-mat_in[i,i],-mat_in[j,j]*probs[j]/probs[i]) 
  delta_max=mat_in[i,j]
  delta = np.random.uniform(low=delta_min, high=delta_max, size=1)

  if verbose:
    print(f'(i,j) = {(i,j)}, delta= {delta} in ({delta_min},{delta_max})')

  mat_out = np.copy(mat_in)
  mat_out[i,i] = mat_in[i,i]+delta
  mat_out[j,j] = mat_in[j,j]+delta* probs[i]/probs[j]
  mat_out[i,j] = mat_out[i,j]-delta
  mat_out[j,i] = mat_out[j,i]-delta* probs[i]/probs[j]

  aux_in  = np.array([mat_in[i,i], mat_in[j,j], mat_in[i,j], mat_in[j,i]])
  aux_out = np.array([mat_out[i,i], mat_out[j,j], mat_out[i,j], mat_out[j,i]])
  aux_probs  = np.array([probs[i], probs[j], probs[i], probs[j]])
    
  ll_in = np.sum( np.log(aux_in )*(-1/2 + aux_probs*aux_in) ) 
  ll_out = np.sum( np.log(aux_out)*(-1/2 + aux_probs*aux_out) ) 
  
  return mat_out , ll_in - ll_out 


def llike_obs_dyn_point (pi0in,  obms_list, obsexpS, T0, lags, sigmaexp=.1):
  llike = []

  for i,lag in enumerate(lags):
    pi_t = np.matmul(pi0in, matrix_power(T0, lags[i]))
    obscalc = forward(obs_list, 0.0, pi_t)[0]
    llike1=[]
    for sample in obsexpS:
      llike1.append(np.sum(-(obscalc - sample)**2/(2*sigmaexp**2)))
    #print(llike1)
  llike.append(np.sum(np.array(llike1)))
  return np.sum(np.array(llike))

len(beta0_in)

#init probs
print(probs_t0[3])

partial_sums = (
    sum(probs_t0[3][0][0:len(set1)]),
    sum(probs_t0[3][0])
)

print(partial_sums)

#inf probs
print(probs_inf[3])

partial_sums = (
    sum(probs_inf[3][0][0:len(set1)]),
    sum(probs_inf[3][0])
)

print(partial_sums)

mat

np.round(eqprobs(mat0)[0][set2],2)

sum(eqprobs(mat0)[0][set2])

lagsp = time_powers

acceptance_pi0 = metropolis()
acceptance_pinf = metropolis()
acceptance_T = metropolis()
acceptance_theta = metropolis()

# to collect the sampling results
pi0_sample=[]
piinf_sample=[]
T_sample=[]
theta_sample = []
beta0_sample= []
obs_sample = []

# priors for probabilities
prior_centers_t0 = (ff_mean_t0,ff_std_t0)
prior_centers_inf = (ff_mean_inf,ff_std_inf)

# initial values
T_in = mat
pi0_in = probs_t0[3]
piinf_in = eqprobs(T_in)
theta_in = theta_s3[0]   
beta0_in = beta0_s3[0]
obs_list_in = obs_list_function_base( theta_in, beta0_in )

fesstep = .8
nepochs = 40 + 5*200
for isample in range(nepochs):

  # sample pi0 ------------------------------------ 
  llike_in = llike_pi_static2(pi0_in, sigmafes=2, fes_priors=fespriors_t0) + \
     llike_obs_static(   pi0_in, obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
     llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_in, lagsp, sigmaexp=sigmaexp) + \
     llike_TrQI_static(pi0_in, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

  for ipi0 in range(nstates):
    pi0_prop = propose_pi(pi0_in,fes_step=fesstep)
    llike_out = llike_pi_static2(pi0_prop, sigmafes=2, fes_priors=fespriors_t0) + \
      llike_obs_static(   pi0_prop, obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
      llike_obs_dyn_point(pi0_prop, obs_list_in, exp_time, T_in, lagsp, sigmaexp=sigmaexp) + \
      llike_TrQI_static(pi0_prop, QI_sw, trqi_exp_t0 ,sigmaexp=.1)

    if acceptance_pi0.accept_proposal(llike_out,llike_in):
       pi0_in = pi0_prop
       llike_in=llike_out
  # sample pi0 ------------------------------------ 

  # sample piinf ------------------------------------
  llike_in = llike_pi_static2(piinf_in, sigmafes=2, fes_priors=fespriors_inf) + \
     llike_obs_static(piinf_in, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
     llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_in, lagsp, sigmaexp=sigmaexp) + \
     llike_TrQI_static(piinf_in, QI_sw, trqi_exp_inf ,sigmaexp=.005) 
     #sigmaexp was 0.01 for sampling_f05/ chain 10

  for ipiinf in range(nstates):
    piinf_prop = propose_pi(piinf_in,fes_step=fesstep)
    T_prop, Tllike = propose_Tmat(T_in,pout=piinf_prop)
    llike_out = llike_pi_static2(piinf_prop, sigmafes=2, fes_priors=fespriors_inf) + \
      llike_obs_static(piinf_prop, obs_list_in, exp_inf,sigmaexp=sigmaexp) + \
      llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_prop, lagsp, sigmaexp=sigmaexp) + \
      Tllike + \
      llike_TrQI_static(piinf_prop, QI_sw, trqi_exp_inf ,sigmaexp=.005)

    if acceptance_pinf.accept_proposal(llike_out,llike_in):
      piinf_in = piinf_prop
      T_in = T_prop
      llike_in=llike_out
  # sample piinf ------------------------------------

  # sample T ------------------------------------
  llike_in = llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_in, lagsp, sigmaexp=sigmaexp)
  for iT in range(nstates*nstates):
     T_prop, Tllike = propose_Tmat(T_in)
     llike_out = llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_prop, lagsp, sigmaexp=sigmaexp) + Tllike
     if acceptance_T.accept_proposal(llike_out,llike_in):
        T_in = T_prop
        llike_in=llike_out
  # sample T ------------------------------------
  


  # sample theta ------------------------------------
  for j in range(10):
    try:
      likein = llike_obs_static(pi0_in,   obs_list_in, exp_t0, sigmaexp=sigmaexp) + \
               llike_obs_static(piinf_in, obs_list_in, exp_inf, sigmaexp=sigmaexp) + \
               llike_obs_dyn_point(pi0_in, obs_list_in, exp_time, T_in, lagsp, sigmaexp=sigmaexp)
      theta_out, llog_theta = vpriors.propose_theta(theta_in)
      beta0_out, llog_beta0 =      vpriors.propose_beta0(beta0_in)
      obs_list_out = obs_list_function_base( theta_out,beta0_out )
      likeout = llike_obs_static(   pi0_in, obs_list_out, exp_t0,  sigmaexp=sigmaexp) + \
                llike_obs_static( piinf_in, obs_list_out, exp_inf,   sigmaexp=sigmaexp) + \
                llike_obs_dyn_point(pi0_in, obs_list_out, exp_time, T_in, lagsp, sigmaexp=sigmaexp) + \
                llog_theta + llog_beta0
    except Exception as e:
      print(f'{e}')
    if acceptance_theta.accept_proposal(likeout,likein):
      theta_in = theta_out
      beta0_in = beta0_out
      obs_list_in = obs_list_out
  # sample theta ------------------------------------

  if isample % 10 == 0:
    probchecks = (np.sum(eqprobs(T_in)[0][set2]), np.sum(piinf_in[0][set2]), np.sum(pi0_in[0][set1]))
    print(probchecks)

    print(f'{isample} acceptance' +
       f'T: {acceptance_T.get_fraction():.3f}, '+
       f'pi0: {acceptance_pi0.get_fraction():.3f}, '+
       f'piinf: {acceptance_pinf.get_fraction():.3f}, ' +
       f'theta: {acceptance_theta.get_fraction():.3f}')

  pi0_sample.append(pi0_in)
  piinf_sample.append(piinf_in)
  T_sample.append(T_in)
  theta_sample.append(theta_in)
  beta0_sample.append(beta0_in)
  obs_sample.append(obs_list_in)

## stride and burnin
probs_t0 = pi0_sample[40::5]
probs_inf = piinf_sample[40::5]
theta_s4str = theta_sample[40::5]
beta0_s4str = beta0_sample[40::5]
obs_s4str = obs_sample[40::5]
mat_s4str = T_sample[40::5]

len(mat_s4str)

## check TrQI vals at inf
print(np.round(probs_inf[0] ,2))
forward([QI_sw], 0, probs_inf[0::10])

## check TrQI vals at t0
print(np.round(probs_t0 [0],2))
forward([QI_sw], 0, probs_t0[0::10])

#save_chain_number = 10
#f'save_{save_chain_number:03}.npy'

np.array(probs_t0).shape

"""### Save sampled basic data"""

save_chain_number = 3

with open(f'{basedir}/sampling_g02/save_{save_chain_number:03}.npy', 'wb') as f:
  for data in [probs_t0,probs_inf,theta_s4str,beta0_s4str,obs_s4str,mat_s4str]:
    data2 = np.array(data)
    print(data2.shape)
    np.save(f, data2)

### read all
#if False:
#  saved_data = [f'save_00{i}.npy' for i in [1,2,3,4,5,6,7]]
#  saved_data = [f'{basedir}/sampling_f05/save_{10:03}.npy']
#  for filename in saved_data:
#    with open(filename, 'rb') as f:
#      r_probs_t0=np.load(f)
#      rprobs_inf=np.load(f)
#      rtheta_s4str=np.load(f)
#      rbeta0_s4str=np.load(f)
#      robs_s4str=np.load(f)
#      rmat_s4str=np.load(f)
#      print(f'{filename} -> {r_probs_t0.shape}')



import yaml

ydump = {}
ydumppep = {}
for i,kp in enumerate(keep_peptides):
  ydumppep[f'peptide{i}']=kp
ydump['peptides'] = ydumppep

ydumpclusters={}
for i,trn in enumerate(tr_names):
  if i in keep_states:
    itm = {'index':i}
    if isinstance(trn,tuple):
      itm['filename'] = trn[0]
      itm['frames'] = trn[1]
    else:  
      itm['filename'] = trn
    ydumpclusters[f'cluster{i}'] = itm

ioffset=i
for i,trn in enumerate(tr_names2):
  if i in keep_states2:
    itm = {'index':i+ioffset}
    if isinstance(trn,tuple):
      itm['filename'] = trn[0]
      itm['frames'] = trn[1]
    else:  
      itm['filename'] = trn
    ydumpclusters[f'cluster{i+ioffset}'] = itm


ydump['clusters']=ydumpclusters
ydump['fragrecords'] = df_frags.reset_index().to_dict(orient='records')

ydobs = {}
ydob = {'label':lab_t0, 'power':itimes[0]}
for i,o in enumerate(exp_t0):
  ydob[f'observation{i}'] = ';'.join([f'{v}' for v in o]) 
ydobs['timepoint_0'] = ydob

ydob = {'label':lab_inf, 'power':itimes[-1]}
for i,o in enumerate(exp_inf):
  ydob[f'observation{i}'] = ';'.join([f'{v}' for v in o]) 
ydobs['timepoint_inf'] = ydob

for it in range(1, len(itimes)-1):
  ydob = {'label':timelags[it], 'power':itimes[it]}
  for i,o in enumerate(exp_time[it-1]):
    ydob[f'observation{i}'] = ';'.join([f'{v}' for v in o]) 
  ydobs[f'timepoint_{it}'] = ydob


ydump['experimental_data'] = ydobs

with open(f'{basedir}sampling_g02/data.yaml', 'w') as f:
  yaml.dump(ydump, f, explicit_start=True, default_flow_style=False)





"""### Other analysis"""

# use this to read basic data to perform the full analsysis
#saved_data = [f'{basedir}/sampling_011/save_00{i}.npy' for i in [4,5,6,7]]
#
#for filename in [saved_data[3]]:
#  with open(filename, 'rb') as f:
#    probs_t0=np.load(f)
#    probs_inf=np.load(f)
#    theta_s4str=np.load(f)
#    beta0_s4str=np.load(f)
#    obs_s4str=np.load(f)
#    mat_s4str=np.load(f)
#    print(f'{filename} -> {probs_t0.shape}')

#sample evolution

# check that all the mats are ok...
for tmat in mat_s4str:
  px = eqprobs(tmat)
  if np.abs(np.sum(px)-1)> 1e-8:
    print(np.sum(px))

times = np.arange(0,120)

all_evols_p = []
all_evols_o = []

for i in range(len(mat_s4str)):
  prob10trjs = [probs_t0[i]]
  mat = mat_s4str[i]
  obs = obs_s4str[i]
  
  for it in times[1:]:
    prt = np.matmul(prob10trjs[it-1], mat)
    prob10trjs.append(prt)
    if np.abs(prt.sum()-1)>1e-5:
      print(prt.sum())
      
  pred_time0 = forward(obs,0, prob10trjs)
  pred_time = np.vstack(pred_time0)

  all_evols_p.append(prob10trjs)
  all_evols_o.append(pred_time)

beta0_means = [np.mean(bt0) for bt0 in beta0_s4str]
plt.plot(beta0_means)

"""#### Plot obs evolution"""

#int(nobs/3)
#print(nobs)
#chx = np.array([-7,-6,-1,-2,0,-7,10,8,7,5])
#chx = np.array([0,0,0,-5,0,-5,-5,0,0,0,0,0,0,-5])
#len([0,0,0,-5,0,-5,-5,0,0,0,0,0,0,-5])
chx = [0 for i in range(nobs)]
chx = [0 for i in range(nobs)]
#for i in range(4,6): chx[i] = -5

## stride and burnin
probs_t0 = pi0_sample[40::5]
probs_inf = piinf_sample[40::5]
theta_s4str = theta_sample[40::5]
beta0_s4str = beta0_sample[40::5]
obs_s4str = obs_sample[40::5]
mat_s4str = T_sample[40::5]

f, ax = plt.subplots(1,nobs,figsize=(15,5))

for iobs,obs in  enumerate(keep_peptides):
  for rep in range(len(all_evols_o)):
    ax[iobs].plot(times, all_evols_o[rep][:,iobs]+1*chx[iobs] )
    #ax[iobs].plot(times, all_evols_o[rep][:,iobs] )


  exppts0 = [dp[iobs] for dp in exp_t0]
  exppts_inf = [dp[iobs] for dp in exp_inf]
  ax[iobs].plot(np.zeros((len(exppts0),)), exppts0,marker='o', label=obs)
  ax[iobs].plot(0, np.mean(exppts0),marker='x')
  ax[iobs].plot(np.zeros((len(exppts0),))+np.max(times), exppts_inf ,marker='o')
  ax[iobs].plot(np.max(times), np.mean(exppts_inf) ,marker='x')

  for ix, timex in enumerate(time_powers):
    exppt_time = [dp[iobs] for dp in exp_time[ix]]
    ax[iobs].plot(np.zeros((len(exppt_time),))+timex, exppt_time,marker='o')
    ax[iobs].plot(timex, np.mean(exppt_time),marker='x')

  ax[iobs].set_ylim([0, 40])
  ax[iobs].legend()

keep_peptides

"""#### Observation evolution of other peptides"""

extra_peptides = ['248-256', '9-31']
extra_frags  = [pep_to_range(peplab) for peplab in extra_peptides]
n_extra_obs = len(extra_peptides)

expdata_mod = expdata.copy()
expdata_mod.loc[(expdata['Fragment ']=="8-31"),'Fragment ']='9-31'

extra_exp_time=[]
extra_exp_t0 = get_obs(expdata_mod,lab_t0,extra_peptides)
extra_exp_inf = get_obs(expdata_mod,lab_inf,extra_peptides)

for timelag, timepower in timelags:
  tp_t0 = [tp for tl,tp in timelags if tl==lab_t0][0]
  tp_inf = [tp for tl,tp in timelags if tl==lab_inf][0]  
  if timepower> tp_t0 and timepower < tp_inf:
    extra_exp_time.append(get_obs(expdata_mod,timelag,extra_peptides))


extra_dists_all_clusters = []
for i, tr in enumerate(tr_names):
  if i in keep_states:
    extra_dists_all_clusters.append([get_distances(tr, topol, fr) for fr in extra_frags])

for i, tr2 in enumerate(tr_names2):
  if i in keep_states2:
    extra_dists_all_clusters.append([get_distances(tr2, topol2, fr) for fr in extra_frags])


def extra_obs_list_function_base(theta,beta0):
  b, xc, xh, bc, bh, b0 = theta
  obs_byclus = []
  for cluster_dists in extra_dists_all_clusters:
    obs_byclus.append([get_pf(dist, b, xc, xh, bc, bh, b0+beta0[i]) for i,dist in enumerate(cluster_dists)])
  # this is a list in which the first element are the values of obs1 for all clusters, the second is obs2 for all clusters, ...
  obs_byobs = [f for f in np.array(obs_byclus).transpose()]
  return obs_byobs


theta_s4str = theta_sample[40::5]
beta0_s4str = beta0_sample[40::5]


extra_obs_s4str = []
for i in range(len(beta0_s4str)):
  theta_ = theta_s4str[i]
  beta0_ = beta0_s4str[i]
  extra_obs_list_out = extra_obs_list_function_base( theta_,beta0_ )
  extra_obs_s4str.append(extra_obs_list_out)



#sample evolution

all_evols_extra_o = []

for i in range(len(mat_s4str)):
  prob10trjs = [probs_t0[i]]
  mat = mat_s4str[i]
  obs = extra_obs_s4str[i]
  
  for it in times[1:]:
    prt = np.matmul(prob10trjs[it-1], mat)
    prob10trjs.append(prt)
    if np.abs(prt.sum()-1)>1e-5:
      print(prt.sum())
      
  pred_time0 = forward(obs,0, prob10trjs)
  pred_time = np.vstack(pred_time0)
  all_evols_extra_o.append(pred_time)

f, ax = plt.subplots(1,n_extra_obs,figsize=(15,5))

for iobs,obs in  enumerate(extra_peptides):
  for rep in range(len(all_evols_extra_o)):
    ax[iobs].plot(times, all_evols_extra_o[rep][:,iobs])
    #ax[iobs].plot(times, all_evols_o[rep][:,iobs] )


  exppts0 = [dp[iobs] for dp in extra_exp_t0]
  exppts_inf = [dp[iobs] for dp in extra_exp_inf]
  ax[iobs].plot(np.zeros((len(exppts0),)), exppts0,marker='o', label=obs)
  ax[iobs].plot(0, np.mean(exppts0),marker='x')
  ax[iobs].plot(np.zeros((len(exppts0),))+np.max(times), exppts_inf ,marker='o')
  ax[iobs].plot(np.max(times), np.mean(exppts_inf) ,marker='x')

  for ix, timex in enumerate(time_powers):
    exppt_time = [dp[iobs] for dp in extra_exp_time[ix]]
    ax[iobs].plot(np.zeros((len(exppt_time),))+timex, exppt_time,marker='o')
    ax[iobs].plot(timex, np.mean(exppt_time),marker='x')

  ax[iobs].set_ylim([0, 40])
  ax[iobs].legend()

"""#### Probability evolution"""

f, ax = plt.subplots(1,nstates,figsize=(10,6))

for ist in  range(nstates):
  for rep in range(len(all_evols_p)):
    ax[ist].plot(times, np.vstack(all_evols_p[rep])[:,ist])

  ax[ist].set_ylim([0, 1.0])

set1

#set1 = [0,1,2,3,4,5,6,7,8]
f, ax = plt.subplots(1,len(set1),figsize=(10,4))
for ist in  set1:
  for rep in range(len(all_evols_p)):
    ax[ist].plot(times, np.vstack(all_evols_p[rep])[:,ist])

  ax[ist].set_ylim([0, 1.0])

#set2 = [9,10,11,12,13,14,15]
f, ax = plt.subplots(1,len(set2),figsize=(10,4))
for ist in set2:
  for rep in range(len(all_evols_p)):
    ax[ist-min(set2)].plot(times, np.vstack(all_evols_p[rep])[:,ist])

  ax[ist-min(set2)].set_ylim([0, 1.0])

"""#### predict TrIQ"""

print(set1)
print(set2)

i=1
print( np.round(probs_t0[i][0][set1],2) )
print( np.round(probs_t0[i][0][set2],2) )


print( np.round(probs_inf[i][0][set1],2))
print( np.round(probs_inf[i][0][set2],2))

#all_evols_p[0] is 0th replica
#all_evols_p[0][:,ist] state ist for all times
#trQI_sw = [1-np.mean(yy) for yy in [((xx/1.6)**8+1)/((xx/1.6)**12+1) for xx in trQI_dst]]


for irep in range(len(all_evols_p)):
  trq = np.sum(np.vstack([ np.vstack(all_evols_p[irep])[:,ist] * (QI_sw[ist]) for ist in range(len(QI_sw))]), axis=0)
  plt.plot((times), trq)

trq0 = []
for irep in range(len(all_evols_p)):
  trq0.append(np.sum(np.vstack([ np.vstack(all_evols_p[irep])[:,ist] * (QI_sw[ist]) for ist in range(len(QI_sw))]), axis=0))

trq_mn = np.max(np.vstack(trq0), axis=0)
plt.plot((times), trq_mn)

hdr = 'teim,trqi'
np.savetxt(f'trqi_pred-{save_chain_number:03}.csv',
        np.transpose(np.vstack((times,trq_mn))),
        header=hdr,delimiter=',',comments='')

"""#### Calculate MFPT times"""

fpts = []
for tmat in mat_s4str:
  msm1=pyemma.msm.MSM(mat_s4str[0])
  f1 = [msm1.mfpt(ia,ib) for ia in range(nstates) for ib in range(nstates) if not ia==ib]
  fpts.append(f1)

np.array(fpts).shape

"""#### Calculate Deer/opening predictions"""

tr_x3=md.load(tr_names[3], top=topol )

print(tr_x3.topology.select('residue 112 and name CA and chainid 0'))
print(tr_x3.topology.select('residue 261 and name CA  and chainid 0'))

ifr = 1
delta = tr_x3[ifr].xyz[0,1439,:] - tr_x3[ifr].xyz[0,3847,:]
np.sqrt(np.sum(delta**2))

pairs = [[1439,3847]]
md.compute_distances(tr_x3, pairs,  periodic=False)

def get_distances2(cluster_tr_filename, topology_filename, res_pair, verbose=False):

  '''
  res_pair = [('ALA',1),('ALA',2)]

  '''

  if isinstance(cluster_tr_filename, tuple):
    # only keep some frames
    tr0 = md.load(cluster_tr_filename[0], top=topology_filename )
    tr1 = tr0[cluster_tr_filename[1]]
  else:
    tr1 = md.load(cluster_tr_filename, top=topology_filename )

  reslist = [res_pair[0][1], res_pair[1][1]]
  for res in tr1.topology.residues:
    if (res.resSeq in reslist):
      if verbose: print(res)

  at1 = tr1.topology.select(f'residue {res_pair[0][1]} and name CA')
  at2 = tr1.topology.select(f'residue {res_pair[1][1]} and name CA')
  if verbose: print((at1,at2))

  
  pairs = [[at1[0],at2[0]]]
  distances = md.compute_distances(tr1, pairs,  periodic=False)
  distances_flat = distances.reshape((distances.shape[0],))
  return distances_flat

#example fro a given state
#get_distances2(tr_names[0], topol, [('ASN',112),('ASN',261)],verbose=True)

# these are already defined
#topol = '/content/drive/MyDrive/21_12_HDX/traj/Ga/Ga.pdb'
#tr_names = [f'/content/drive/MyDrive/21_12_HDX/traj/Ga/cluster_rep_0_{i}.dcd' for i in range(0,10)]
#keep_states = [0,1,2,3,4,5,6,7,8]
#topol2 = '/content/drive/MyDrive/21_12_HDX/traj/Ga_nolig/Ga-no-lig.pdb'
#tr_names2 = [f'/content/drive/MyDrive/21_12_HDX/traj/Ga_nolig/cluster_rep_2_{i}.dcd' for i in [0,1,2,3,5,6,9]]
#tr_names2[5] = (tr_names2[5], [0,1,2,4,5,8,9])
#keep_states2 = [0,1,2,3,4,5,6]

res_pair = [('ASN',112),('ASN',261)]
dist_deer = []

for i, tr in enumerate(tr_names):
  if i in keep_states:
    dist_deer.append(get_distances2(tr, topol, res_pair))

for i, tr2 in enumerate(tr_names2):
  if i in keep_states2:
    dist_deer.append(get_distances2(tr2, topol2, res_pair))

#dist_deer[0] is the distances for state 0 
state = 13
print(np.mean(dist_deer[state]))
print(dist_deer[state])

#all_evols_p[frame][time] is the nstates probabilities for each state at time 0, for frame 0 

def reweigh_dist(all_evols_p, dist_deer):
  rewei = []
  for frame, pstate_frame in enumerate(all_evols_p):
    for time, pstate_time in enumerate(pstate_frame):
      #print(pstate_time.shape)
      #pstate_time are the nstate probs for this frame, this time
      for state, dists in enumerate(dist_deer):
        w = pstate_time[0][state]
        mdist = np.mean(dists)
        rewei.append([frame, time, state, mdist, w])
        #assert abs(wtot-1)<1e-8
  
  rewei2 = np.vstack(rewei)
  print(rewei2.shape)
  return rewei2

rej_trj = reweigh_dist(all_evols_p, dist_deer)

rej_trj.shape

dist_deer

# frame, time, state, mdist, w for all frames, times, states
# representatives are less than 10 for some states
print(rej_trj.shape)
print(200*120*16)

#np.histogram( rej_trj[rej_trj[:,1]==0, 4], weights=rej_trj[rej_trj[:,1]==0, 5])
#plt.hist(rej_trj[rej_trj[:,1]==119, 3], weights=rej_trj[rej_trj[:,1]==119,4], alpha=0.5,label='time=119', bins=20)
#plt.hist(rej_trj[rej_trj[:,1]==0. , 3], weights=rej_trj[rej_trj[:,1]==  0,4], alpha=0.5,label='time=0', bins=20)
#plt.legend()
