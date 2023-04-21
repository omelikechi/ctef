# Simulate data from Ellipsoid-Gaussian model
"""
References:
  * https://arxiv.org/abs/2201.08502
"""

#---------- Imports ----------#
# third party
import numpy as np
from sklearn.decomposition import PCA
from numpy.random import normal, uniform
from scipy.stats import special_ortho_group

# local
from vmf import vmf


#---------- Generate data ----------#

# generate parameters of true ellipsoid
def generate_truth(p, k=None, tau=0, axis_ratio=2, center=None, method='axis_ratio'):
  if k == None:
    k = p
  if type(center) == type(None):
    center = normal(0,10,p)
    
  mu = normal(0,1,k)
  mu = mu/np.linalg.norm(mu)

  if method == 'normal':
    Lambda = normal(0,1,(p,k))

  if method == 'axis_ratio':
    if not axis_ratio:
      print('Error: argument "axis_ratio" is missing in generate_truth function')
    U = special_ortho_group.rvs(p)
    V = special_ortho_group.rvs(k)
    sigma = np.ones(k)
    sigma[0] = axis_ratio
    sigma[1:-1] = uniform(sigma[-1], sigma[0], k-2)
    sigma /= np.prod(sigma)**(1/k)
    S = np.zeros((p,k))
    S[:k,:k] = np.diag(sigma)

    Lambda = U @ S @ V

  return {'p': p, 'k': k, 'center': center, 'Lambda': Lambda, 'mu': mu, 'tau': tau, 'major_axis_length': sigma[0]}

# von Mises-Fisher samples
"""
:Returns n samples from von Mises-Fisher distribution with parameters mu and tau
"""
def simulate_data(n_samples, noise_level, truth):
  p = truth['p']
  Eta = vmf(truth['mu'], truth['tau'], n_samples)
  noise = 2 * noise_level * truth['major_axis_length']
  X = np.empty((n_samples, p))
  for i in range(n_samples):
    X[i,:] = truth['center'] + truth['Lambda'] @ Eta[i,:] + normal(0, noise, p)

  return X

















