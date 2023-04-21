# Simulate data from Ellipsoid-Gaussian model
"""
References:
  * https://arxiv.org/abs/2201.08502
"""


#---------- Imports ---------#
# third party
import numpy as np
import numpy.matlib
from numpy.random import normal, uniform
import scipy as sc
from scipy.linalg import null_space
from scipy.stats import special_ortho_group
from sklearn.decomposition import PCA


###############################################################
########## Sample from von Mises-Fisher distribution ##########
###############################################################
"""
This code is taken from
  * https://dlwhittenbury.github.io/ds-2-sampling-and-visualising-the-von-mises-fisher-distribution-in-p-dimensions.html
"""

def rand_uniform_hypersphere(N,p):
    if (p<=0) or (type(p) is not int):
        raise Exception('p must be a positive integer.')
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")

    v = np.random.normal(0,1,(N,p))

    v = np.divide(v,np.linalg.norm(v,axis=1,keepdims=True))

    return v

a = rand_uniform_hypersphere(10,3)

def rand_t_marginal(kappa,p,N=1):
    """
      rand_t_marginal(kappa,p,N=1)
      ============================

      Samples the marginal distribution of t using rejection sampling of Wood [3].

      INPUT:

          * kappa (float) - concentration        
          * p (int) - The dimension of the generated samples on the (p-1)-dimensional hypersphere.
              - p = 2 for the unit circle $\mathbb{S}^{1}$
              - p = 3 for the unit sphere $\mathbb{S}^{2}$
          Note that the (p-1)-dimensional hypersphere $\mathbb{S}^{p-1} \subset \mathbb{R}^{p}$ and the
          samples are unit vectors in $\mathbb{R}^{p}$ that lie on the sphere $\mathbb{S}^{p-1}$.
          * N (int) - number of samples

      OUTPUT:

          * samples (array of floats of shape (N,1)) - samples of the marginal distribution of t
    """

    # Check kappa >= 0 is numeric
    if (kappa < 0):
        raise Exception("kappa must be a non-negative number.")

    if (p<=0) or (type(p) is not int):
        raise Exception("p must be a positive integer.")

    # Check N>0 and is an int
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")


    # Start of algorithm
    b = (p - 1.0) / (2.0 * kappa + np.sqrt(4.0 * kappa**2 + (p - 1.0)**2 ))    
    x0 = (1.0 - b) / (1.0 + b)
    c = kappa * x0 + (p - 1.0) * np.log(1.0 - x0**2)

    samples = np.zeros((N,1))

    # Loop over number of samples
    for i in range(N):

        # Continue unil you have an acceptable sample
        while True:

            # Sample Beta distribution
            Z = np.random.beta( (p - 1.0)/2.0, (p - 1.0)/2.0 )

            # Sample Uniform distribution
            U = np.random.uniform(low=0.0,high=1.0)

            # W is essentially t
            W = (1.0 - (1.0 + b) * Z) / (1.0 - (1.0 - b) * Z)

            # Check whether to accept or reject
            if kappa * W + (p - 1.0)*np.log(1.0 - x0*W) - c >= np.log(U):

                # Accept sample
                samples[i] = W
                break

    return samples

def vmf(mu,kappa,N=1):
    """
      rand_von_mises_fisher(mu,kappa,N=1)
      ===================================

      Samples the von Mises-Fisher distribution with mean direction mu and concentration kappa.

      INPUT:

          * mu (array of floats of shape (p,1)) - mean direction. This should be a unit vector.
          * kappa (float) - concentration.
          * N (int) - Number of samples.

      OUTPUT:

          * samples (array of floats of shape (N,p)) - samples of the von Mises-Fisher distribution
          with mean direction mu and concentration kappa.
    """


    # Check that mu is a unit vector
    eps = 10**(-8) # Precision
    norm_mu = np.linalg.norm(mu)
    if abs(norm_mu - 1.0) > eps:
        raise Exception("mu must be a unit vector.")

    # Check kappa >= 0 is numeric
    if (kappa < 0):
        raise Exception("kappa must be a non-negative number.")

    # Check N>0 and is an int
    if (N<=0) or (type(N) is not int):
        raise Exception("N must be a non-zero positive integer.")

    # Dimension p
    p = len(mu)

    # Make sure that mu has a shape of px1
    mu = np.reshape(mu,(p,1))

    # Array to store samples
    samples = np.zeros((N,p))

    #  Component in the direction of mu (Nx1)
    t = rand_t_marginal(kappa,p,N)

    # Component orthogonal to mu (Nx(p-1))
    xi = rand_uniform_hypersphere(N,p-1)

    # von-Mises-Fisher samples Nxp

    # Component in the direction of mu (Nx1).
    # Note that here we are choosing an
    # intermediate mu = [1, 0, 0, 0, ..., 0] later
    # we rotate to the desired mu below
    samples[:,[0]] = t

    # Component orthogonal to mu (Nx(p-1))
    samples[:,1:] = np.matlib.repmat(np.sqrt(1 - t**2), 1, p-1) * xi

    # Rotation of samples to desired mu
    O = null_space(mu.T)
    R = np.concatenate((mu,O),axis=1)
    samples = np.dot(R,samples.T).T

    return samples


##############################################
########## Ellipsoid-Gaussian model ##########
##############################################

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
Output: n samples from von Mises-Fisher distribution with parameters mu and tau
"""
def simulate_data(n_samples, noise_level, truth):
  p = truth['p']
  Eta = vmf(truth['mu'], truth['tau'], n_samples)
  noise = 2 * noise_level * truth['major_axis_length']
  X = np.empty((n_samples, p))
  for i in range(n_samples):
    X[i,:] = truth['center'] + truth['Lambda'] @ Eta[i,:] + normal(0, noise, p)

  return X
