# Cayley tranform ellipsoid fitting
"""
Code implementing the Cayley tranform ellipsoid fitting algorithm.
Reference: 
"""

#---------- Imports ----------#
import numpy as np
from scipy.optimize import least_squares
from sklearn.decomposition import PCA


#---------- Loss function ----------#
# Cayley transform
"""
Input: s in R^(k*(k-1)/2)
Output: (I + S(s))^{-1}(I - S(s))
"""
def cayley_transform(s,k,triu_idx):
    S = np.zeros((k,k))
    S[triu_idx] = s
    S = S - S.T
    return np.linalg.inv(np.eye(k) + S) @ (np.eye(k) - S)

# Loss function for one point x
"""
Input: x in R^k
Output: ||diag(a) * R * (x - c)|| - 1 where
    * diag(a) is the diagonal matrix with entries a_1,...,a_k
    * c is the approximate center of the ellipsoid
    * R is the rotation corresponding to s
"""
def loss_individual(theta,x,triu_idx):
    k = x.shape[0]
    a = theta[:k]
    c = theta[k:2*k]
    s = theta[2*k:]
    R = cayley_transform(s,k,triu_idx)

    y = np.diag(a) @ (R @ (x-c))
    return np.dot(y,y) - 1

# n-dimensional array containing the loss function for each data point
def loss(theta,X,triu_idx):
    n = X.shape[0]
    output = np.empty(n)
    for i in range(n):
        output[i] = loss_individual(theta,X[i,:],triu_idx)
    return .5 * output

# gradient of loss
def grad_loss(theta,X,triu_idx):
    n, k = X.shape
    a = theta[:k]
    b = theta[k:2*k]
    s = theta[2*k:]
    R = cayley_transform(s,k,triu_idx)
    S = np.zeros((k,k))
    S[triu_idx] = s
    S = S - S.T

    a_sqr = a * a
    M = (a_sqr * R.T).T
    RtM = R.T @ M
    C1 = np.linalg.inv(np.eye(k)-S)
    C2 = np.eye(k) + R.T

    grad = np.empty((n, int(2*k + k*(k-1)/2)))
    for i in range(n):
        zi = X[i,:] - b
        ri = R @ zi
        B = C1 @ M @ np.outer(zi,zi) @ C2
        grad[i,:k] = a * ri * ri
        grad[i,k:2*k] = -ri.T @ M
        grad[i,2*k:] = (B.T - B)[np.triu_indices(k, 1)]

    return grad


#---------- Parameters for trust region reflective (trr) algorithm ----------#
"""
Input: 
    * Data matrix Y in R^(n by k) (Y is typically the PCA transformation of X)
    * weight w
    * (optional) axis weight w_axis (default is 10)
Output:
    * Initial value theta_0 = (a_0, b_0, s_0)
    * Parameter bounds (lb, ub)
"""
def trr_init(Y, w, w_axis=10):
    _, k = Y.shape

    # Min and max along principal axis
    mins = np.empty(k)
    maxs = np.empty(k)
    for i in range(k):
        mins[i] = np.min(Y[:,i])
        maxs[i] = np.max(Y[:,i])
    M = maxs - mins

    a_lb = 1/(w_axis*max(M))
    a_ub = np.inf
    c_lb = mins
    c_ub = maxs
    s_lb = -5
    s_ub = 5

    a_0 = np.ones(k)
    c_0 = (c_lb + c_ub)/2
    s_0 = np.zeros(int(k*(k-1)/2))

    c_width = (c_ub - c_lb)/2

    # set initial value
    theta_0 = np.concatenate((a_0, c_0, s_0))

    # lower and upper bounds (lb and ub) for trust region reflective algorithm
    lb = np.empty(len(theta_0))
    for i in range(len(theta_0)):
        if i < k:
            lb[i] = a_lb
        elif i < 2*k:
            lb[i] = c_0[i-k] - w*c_width[i-k]
        else:
            lb[i] = s_lb
    ub = np.empty(len(theta_0))
    for i in range(len(theta_0)):
        if i < k:
            ub[i] = a_ub
        elif i < 2*k:
            ub[i] = c_0[i-k] + w*c_width[i-k]
        else:
            ub[i] = s_ub

    return theta_0, lb, ub


#---------- Cayley transform ellipsoid fitting (CTEF) ----------#
"""
Input: 
    * X: n-by-p data matrix
    * k (default None): latent factor dimension. If 'None' k is set to p.
    * w (default 1/2): weight to apply to trust region bounds
    * w_axis (default 10): weight to apply to axis length upper bound 
    * ellipsoid_dimensions (default None): k dimensions in which to fit ellipsoid. If
    	'None' uses first k principal components
    * trr_params (default None): initial conditions for trust region algorithm. If 'None'
        applies trr_init function above.
Output:
    * center: center of fitted ellipsoid
    * Lambda: factor loading matrix of fitted ellipsoid; equal to Vk @ R^T @ A^{-1}
    * Lambda_inv: 'inverse' of Lambda; equal to A @ R @ Vk^T
    * result: output of STIR algroithm 
        ref: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html
"""

def ctef(X, k=None, w=0.5, w_axis=10, ellipsoid_dimensions=None, trr_params=None):

    n, p = X.shape
    if k == None:
    	k = p
    triu_idx = np.triu_indices(k,1)

    # pca
    avg = np.average(X,axis=0)
    X = X - avg
    pca = PCA(n_components = p)
    pca = pca.fit(X)
    V = pca.components_.T
    if k == p:
    	pass
    elif ellipsoid_dimensions == None:
        V = V[:,:k]
    else:
        V = V[:,ellipsoid_dimensions]

    Y = X @ V

    # feasible set (trust region bounds)
    if trr_params == None:
        theta_0, lb, ub = trr_init(Y, w, w_axis)
    else:
        theta_0, lb, ub = trr_params

    # STIR to minimize loss   
    result = least_squares(loss, theta_0, jac=grad_loss, args=([Y,triu_idx]), method = 'trf', bounds=(lb, ub))
    a = result.x[:k]
    c = result.x[k:2*k]
    s = result.x[2*k:]
    R = cayley_transform(s,k,triu_idx)

    return {'center': avg + V @ c, 'Lambda': V @ R.T @ np.diag(1/a), 'Lambda_inv': np.diag(a) @ R @ V.T, 'result': result}


#---------- Estimate von Mises-Fisher parameters ----------#
"""
Input: X, Lambda^{-1}=ARV' ('pseudoinverse' if k < p), center
Output: MLE of mu and tau
"""
def vmf_mle(X, Lambda_inv, center):
    n, _ = X.shape
    k, _ = Lambda_inv.shape

    X_sphere = np.empty((n,k))
    for i in range(n):
        v = Lambda_inv @ (X[i,:] - center)
        X_sphere[i,:] = v/np.linalg.norm(v)

    # MLE of direction paprameter
    mu_hat = np.average(X_sphere, axis=0)
    r_hat = np.linalg.norm(mu_hat)
    mu_hat = mu_hat/r_hat

    # MLE of concentration parameter
    tau_hat = (r_hat * (k - r_hat**2)) / (1 - r_hat**2)

    return mu_hat, tau_hat


# #---------- Draw samples ----------#
# """
# Input: center, mu, tau, ax_lengths, U, n_samples
# Output: samples from center + Lambda @ eta where Lambda = U^T @ diag(1/ax_lengths)
# """
# def samples(center, mu, tau, Lambda, n_samples):
#     samples = vmf(mu, tau, n_samples).T

#     return (Lambda @ samples).T + center


#---------- Error ----------#
"""
Input: true center, approximate center, true Lambda matrix, approximate Lambda matrix
Output: offset error and shape error
"""
def error(c_true, c_approx, Lambda_true, Lambda_approx):
    center_error = np.linalg.np.linalg.norm(c_true - c_approx)
    singular_vals = np.linalg.svd(np.linalg.inv(Lambda_approx) @ Lambda_true)[1]
    shape_error = singular_vals[0]/singular_vals[-1] - 1

    return offset_error, shape_error








