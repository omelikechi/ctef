from ctef.ctef import ctef
import matplotlib.pyplot as plt
from examples.helpers import generate_truth, simulate_data
import numpy as np

np.random.seed(0)

p = 2
tau = 1
axis_ratio = 2
noise_level = .01
n_samples = 50

truth = generate_truth(p, tau=tau, axis_ratio=axis_ratio)
Lambda_true = truth['Lambda']
center_true = truth['center']

X = simulate_data(n_samples, noise_level, truth)

#---------- Fit ellipsoid to X with CTEF ----------#
fit = ctef(X)
Lambda = fit['Lambda']
center = fit['center']


#---------- Plot results ----------#
n_mesh = 100

fig = plt.figure(figsize=(8, 6))

if p == 2:
  ax = fig.add_subplot(1, 1, 1)

  # data
  ax.scatter(X[:,0], X[:,1], color='black', s=5, alpha=1)

  # fit
  theta = np.linspace(0, 2*np.pi, n_mesh)
  x = np.cos(theta)
  y = np.sin(theta)
  for j in range(n_mesh):
      eta = np.array([x[j],y[j]])
      [x[j],y[j]] = center + Lambda @ eta
  ax.plot(x, y, alpha=1, lw=2, color='deepskyblue')


plt.show()



