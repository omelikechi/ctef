# Human cell cycle data
"""
Data source: https://zenodo.org/record/4525425#.ZAvqHOzMLDJ
Paper: https://www.biorxiv.org/content/10.1101/2021.02.11.430845v1.full
"""

#---------- Imports ----------#
# standard
import sys
import time

# third party
from ctef.ctef import ctef
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale


#---------- Data ----------#
phase = np.load("datasets/cell_phase.npy", allow_pickle=True)

# 40 indicates the 40 core cell cycle features have been selected
X = np.load("datasets/cell_cycle.npy")

# column 36 of X has nan values
X = np.delete(X, obj=36, axis=1)
X = scale(X)
n, p = X.shape

# pca
p = 3
pca = PCA(n_components=p)
pca.fit(X)
V = pca.components_
V = V.T
X = X @ V

# remove outliers
eps0 = .05
eps1 = .015
norm = np.empty(n)
for i in range(n):
	norm[i] = np.linalg.norm(X[i,:])
q0 = np.quantile(norm, [eps0, 1-eps0])
outlier_idx = []
for i in range(n):
	if norm[i] > q0[1]:
		outlier_idx.append(i)
phase_new = np.delete(phase, outlier_idx, axis=0)
X_outliersRemoved = np.delete(X, outlier_idx, axis=0)


#---------- Ellipsoid fit ----------#
fit = ctef(X_outliersRemoved)
# np.save('fit.npy', fit)
# fit = np.load('fit.npy', allow_pickle=True).item()

result = fit['result']
center = fit['center']
Lambda = fit['Lambda']
Lambda_inv = fit['Lambda_inv']


#---------- Plot ----------#
cdict = {'G1': 'deepskyblue', 'S': 'lime', 'G2': 'red', 'M': 'purple'}

X_ellipsoid = np.empty((n,p))
for i in range(n):
	x = Lambda_inv @ (X[i,:] - center)
	x = x/np.linalg.norm(x)
	X_ellipsoid[i,:] = Lambda @ x + center
		
fig = plt.figure(figsize=(8, 6))

# first plot
ax = fig.add_subplot(1, 2, 1, projection='3d')
s = 5

for g in np.unique(phase):
	ix = np.where(phase == g)
	if g == 'S':
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=.25)
	elif g == 'M':
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=1)
	else:
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=.5)

m = 100
phi = np.linspace(0, 2*np.pi, m)
theta = np.linspace(0, np.pi, m)
x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))
z = np.outer(np.ones_like(phi), np.cos(theta))
for i in range(m):
  for j in range(m):
    eta = np.array([x[i,j],y[i,j],z[i,j]])
    [x[i,j],y[i,j],z[i,j]] = center + Lambda @ eta
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='gray', alpha=.1)

ax.view_init(elev=10, azim=238)

# second plot
ax = fig.add_subplot(1, 2, 2, projection='3d')
s = 5

for g in np.unique(phase):
	ix = np.where(phase == g)
	if g == 'S':
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=.25)
	elif g == 'M':
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=1)
	else:
		ax.scatter(X_ellipsoid[ix,0], X_ellipsoid[ix,1], X_ellipsoid[ix,2], c = cdict[g], label = g, s = s, alpha=.5)

m = 100
phi = np.linspace(0, 2*np.pi, m)
theta = np.linspace(0, np.pi, m)
x = np.outer(np.cos(phi), np.sin(theta))
y = np.outer(np.sin(phi), np.sin(theta))
z = np.outer(np.ones_like(phi), np.cos(theta))
for i in range(m):
  for j in range(m):
    eta = np.array([x[i,j],y[i,j],z[i,j]])
    [x[i,j],y[i,j],z[i,j]] = center + Lambda @ eta
ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='gray', alpha=.1)

leg = ax.legend(loc=(.8,.5), markerscale=2)
for lh in leg.legend_handles: 
	lh.set_alpha(1)
ax.view_init(elev=5, azim=30)

plt.tight_layout()
plt.show()

















