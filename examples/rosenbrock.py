# 3d hypbrid rosenbrock data


#---------- Imports ----------#
# third party
from ctef.ctef import ctef
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#---------- Import data ----------#
X = pd.read_csv("datasets/rosenbrock.csv")
X = X.to_numpy()[:,1:]
n, p = X.shape

#---------- Parameters ----------#
k = 2
w = 2
ellipsoid_dims = [0,2]

#---------- Results ----------#
fit = ctef(X, k, w=w, ellipsoid_dimensions=ellipsoid_dims)
center = fit['center']
Lambda = fit['Lambda']


#---------- Plots ----------#
fig = plt.figure(figsize=(8, 6))

ax = fig.add_subplot(1, 1, 1, projection='3d')

ax.scatter(X[:,0], X[:,1], X[:,2], color='black', s=5, alpha=.25)
ax.scatter(center[0], center[1], center[2], color='red', s=30)

# Number of points to construct ellipse
m = 100

if k == 2:
  theta = np.linspace(0, 2*np.pi, m)
  x = np.cos(theta)
  y = np.sin(theta)
  z = np.empty(m)
  for i in range(m):
      eta = np.array([x[i],y[i]])
      [x[i],y[i],z[i]] = center + Lambda @ eta
  ax.plot(x, y, z, color='orange', alpha=1)

elif k == 3:
  phi = np.linspace(0, 2*np.pi, m)
  theta = np.linspace(0, np.pi, m)
  x = np.outer(np.cos(phi), np.sin(theta))
  y = np.outer(np.sin(phi), np.sin(theta))
  z = np.outer(np.ones_like(phi), np.cos(theta))
  for i in range(m):
    for j in range(m):
      eta = np.array([x[i,j],y[i,j],z[i,j]])
      [x[i,j],y[i,j],z[i,j]] = center + Lambda @ eta
  ax.plot_surface(x, y, z,  rstride=4, cstride=4, color='orange', alpha=.5)


plt.tight_layout()
plt.show()



