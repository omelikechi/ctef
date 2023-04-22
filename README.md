# Cayley transform ellispoid fitting (CTEF)

Python implementation of Cayley tranform ellipsoid fitting (CTEF).

## Installation
To install from PyPI:
```
pip install ctef
```

## Usage
To fit an ellipsoid to data arranged in an n-by-p numpy array X (n = number of samples, p = dimension):
```python
from ctef.ctef import ctef

fit = ctef(X)
```
With all arguments (and their default values) explicitly expressed:
```python
fit = ctef(X, k=None, w=0.5, w_axis=10, ellipsoid_dimensions=None, trr_params=None)
```

## Example
Several examples are available in the examples folder. To open a Google Colab version of ellipsoid_gaussian.py, click here: <a target="_blank" href="https://colab.research.google.com/github/omelikechi/ctef/blob/main/examples/ellipsoid_gaussian.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

```python
from ctef.ctef import ctef
from examples.helpers import generate_truth, simulate_data
import matplotlib.pyplot as plt
import numpy as np

p, tau, axis_ratio, noise_level, n_samples = 2, 2, 3, .01, 50

# generate data from Ellipsoid-Gaussian model
truth = generate_truth(p, tau, axis_ratio)
X = simulate_data(n_samples, noise_level, truth)

# fit ellipsoid to data X
fit = ctef(X)
Lambda, center = fit['Lambda'], fit['center']

# plot result
n_mesh = 100
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(1, 1, 1)
ax.scatter(X[:,0], X[:,1], color='black', s=5, alpha=1)
theta = np.linspace(0, 2*np.pi, n_mesh)
x = np.cos(theta)
y = np.sin(theta)
for j in range(n_mesh):
    eta = np.array([x[j],y[j]])
    [x[j],y[j]] = center + Lambda @ eta
ax.plot(x, y, alpha=1, lw=2, color='deepskyblue')

plt.show()
```
Lambda and center yield the best fit ellipsoid $\\{\Lambda\eta+c : \lVert\eta\rVert=1\\}$ pictured below:

![example](https://user-images.githubusercontent.com/85212572/233739931-876fc8b3-467f-4499-815e-ad9f713f2c6d.png)

## Clustering
ctef_clustering.py in the ctef folder implements our ellipsoid clustering algorithm. This algorithm is tested against other clustering algorithms on two toy examples in the compare.py file. Here is its output.

![compare](https://user-images.githubusercontent.com/85212572/233740865-d516c1d9-9d43-4234-8a47-d33a4f67f052.png)
