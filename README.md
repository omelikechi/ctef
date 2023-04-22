# Cayley tranform ellispoid fitting (CTEF)

Python implementation of Cayley tranform ellipsoid fitting (CTEF).

## Installation
To install from PyPI:
```
pip install ctef
```

## Basic usage
To fit an ellipsoid to data arranged in an n-by-p matrix X (n = number of samples, p = dimension):
```python
fit = ctef(X)
```
With all default arguments explicitly expressed:
```python
fit = ctef(X, k=None, w=0.5, w_axis=10, ellipsoid_dimensions=None, trr_params=None)
```

## Simple example
From the ellipsoid_gaussian.py file in the examples folder:
```python
p, tau, axis_ratio, noise_level, n_samples = 2, 2, 3, .01, 50

# generate data from Ellipsoid-Gaussian model
truth = generate_truth(p, tau, axis_ratio)
X = simulate_data(n_samples, noise_level, truth)

# fit ellipsoid to data X
fit = ctef(X)
Lambda, center = fit['Lambda'], fit['center']
```
Lambda and center yield the best fit ellipsoid $\\{\Lambda\eta+c : \lVert\eta\rVert=1\\}$ pictured below:

![example](https://user-images.githubusercontent.com/85212572/233739931-876fc8b3-467f-4499-815e-ad9f713f2c6d.png)

## Clustering
ctef_clustering.py in the ctef folder implements our ellipsoid clustering algorithm. This algorithm is tested against other clustering algorithms on two toy examples in the compare.py file. Here is its output.

![compare](https://user-images.githubusercontent.com/85212572/233740865-d516c1d9-9d43-4234-8a47-d33a4f67f052.png)
