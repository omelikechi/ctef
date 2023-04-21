# Cayley tranform ellispoid fitting (CTEF)

Python implementation of Cayley tranform ellipsoid fitting (CTEF).

The main algorithm is in the ctef.py file.

The main function is ctef in the ctef.py file.

## Basic usage
Given an n-by-p data matrix X with n the number of samples and p the dimension, the ellipsoid fitting is
```python
fit = ctef(X)
```
Expanding ctef with all default arguments explicitly expressed,
```python
fit = ctef(X, k=None, w=0.5, w_axis=10, ellipsoid_dimensions=None, trr_params=None)
```
The ctef function returns the following dictionary.
```python
{'center': avg + V @ c, 'Lambda': V @ R.T @ np.diag(1/a), 'Lambda_inv': np.diag(a) @ R @ V.T, 'result': result}
```
For details see the ctef.py file and/or our paper.

## Simple example
This example is from the ellipsoid_gaussian.py file in the examples folder. See that file for parameter details.
```python
p, tau, axis_ratio, noise_level, n_samples = 2, 2, 3, .01, 50

# generate data from Ellipsoid-Gaussian model
truth = generate_truth(p, tau, axis_ratio)
X = simulate_data(n_samples, noise_level, truth)

# fit ellipsoid to data X
fit = ctef(X)
Lambda = fit['Lambda']
center = fit['center']
```
Lambda and center yield the ellipsoid $\\{\Lambda\eta+c : \lVert\eta\rVert=1\\}$ fitted to X pictured below. Note simulate_data is random, so rerunning this code will produce different data and hence a different ellipsoid of best fit.

![example](https://user-images.githubusercontent.com/85212572/233739931-876fc8b3-467f-4499-815e-ad9f713f2c6d.png)

## Clustering
ctef_clustering.py in the clustering folder implements our ellipsoid clustering algorithm. This algorithm is tested against other clustering algorithms on two toy examples in the compare.py file. Here is its output.

![compare](https://user-images.githubusercontent.com/85212572/233740865-d516c1d9-9d43-4234-8a47-d33a4f67f052.png)
