# Cayley tranform ellispoid fitting (CTEF)

Python implementation of Cayley tranform ellipsoid fitting (CTEF).

The main algorithm is in the ctef.py file.

The main function is ctef in the ctef.py file.

## Basic usage
Given an $n$-by-$p$ data matrix $X$ with $n$ the number of samples and $p$ the dimension, the ellipsoid fitting is
```python
fit = ctef(X)
```
Expanding ctef with all default arguments explicitly expressed,
```python
fit = ctef(X, k=None, w=0.5, w_axis=10, ellipsoid_dimensions=None, trr_params=None)
```
For details see the ctef.py file or our paper.

## Simple example
The following example is from the ellipsoid_gaussian.py file in the examples folder. See that file for parameter details.
```python
p, tau, axis_ratio, noise_level, n_samples = 2, 2, 3, .01, 50

truth = generate_truth(p, tau, axis_ratio)
X = simulate_data(n_samples, noise_level, truth)

# fit ellipsoid to data X
fit = ctef(X)
Lambda = fit['Lambda']
center = fit['center']
```
This produces the following image.



