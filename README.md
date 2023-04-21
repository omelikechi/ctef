# Cayley tranform ellispoid fitting (CTEF)

Python implementation of Cayley tranform ellipsoid fitting (CTEF).

The main algorithm is in the ctef.py file.

The main function is ctef in the ctef.py file.

## Basic usage
Given an $n$-by-$p$ data matrix $X$ with $n$ the number of samples and $p$ the dimension, the basic fitting is

```python
fit = ctef(X)
```


