# Clustering

#---------- Imports ----------#
from ctef import ctef, trr_init
import numpy as np
from sklearn.cluster import KMeans


#---------- Residual distance ----------#
def residual(x,c,L):
	a = L @ (x - c)
	return (np.dot(a,a) - 1)**2


#---------- Clustering algorithm ----------#
def ctef_clustering(X, k, n_clusters, n_steps, w):
	n, p = X.shape

	# Initiate clusters with k-means
	kmeans = KMeans(n_clusters=n_clusters, n_init=10).fit(X)
	cluster = kmeans.labels_
	np.random.shuffle(cluster)

	cluster_size = [np.count_nonzero(i == 0) for i in range(n_clusters)]

	losses = np.ones(n_steps)

	trr_params = trr_init(X,w)

	for h in range(n_steps):
		center = np.empty((n_clusters,k))
		Lambda_inv = np.empty((n_clusters,k,k))

		# fit ellipsoid
		for j in range(n_clusters):
			Xj = np.copy(X[cluster == j])
			cluster_size[j] = len(Xj)
			if Xj.shape[0] < 3:
				print('Error: Cluster size cannot be less than 3')
				break
			
			fit = ctef(Xj, k, w, trr_params=trr_params)
			center[j,:] = fit['center']
			Lambda_inv[j,:,:] = fit['Lambda_inv']

		# check membership
		residuals = np.empty(n)
		for i in range(n):
			r = np.inf
			for j in range(n_clusters):
				r_new = residual(X[i,:], center[j,:], Lambda_inv[j,:,:])
				if r_new < r:
					r = r_new
					cluster[i] = j
			residuals[i] = r

		loss = 0
		for j in range(n_clusters):
			loss += np.sum(residuals[cluster == j])/np.count_nonzero(cluster == j)
		losses[h] = loss

	return cluster, losses, center, Lambda_inv










