import numpy as np
import matplotlib.pyplot as plt  
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import KMeans 
from sklearn import metrics

with open("proj_m2.dat") as read_file:
	m_list = []
	for line in read_file:
		if not line:
			continue
		fields = [float(x) for x in line.split()]
		m_list.append(fields)

X = np.array(m_list)
sample_size = len(m_list)

for ncl in range(2, min(int(len(m_list)/4)+1, 11)):
	kmeans = KMeans(n_clusters=ncl)
	kmeans.fit(X)

	score = metrics.silhouette_score(X, kmeans.labels_,
	                                      metric='euclidean',
	                                      sample_size=sample_size)
	print(ncl, score)
	print(kmeans.labels_)

#	ax = Axes3D(fig)
	for i in range(len(m_list[0])):
		for j in range(i+1, len(m_list[0])):
			fig = plt.figure()
			plt.scatter(X[:,i],X[:,j], c=kmeans.labels_, cmap='rainbow')
			plt.savefig(str(ncl) + "-means_{0}_{1}.png".format(i,j))
			plt.close(fig)
