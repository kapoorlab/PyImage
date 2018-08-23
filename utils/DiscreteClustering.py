from sklearn.feature_extraction import image as imp
import numpy as np
import time
from sklearn.cluster import spectral_clustering
def Discretize_Clustering(twoDimg):
      """Put clustering code"""
      graph = imp.img_to_graph(twoDimg)
      beta = 5
      eps = 1e-6
      graph.data = np.exp(-beta * graph.data / twoDimg.std()) + eps 
      N_REGIONS = 100
    
      t0 = time.time()
      labels = spectral_clustering(graph, n_clusters=N_REGIONS,
                                 assign_labels='discretize',
                                 random_state=1)
      t1 = time.time()
      labels = labels.reshape(twoDimg.shape)
      print('time taken', t1 - t0)
      
      return labels, N_REGIONS     
   