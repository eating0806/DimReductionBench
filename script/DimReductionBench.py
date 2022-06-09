import numpy as np
from sklearn.decomposition import PCA, FastICA, KernelPCA
from sklearn.manifold import TSNE, Isomap
import umap
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_rand_score, silhouette_score, calinski_harabasz_score, davies_bouldin_score
import warnings

warnings.filterwarnings('ignore')

def pca_py(X,n_components):
  pca = PCA(n_components = int(n_components))
  X_pca = pca.fit_transform(X)
  return X_pca

def ica_py(X,n_components):
  ica = FastICA(n_components=int(n_components))
  X_ica = ica.fit_transform(X)
  return X_ica

def poly_py(X,n_components):
  poly = KernelPCA(n_components=int(n_components), kernel='poly')
  X_poly = poly.fit_transform(X)
  return X_poly

def rbf_py(X,n_components):
  rbf = KernelPCA(n_components=int(n_components), kernel='rbf')
  X_rbf = rbf.fit_transform(X)
  return X_rbf

def tsne_py(X,n_components):
  tsne = TSNE(n_components = int(n_components))
  X_tsne = tsne.fit_transform(X)
  return X_tsne

def isomap_py(X,n_components):
  isomap = Isomap(n_components = int(n_components))
  X_isomap = isomap.fit_transform(X)
  return X_isomap

def umap_py(X,n_components):
  X_umap = umap.UMAP(n_components=int(n_components)).fit_transform(X)
  return X_umap

def cluster_py(X,n_clusters):
  X_cluster = KMeans(n_clusters=int(n_clusters)).fit(X)
  return X_cluster






