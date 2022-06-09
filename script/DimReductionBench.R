
library(reticulate)
library(SingleCellExperiment)

source_python("DimReductionBench.py")

pca_r <- function(sce,n_components){
  metadata(sce)$dr_method = "pca"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_pca = pca_py(r_to_py(mat),n_components)
    reducedDim(sce,"PCA") <- mat_pca
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

ica_r <- function(sce,n_components){
  metadata(sce)$dr_method = "ica"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_ica = ica_py(r_to_py(mat),n_components)
    reducedDim(sce,"ica") <- mat_ica
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

poly_r <- function(sce,n_components){
  metadata(sce)$dr_method = "poly"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_poly = poly_py(r_to_py(mat),n_components)
    reducedDim(sce,"poly") <- mat_poly
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

rbf_r <- function(sce,n_components){
  metadata(sce)$dr_method = "rbf"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_rbf = rbf_py(r_to_py(mat),n_components)
    reducedDim(sce,"rbf") <- mat_rbf
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

tsne_r <- function(sce,n_components){
  metadata(sce)$dr_method = "tsne"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_tsne = tsne_py(r_to_py(mat),n_components)
    reducedDim(sce,"TSNE") <- mat_tsne
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

isomap_r <- function(sce,n_components){
  metadata(sce)$dr_method = "isomap"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_isomap = isomap_py(r_to_py(mat),n_components)
    reducedDim(sce,"Isomap") <- mat_isomap
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

umap_r <- function(sce,n_components){
  metadata(sce)$dr_method = "umap"
  tp = system.time({
    mat <- t(logcounts(sce))
    mat_umap = umap_py(r_to_py(mat),n_components)
    reducedDim(sce,"umap") <- mat_umap
  })
  metadata(sce)$running_time = unname(tp)[1]
  return(sce)
}

cluster_r <- function(sce,n_clusters){
  mat <- reducedDim(sce)
  mat_cluster <- cluster_py(r_to_py(mat),n_clusters)
  sce$clustering_res <- mat_cluster$labels_
  return(sce)
}

nmi_r <- function(sce){
  score <- normalized_mutual_info_score(r_to_py(sce$clustering_res),r_to_py(sce$group))
  return(round(score,3))
}

ari_r <- function(sce){
  score <- adjusted_rand_score(r_to_py(sce$clustering_res),r_to_py(sce$group))
  return(round(score,3))
}

silhouette_r <- function(sce){
  mat <- reducedDim(sce)
  score <- silhouette_score(r_to_py(mat),r_to_py(sce$clustering_res))
  return(round(score,3))
}

calinski_harabasz_r <- function(sce){
  mat <- reducedDim(sce)
  score <- calinski_harabasz_score(r_to_py(mat),r_to_py(sce$clustering_res))
  return(round(score))
}

davies_bouldin_r <- function(sce){
  mat <- reducedDim(sce)
  score <- davies_bouldin_score(r_to_py(mat),r_to_py(sce$clustering_res))
  return(round(score,3))
}

all_methods <- list(pca_r,tsne_r,isomap_r,ica_r,poly_r,rbf_r,umap_r)

names(all_methods) <- c("pca","tsne","isomap","ica","poly","rbf","umap")

call_method <- function(method){
  all_methods[method]
}

DimReduction <- function(datasets,methods,n_components){
  if (methods=="all"){
    dr_methods <- all_methods
  }
  else{
    dr_methods <- lapply(methods,call_method)
  }
  for (data in datasets) {
    for (method in dr_methods) {
      res <- method(data,n_components)
      if (exists("DRsce")) {
        DRsce = c(DRsce,list(res))
      }
      else {
        DRsce = list(res)
      }
    }
  }
  return(DRsce)
}

DimReductionEval <- function(DRsce,n_clusters) {
  Clres = lapply(DRsce,cluster_r,n_clusters=n_clusters)
  for (sce in DRsce) {
    sce_df <- data.frame(seq_method=metadata(sce)$seq_method,
                         dr_method=metadata(sce)$dr_method,
                         time=metadata(sce)$running_time)
    if (exists("DRbench")) {
      DRbench = rbind(DRbench,sce_df)
    }
    else {
      DRbench = sce_df
    }
  }
  DRbench$NMI = lapply(Clres,nmi_r)
  DRbench$ARI = lapply(Clres,ari_r)
  DRbench$silhouette = lapply(Clres,silhouette_r)
  DRbench$VRC = lapply(Clres,calinski_harabasz_r)
  DRbench$DBI = lapply(Clres,davies_bouldin_r)
  return(DRbench)
}


