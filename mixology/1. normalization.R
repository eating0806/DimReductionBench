# load all data 
load("mRNAmix_qc.RData")
sce2_qc$group = paste(sce2_qc$H2228_prop,sce2_qc$H1975_prop,sce2_qc$HCC827_prop)
sce8_qc$group = paste(sce8_qc$H2228_prop,sce8_qc$H1975_prop,sce8_qc$HCC827_prop)
datasets <- list(
  RNAmix_CELseq2=sce2_qc,
  RNAmix_Sortseq=sce8_qc
)

# set gene filter
gene_filter = function(sce){
  keep1 = (apply(counts(sce), 1, function(x) mean(x[x>0])) > 1.1)  # average count larger than 1.1
  keep2 = (rowSums(counts(sce)>0) > 10) # expressed in more than 10 cells
  sce = sce[(keep1 & keep2), ]
  return(sce)
}

# set normalization methods
library(Linnorm)
linnorm_norm = function(sce){
  tp = system.time({
    logcounts(sce) = Linnorm(counts(sce))
  })
  method_name = "Linnorm"
  method_type = "norm"
  if (!is.null(metadata(sce)$running_time)){
    metadata(sce)$running_time = rbind(metadata(sce)$running_time, data.frame(method=method_name, method_type=method_type, time=unname(tp)[1]))
  }else{
    metadata(sce)$running_time = data.frame(method=method_name,method_type=method_type,time=unname(tp)[1])
  }
  return(sce)
}
norm_method <- list(
  Linnorm=linnorm_norm
)

# apply gene filtering
datasets = lapply(datasets,gene_filter)

# apply normalization
res <- lapply(datasets,linnorm_norm)

# save data
saveRDS(res, file="mRNAmix_norm.Rds")
