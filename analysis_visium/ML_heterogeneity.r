library(Seurat)
library(dplyr)
library(glmnet)
library(plot3D)
library(gplots)

bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))

# merge three slices
seurat_slices <- list(readRDS("data/visium/seurat/SLICE_1.seurat.rds"),
                      readRDS("data/visium/seurat/SLICE_2.seurat.rds"),
                      readRDS("data/visium/seurat/SLICE_3.seurat.rds"))
seurat_slices_spatial <- lapply(seurat_slices, function(x){
  DefaultAssay(x) <- "Spatial"
  x[['SCT']] <- NULL
  return(x)
})
seurat_slices_merged <- merge(seurat_slices_spatial[[1]], seurat_slices_spatial[-1])
spatial_coord <- do.call(rbind, lapply(seurat_slices_merged@images, function(x) x@coordinates))
rownames(spatial_coord) <- sapply(strsplit(rownames(spatial_coord), "\\."), "[", 2)
spatial_coord <- spatial_coord[colnames(seurat_slices_merged),]
spatial_coord <- cbind(spatial_coord[,3], -spatial_coord[,2], as.numeric(gsub("^S", "", seurat_slices_merged$slice)))
seurat_slices_merged[['spatial_coord']] <- CreateDimReducObject(spatial_coord, key="SPATIALCOORD_", assay = DefaultAssay(seurat_slices_merged))
seurat_slices_merged[['spatial_coord_3d']] <- CreateDimReducObject(as.matrix(setNames(data.frame(meta_slices[,c("row.sc","col.sc","ZZ")], row.names=gsub("_S", "-1_", meta_slices$X))[colnames(seurat_slices_merged),], c(paste0("SPATIALCOORD3D_",1:3)))), key="SPATIALCOORD3D_", assay = DefaultAssay(seurat_slices_merged))

## the merged data is provided as following
seurat_slices_merged <- readRDS("data/visium/seurat/merged.seurat.rds")

# normalization and dimension reduction of the merged data
seurat_slices_merged <- SCTransform(seurat_slices_merged, assay = "Spatial")
seurat_slices_merged <- RunPCA(seurat_slices_merged, npcs = 20) %>%
  RunUMAP(dims = 1:20)

# train glmnet model for single-cell cluster label prediction
seurat_sc <- readRDS("data/whole_organoids/seurat.rds")
cl_fates <- unique(data.frame(cl = seurat_sc$css_graph_res.1,
                              fate = seurat_sc$annot_region, row.names = NULL))
cl_fates <- setNames(cl_fates$fate, cl_fates$cl)

## only use a subset of single cells for training to speed out
idx_sc_train <- unlist(lapply(sort(unique(seurat$annot_region)), function(x){
  set.seed(30)
  num <- min(c(sum(seurat$annot_region == x), 400))
  sample(which(seurat$annot_region == x), replace = F)[1:num]
}))
input_genes <- intersect(VariableFeatures(seurat_slices_merged), VariableFeatures(seurat_sc))
model_cl <- cv.glmnet(x = as.matrix(t(seurat_sc@assays$RNA@data[input_genes,idx_sc_train])), # log-transformed input
                      y = seurat_sc$css_graph_res.1[idx_sc_train],
                      family = "multinomial")

## predict cluster labels and convert to cell fate identity
pred_cl <- predict(model_cl, t(as.matrix(seurat_slices_merged@assays$SCT@data[input_genes,])), type="response")[,,1]

# convert predicted cluster labels to cell fate identities
pred_cl_ident <- sapply(unique(cl_fates), function(x){
  idx <- which(cl_fates==x)
  if(length(idx)==1) return(pred_cl[,names(cl_fates)[idx]])
  return(rowSums(pred_cl[,names(cl_fates)[idx]]))
})
seurat_slices_merged$pred_cl_ident <- ifelse(apply(pred_cl_ident,1,function(x) max(x)-sort(x,decreasing=T)[2])>0.05, colnames(pred_cl_ident)[apply(pred_cl_ident,1,which.max)], "low_maxscore")

# NPC-neuron identity
DE_neuron_NPC <- read.table("data/ext/Kanton_DE_NPC_neurons.tsv", sep="\t", stringsAsFactors=F)
high_neuron <- intersect(rownames(seurat_slices_merged), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "neuron" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
high_NPC <- intersect(rownames(seurat_slices_merged), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "NPC" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
seurat_slices_merged$neuron_vs_NPC <- colMeans(seurat_slices_merged@assays$SCT@data[high_neuron,]) - colMeans(seurat_slices_merged@assays$SCT@data[high_NPC,])
seurat_slices_merged$neuron_score <- colMeans(seurat_slices_merged@assays$SCT@data[high_neuron,])
seurat_slices_merged$NPC_score <- colMeans(seurat_slices_merged@assays$SCT@data[high_NPC,])

# visualize
plotMultiFeatures(Embeddings(seurat_slices_merged, "spatial_coord")[which(seurat_slices_merged$slice=="S1"),1:2],
                  t(pred_cl_ident[which(seurat_slices_merged$slice=="S1"),c(1:5,8,6,7)])+0.1,
                  colorPal = bluered_colscheme, cex=2, pt_border=T, lwd_border = 0.5, ncol=4, cex.main = 2)

scatter3D(x = Embeddings(seurat_slices_merged, "spatial_coord_3d")[,1], y = Embeddings(seurat_slices_merged, "spatial_coord_3d")[,2], 
          z = Embeddings(seurat_slices_merged, "spatial_coord_3d")[,3],  
          colvar = as.numeric(factor(seurat_slices_merged$pred_cl_ident, levels=c("tel","dien-mesen","hind","nc","epi","mesen","low_maxscore"))),
          col = c("#932E20","#936CA8","#277C69","#D3AE2A","#1F2E3A","#CFCFCF","#FFFFFF"),
          bty = "n", pch = 19, cex = 1.3, theta = 75, phi = 30, colkey = FALSE)

freq_decov_ml <- table(seurat_slices_merged$pred_cl_ident, seurat_slices_merged$deconv_ident)[c(7,1,3,6,2,5,4),c(7,1,3,6,2,5,4)]
heatmap.2(freq_decov_ml, Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none", col=greyscale_colscheme(30), key=F, keysize=0.5, margins = c(7,7), colsep=0:7, rowsep=0:7, sepcolor = "#909090", sepwidth=c(0.01,0.01))

