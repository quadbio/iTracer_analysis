library(Seurat)
library(Matrix)
library(data.tree)
library(networkD3)

# read data and create seurat object
counts <- readMM("data/microdissected/counts.mtx.gz")
meta <- read.table("data/microdissected/meta.tsv.gz", header=T, row.names=1)
rownames(counts) <- read.table("data/microdissected/features.tsv.gz", header=F, stringsAsFactors=F)[,1]
colnames(counts) <- rownames(meta)
seurat_regional <- CreateSeuratObject(counts = counts, meta.data = meta)
seurat_regional <- FindVariableFeatures(seurat_regional, nfeatures = 5000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50) %>%
  RunUMAP(dims = 1:20) %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters(resolution = 0.6)

# annotation: cell type markers expression
library(gplots)
cols_regions <- setNames(rep(c("#303030","#bdbdbd"), each=2), levels(seurat_regional$orig.ident))
cols_cl <- setNames(c("#6A51A3","#278F48","#9E9AC8","#BCBDDC","#006D2C","#87CD85","#c994c7","#A1D99B","#147E3A","#225EA8","#253494","#807DBA","#4FB264","#7FCDBB","#DADAEB","#1D91C0","#41B6C4","#38A156","#696969"), levels(seurat_regional$RNA_snn_res.0.6))
orders_cl <- c(15,7,4,3,12,1,14,17,8,6,16,13,18,10,11,2,9,5,19)
markers <- c("VIM","NES","HES1","SOX9", # NPC
             "GRIA2","DCX","NRXN1","SYP","ELAVL4", # neuron
             "RSPO3","OTX2", # diencephalon/mesencephalon
             "TTR", # choroid plexus
             "LHX9","NR2F1","NR2F2","PCDH9", # diencephalon/mesencephalon
             "TFAP2B","TFAP2A", # mid-hind boundary
             "HOXB2","HOXB5","HOXA2","HOXA3","SLC1A3", # hindbrain
             "DCN","COL1A2", # mesenchyme
             "ACTG2","TAGLN","ACTA2", # smooth muscle
             "DDIT4","HILPDA", #hypoxia
             "MKI67","CDK1","TPX2" # G2M phase
)
heatmap.2(avg_expr_regional_cl[markers,orders_cl], Rowv=NA, Colv=NA, dendrogram="none", scale="row", trace="none", key=F, keysize=0.5, ColSideColors = cols_cl[orders_cl], col = rev(brownpurple_colscheme(30)))
barplot(apply(table(seurat_regional$orig.ident, seurat_regional$RNA_snn_res.0.6)[,orders_cl],2,function(x) x/sum(x)), border=NA, space = 0, col = cols_regions, axes=F, xlab=NA, ylab=NA)

# annotation: projection to the non-spatial data
library(simspec)
css_model <- readRDS("analysis_non-spatial/css_model_nonspatial.rds")
cl_nonspatial <- read.table("non-spatial/meta.tsv.gz")$css_graph_res.1
seurat_regional <- css_project(seurat_regional, model = css_model, reduction.name = "css_proj", reduction.key = "CSSPROJ_")
proj_org_cl <- transfer_labels(data_ref = css_model$sim2profiles, data_query = Embeddings(seurat_regional, "css_proj"), label_ref = cl_nonspatial, k = 50)
seurat_regional$proj_atlas_cl <- factor(as.numeric(proj_org_cl))

# lineage tree
idx_withscar <- which(!is.na(seurat_regional$ScarMerge))[sapply(strsplit(seurat_regional$ScarMerge[!is.na(seurat_regional$ScarMerge)], "_"), function(x) length(setdiff(x, c("0:91M","0:100M")))>0)]
clone_candidates_stronger <- names(which(table(seurat_regional$GeneBarcodeMerge)>=50))
org_cells <- seurat_regional@meta.data[which(seurat_regional$GeneBarcodeMerge %in% clone_candidates_stronger), c("orig.ident","GeneBarcodeMerge","ScarMerge","RNA_snn_res.0.6")]
org_cells$ScarMerge <- paste0(org_cells$GeneBarcodeMerge, "|", org_cells$ScarMerge)
org_cells <- org_cells[base::order(org_cells$GeneBarcodeMerge, org_cells$ScarMerge, setNames(orders_cl,levels(org_cells$RNA_snn_res.0.6))[as.character(org_cells$RNA_snn_res.0.6)]),]
org_cells$cell <- rownames(org_cells)
org_cells$pathString <- paste("MD",
                              org_cells$GeneBarcodeMerge, 
                              org_cells$ScarMerge,
                              org_cells$cell,
                              sep = "/")
lineage_cells <- as.Node(org_cells)
lineage_cells_network <- ToDataFrameNetwork(lineage_cells, "name")
lineage_cells_list <- ToListExplicit(lineage_cells, unname = TRUE)

idx_node_scars <- grep("\\|", lineage_cells_network[,3])
idx_node_scarred_scars <- idx_node_scars[sapply(strsplit(unlist(lapply(strsplit(grep("\\|", lineage_cells_network[,3], value=T),"\\|"), "[", 2)), "_"), function(x) length(setdiff(x, c("0:91M","0:100M")))>0)]
idx_node_clones <- seq(min(idx_node_scars)-1)
idx_node_scarred_clones <- idx_node_clones[lineage_cells_network[idx_node_clones,3] %in% unique(sapply(strsplit(lineage_cells_network[idx_node_scarred_scars,3], "\\|"), "[", 1))]

## colored by clusters
cols <- c("#303030", cols_cl[as.character(org_cells[lineage_cells_network[,3],"RNA_snn_res.0.6"])])
cols[is.na(cols)] <- "#636363"
cols[setdiff(idx_node_scars,idx_node_scarred_scars)+1] <- "#dedede"
cols[setdiff(idx_node_clones,idx_node_scarred_clones)+1] <- "#dedede"
jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
nodeFillJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
radialNetwork(lineage_cells_list, nodeColour = nodeFillJS, nodeStroke = NA, fontSize = 0, opacity = 0.9)

## colored by regions
cols <- c("#303030", cols_regions[as.character(org_cells[lineage_cells_network[,3],"orig.ident"])])
cols[is.na(cols)] <- "#636363"
cols[setdiff(idx_node_scars,idx_node_scarred_scars)+1] <- "#dedede"
cols[setdiff(idx_node_clones,idx_node_scarred_clones)+1] <- "#dedede"
jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
nodeFillJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
radialNetwork(lineage_cells_list, nodeColour = nodeFillJS, nodeStroke = NA, fontSize = 0, opacity = 0.9)


# barcode family composition similarity
source("r_functions/BarcodeComposSim.r")
shuffled_barcodes <- shuffle_barcodes(clone_lab = seurat_regional$GeneBarcodeMerge, num_shuffle = 1000)
rownames(shuffled_barcodes) <- colnames(seurat_regional)

res_cl_clone_linkage <- count_clone_linkage_and_test(seurat_regional$RNA_snn_res.0.6, seurat_regional$GeneBarcodeMerge, shuffled_barcodes, clone_candidates)
hist(as.numeric(res_cl_clone_linkage$znorm), breaks = 20, xlab="Norm z-score", cex.lab=1.5, main = NA, lwd=2)

hcl <- hclust(as.dist(max(res_cl_clone_linkage$znorm) - res_cl_clone_linkage$znorm), method="ward.D2")
heatmap.2(res_cl_clone_linkage$znorm, Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl), trace="none", scale="none", RowSideColors=col, ColSideColors=col, col = brownpurple_colscheme(30), key=F, keysize=1, cexRow=1.2, cexCol=1.2)
cl_cls <- cutree(hcl, 3)

# scar family composition distances
source("r_functions/ScarComposDist.r")
shuffled_scars_clone <- shuffle_scars(seurat_regional$GeneBarcodeMerge, seurat_regional$ScarMerge)

clone_candidates_for_scarring <- table(seurat_regional$GeneBarcodeMerge, seurat_regional$ScarMerge)
clone_candidates_for_scarring <- clone_candidates_for_scarring[,sapply(strsplit(colnames(clone_candidates_for_scarring), "_"), function(x) sum(x == "0:91M")!=length(x) & sum(x == "0:100M")!=length(x))]
clone_candidates_for_scarring <- clone_candidates_for_scarring[rowSums(clone_candidates_for_scarring) >= 5 & rowSums(clone_candidates_for_scarring > 0) > 1,]
clone_candidates_for_scarring <- clone_candidates_for_scarring[apply(clone_candidates_for_scarring, 1, function(x)sort(x,decreasing=T)[2])/rowSums(clone_candidates_for_scarring) > 0.1,]

res_cl_scar_prop_dist <- calculate_scar_prop_dist_and_test(seurat_regional$RNA_snn_res.0.6, seurat_regional$GeneBarcodeMerge, seurat_regional$ScarMerge, shuffled_scars_clone, rownames(clone_candidates_for_scarring))

dists <- res_cl_scar_prop_dist$znorm - min(res_cl_scar_prop_dist$znorm, na.rm=T)
diag(dists) <- 0

library(gplots)
hcl <- hclust(as.dist(dists[cl_cls==1,cl_cls==1]), method="average")
heatmap.2(dists[cl_cls==1,cl_cls==1], Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl), trace="none", scale="none", col = purplegreen_colscheme(30), key=F, keysize=1, cexRow=2, cexCol=2, RowSideColors = cols_cl[cl_cls==1])
hcl <- hclust(as.dist(dists[cl_cls==2,cl_cls==2]), method="average")
heatmap.2(dists[cl_cls==2,cl_cls==2], Rowv=as.dendrogram(hcl), Colv=as.dendrogram(hcl), trace="none", scale="none", col = purplegreen_colscheme(30), key=F, keysize=1, cexRow=2, cexCol=2, RowSideColors = cols_cl[cl_cls==2])
barplot(apply(table(seurat_regional$seurat_clusters, seurat_regional$orig.ident)[hcl$order,], 1, function(x) x/sum(x)), border=NA, col=cols_regions, space=0)
