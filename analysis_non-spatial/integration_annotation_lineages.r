library(Seurat)
library(Matrix)

# read data and create seurat object
counts <- readMM("data/non-spatial/counts.mtx.gz")
meta <- read.table("data/non-spatial/meta.tsv.gz", header=T, row.names=1)
rownames(counts) <- read.table("data/non-spatial/features.tsv.gz", header=F, stringsAsFactors=F)[,1]
colnames(counts) <- rownames(meta)
seurat <- CreateSeuratObject(counts = counts, meta.data = meta)
seurat <- NormalizeData(seurat) %>% FindVariableFeatures(nfeatures = 5000)

g2m_features_more <- read.table("ext/G2M_genes.txt", stringsAsFactors=F)[,1]
seurat <- CellCycleScoring(seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), union(unlist(cc.genes), g2m_features_more))
seurat <- ScaleData(seurat, vars.to.regress = c("S.Score","G2M.Score")) %>%
  RunPCA(npcs = 50, verbose = F) %>%
  RunUMAP(dims = 1:20)

# integration by CSS representation
library(simspec)
model_css <- cluster_sim_spectrum(seurat, label_tag = "organoid", cluster_resolution = 0.6, return_seuratObj = F)

seurat <- cluster_sim_spectrum(seurat, label_tag = "organoid", cluster_resolution = 0.6)
seurat <- RunUMAP(seurat, reduction = "css", dims = 1:ncol(Embeddings(seurat, "css")), reduction.name = "umap_css", reduction.key = "UMAPCSS_")
seurat[['css_graph']] <- as.Graph(build_knn_graph(t(Embeddings(seurat, "css")), k = 50))
seurat <- FindClusters(seurat, graph.name = "css_graph", resolution = 1)

# annotation: VoxHunt
library(voxhunt)
load_aba_data('ext/voxhunt_rds')
regional_markers <- structure_markers('E13') %>%
  group_by(group) %>%
  top_n(10, auc) %>% 
  {unique(.$gene)}
vox_map <- voxel_map(
  seurat, 
  stage = 'E13', 
  group_name = 'css_graph_res.1', 
  genes_use = regional_markers
)
plot_map(vox_map, nrow=5)

# annotation: heatmap
library(gplots)
cl_orders <- c("11", "18", "8", "22", "17", "1", "16", "13", "14", "10", "15", "4", "12", "9", "26", "6", "2", "5", "3", "0", "19", "21", "20", "27", "28", "23", "7", "25", "24")
markers <- c("VIM","NES","HES1","SOX9", # NPC
             "GRIA2","DCX","NRXN1","SYP","ELAVL4", # neuron
             "GPM6A", # CNS neurons
             "FOXG1","EMX1","LHX2","KCNQ3", # dorsal telencephalon
             "LHX9","RSPO3","OTX2","CDKN1A", # diencephalon/mesencephalon
             "TTR", # choroid plexus
             "TFAP2B","TFAP2A", # mid-hind boundary
             "HOXB2","HOXB5","HOXA2","HOXA3","SLC1A3", # hindbrain
             "SIX6","VSX2", # retina
             "PRPH","SIX1","ISL1", # PNS neuron
             "FOXD3","MPZ","SOX10", # neural crest
             "EPCAM","CDH1", # epithelial
             "AIF1","ITGAM","C1QA","CD40", # microglia
             "WT1","UPK3B","KRT19", # mesothelial
             "DCN","COL1A2", # mesenchyme
             "ACTG2","TAGLN","ACTA2", # smooth muscle
             "MYOG","TNNT2","CDH15", # skeletal muscle
             "MKI67","CDK1","TPX2" # G2M phase
)
avg_expr_markers_cl <- sapply(levels(seurat$css_graph_res.1), function(cl) rowMeans(seurat@assays$RNA@data[markers,which(seurat$css_graph_res.1==cl)]))
avg_expr_markers_cl_scaled <- t(apply(avg_expr_markers_cl, 1, function(x) (x-mean(x))/(max(x)-min(x))))
heatmap.2(avg_expr_markers_cl_scaled[,cl_orders],
          Rowv=NA, Colv=NA, dendrogram="none", scale="none", trace="none",key=F, keysize=0.3, col = rev(brownpurple_colscheme(30)))

# lineage trees
library(data.tree)
library(networkD3)
cl_cols <- setNames(c("#148F77","#CCAFD8","#73C7B7","#2B9D87","#76448A","#43AB97","#8BD5C7","#CFCFCF","#FA9FB5","#DD3497","#8B5E9D","#F2D7D5","#633974","#A179B1","#966CA7","#805193","#B694C4","#D7BDE2","#922B21","#F1C40F","#F7DC6F","#D4AC0D","#C1A2CE","#EFEFEF","#909090","#AFAFAF","#A3E4D7","#212F3C","#5F6A6A"), levels(seurat$css_graph_res.1))
idx_withscar <- which(!is.na(seurat$ScarMerge))[sapply(strsplit(seurat$ScarMerge[!is.na(seurat$ScarMerge)], "_"), function(x) length(setdiff(x, c("0:91M","0:100M")))>0)]
networkd3_lineage_tree_org <- setNames(lapply(levels(seurat$organoid), function(org){
  res <- setNames(lapply(c(TRUE,FALSE), function(scarred_only){
    if (scarred_only){
      org_cells <- seurat@meta.data[intersect(idx_withscar, which(seurat$organoid == org & seurat$GeneBarcodeMerge %in% clone_candidates)), c("GeneBarcodeMerge","ScarMerge","css_graph_res.1")]
    } else{
      org_cells <- seurat@meta.data[which(seurat$organoid == org & seurat$GeneBarcodeMerge %in% clone_candidates), c("GeneBarcodeMerge","ScarMerge","css_graph_res.1")]
    }
    org_cells$ScarMerge <- paste0(org_cells$GeneBarcodeMerge, "|", org_cells$ScarMerge)
    org_cells <- org_cells[base::order(org_cells$GeneBarcodeMerge, org_cells$ScarMerge, setNames(1:length(cl_orders),cl_orders)[as.character(org_cells$css_graph_res.1),"order"]),]
    org_cells$cell <- rownames(org_cells)
    org_cells$pathString <- paste(org, 
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
    cols <- c("#303030", cl_cols[as.character(org_cells[lineage_cells_network[,3],"css_graph_res.1"])])
    cols[is.na(cols)] <- "#636363"
    cols[setdiff(idx_node_scars,idx_node_scarred_scars)+1] <- "#dedede"
    cols[setdiff(idx_node_clones,idx_node_scarred_clones)+1] <- "#dedede"
    jsarray <- paste0('["', paste(cols, collapse = '", "'), '"]')
    nodeStrokeJS <- JS(paste0('function(d, i) { return ', jsarray, '[i]; }'))
    
    res <- list(org_cells = org_cells, lineage_cells = lineage_cells, lineage_cells_list = lineage_cells_list, lineage_cells_network = lineage_cells_network, coljs = nodeStrokeJS)
    return(res)
  }), c("scarred", "all"))
}), levels(seurat$organoid))

org <- "LTv2_NEd46a"
radialNetwork(networkd3_lineage_tree_org[[org]]$all$lineage_cells_list, nodeColour = networkd3_lineage_tree_org[[org]]$all$coljs, nodeStroke = NA, fontSize = 0)
radialNetwork(networkd3_lineage_tree_org[[org]]$scarred$lineage_cells_list, nodeColour = networkd3_lineage_tree_org[[org]]$scarred$coljs, nodeStroke = NA, fontSize = 0)

# scar family composition distances at different time points
source("r_functions/ScarComposDist.r")
shuffled_scars_clone <- shuffle_scars(seurat$GeneBarcodeMerge, seurat$ScarMerge)

clone_candidates_for_scarring <- table(seurat$GeneBarcodeMerge, seurat$ScarMerge)
clone_candidates_for_scarring <- clone_candidates_for_scarring[,sapply(strsplit(colnames(clone_candidates_for_scarring), "_"), function(x) sum(x == "0:91M")!=length(x) & sum(x == "0:100M")!=length(x))]
clone_candidates_for_scarring <- clone_candidates_for_scarring[rowSums(clone_candidates_for_scarring) >= 5 & rowSums(clone_candidates_for_scarring > 0) > 1,]
clone_candidates_for_scarring <- clone_candidates_for_scarring[apply(clone_candidates_for_scarring, 1, function(x)sort(x,decreasing=T)[2])/rowSums(clone_candidates_for_scarring) > 0.1,]

scarring_day <- setNames(setNames(c(15,15,20,20,30,30,4,7,7,7),
                                  c("iDOXd15.1","iDOXd15.2","iDOXd20.2","iDOXd20.3","iDOXd30.1","iDOXd30.2","iEB32","iNE32","iNE46a","iNE46b"))[sapply(strsplit(rownames(clone_candidates_for_scarring), "_"), "[", 1)],
                         rownames(clone_candidates_for_scarring))
res_cl_scar_comp_dist_early <- calculate_scar_prop_dist_and_test(seurat$css_graph_res.1, seurat$GeneBarcodeMerge, seurat$ScarMerge, shuffled_scars_clone, rownames(clone_candidates_for_scarring)[scarring_day <= 7])
res_cl_scar_comp_dist_mid <- calculate_scar_prop_dist_and_test(seurat$css_graph_res.1, seurat$GeneBarcodeMerge, seurat$ScarMerge, shuffled_scars_clone, rownames(clone_candidates_for_scarring)[scarring_day == 15])
res_cl_scar_comp_dist_late <- calculate_scar_prop_dist_and_test(seurat$css_graph_res.1, seurat$GeneBarcodeMerge, seurat$ScarMerge, shuffled_scars_clone, rownames(clone_candidates_for_scarring)[scarring_day > 15])

zscores_early_permut <- lapply(res_cl_scar_comp_dist_early$shuffled_dists, function(obs_dist){
  shuffled_dists <- res_cl_scar_comp_dist_early$shuffled_dists
  zscores <- (obs_dist - apply(simplify2array(shuffled_dists), 1:2, mean, na.rm=T)) / apply(simplify2array(shuffled_dists), 1:2, sd, na.rm=T)
})
zscores_mid_permut <- lapply(res_cl_scar_comp_dist_mid$shuffled_dists, function(obs_dist){
  shuffled_dists <- res_cl_scar_comp_dist_mid$shuffled_dists
  zscores <- (obs_dist - apply(simplify2array(shuffled_dists), 1:2, mean, na.rm=T)) / apply(simplify2array(shuffled_dists), 1:2, sd, na.rm=T)
})
zscores_late_permut <- lapply(res_cl_scar_comp_dist_late$shuffled_dists, function(obs_dist){
  shuffled_dists <- res_cl_scar_comp_dist_late$shuffled_dists
  zscores <- (obs_dist - apply(simplify2array(shuffled_dists), 1:2, mean, na.rm=T)) / apply(simplify2array(shuffled_dists), 1:2, sd, na.rm=T)
})

den_norm_early <- density(res_cl_scar_comp_dist_early$z[!is.na(res_cl_scar_comp_dist_early$z)], bw=0.4, from = -4, to = 4)$y - rowMeans(sapply(zscores_early_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))
den_norm_mid <- density(res_cl_scar_comp_dist_mid$z[!is.na(res_cl_scar_comp_dist_mid$z)], bw=0.4, from = -4, to = 4)$y - rowMeans(sapply(zscores_mid_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))
den_norm_late <- density(res_cl_scar_comp_dist_late$z[!is.na(res_cl_scar_comp_dist_late$z)], bw=0.4, from = -4, to = 4)$y - rowMeans(sapply(zscores_late_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))

den_zscore_early_bg <- rowMeans(sapply(zscores_early_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))
den_norm_early_bg <- sapply(zscores_early_permut, function(zscores){ density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y - den_zscore_early_bg })
den_zscore_mid_bg <- rowMeans(sapply(zscores_mid_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))
den_norm_mid_bg <- sapply(zscores_mid_permut, function(zscores){ density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y - den_zscore_mid_bg })
den_zscore_late_bg <- rowMeans(sapply(zscores_late_permut, function(zscores) density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y))
den_norm_late_bg <- sapply(zscores_late_permut, function(zscores){ density(zscores[!is.na(zscores)], bw=0.4, from = -4, to = 4)$y - den_zscore_late_bg })

plot(rep(seq(-4,4,length.out=512),3), c(den_norm_early, den_norm_mid, den_norm_late), type="n", xlab="Z-score", ylab = "Norm density", main = NA, cex.lab=1.5, bty="n")
abline(h = 0, lty = 2, lwd = 1, col = "#909090")
abline(v = 0, lty = 2, lwd = 1, col = "#909090")
polygon(x = c(seq(-4,4,length.out=512), seq(4,-4,length.out=512)),
        y = c(apply(den_norm_early_bg, 1, quantile, 0.05), rev(apply(den_norm_early_bg, 1, quantile, 0.95))),
        border = "#EDBB99", col = "#EDBB99", density = 20, lwd = 1)
polygon(x = c(seq(-4,4,length.out=512), seq(4,-4,length.out=512)),
        y = c(apply(den_norm_mid_bg, 1, quantile, 0.05), rev(apply(den_norm_mid_bg, 1, quantile, 0.95))),
        border = "#FCF3CF", col = "#FCF3CF", density = 20, lwd = 1)
polygon(x = c(seq(-4,4,length.out=512), seq(4,-4,length.out=512)),
        y = c(apply(den_norm_late_bg, 1, quantile, 0.05), rev(apply(den_norm_late_bg, 1, quantile, 0.95))),
        border = "#D1F2EB", col = "#D1F2EB", density = 20, lwd = 1)

lines(seq(-4,4,length.out=512), den_norm_early, lwd = 3, col = "#F39C12")
lines(seq(-4,4,length.out=512), den_norm_mid, lwd = 3, col = "#F7DC6F")
lines(seq(-4,4,length.out=512), den_norm_late, lwd = 3, col = "#45B39D")


# barcode family composition similarity
source("r_functions/BarcodeComposSim.r")
clone_candidates <- names(which(table(seurat$GeneBarcodeMerge)>=5))
shuffled_barcodes <- shuffle_barcodes(seurat$organoid, seurat$GeneBarcodeMerge)
res_cl_barcode_sim <- count_clone_linkage_and_test(seurat$css_graph_res.1, seurat$GeneBarcodeMerge, shuffled_barcodes, clone_candidates)

hist(as.numeric(res_cl_barcode_sim$zscores_norm), breaks = 20, xlab="Norm z-score", cex.lab=1.5, main = NA, lwd=2)
abline(v = qnorm(c(0.99,0.01)), lwd=3, lty=2, col="#636363")

library(igraph)
adj_cl_clone_linkage_graph <- as.numeric(res_cl_barcode_sim$zscores_norm) > qnorm(0.99)
diag(adj_cl_clone_linkage_graph) <- 0
graph_cl_clone <- graph_from_adjacency_matrix(adj_cl_clone_linkage_graph, mode = "undirected")
components_cl_clone <- components(graph_cl_clone)
plot(graph_cl_clone)

## clustering (suppl fig 5)
### suppl fig 5b
hcl <- hclust(as.dist(max(res_cl_barcode_sim$z)-res_cl_barcode_sim$z), method="ward.D2")
gplots::heatmap.2(sign(res_cl_barcode_sim$z)*sqrt(abs(res_cl_barcode_sim$z)), Rowv = as.dendrogram(hcl), Colv = as.dendrogram(hcl), dendrogram = "both", scale = "none", key=F, keysize=0.9, trace = "none", col = rev(brownpurple_colscheme(30)), cexRow = 0.8, cexCol = 1.5)

### suppl fig 5c
prob_cl_linkage_cells <- sapply(rownames(res_cl_barcode_sim$obs_linkage), function(cl1) sapply(rownames(res_cl_barcode_sim$obs_linkage), function(cl2){
  size_cl1 <- length(intersect(which(seurat$css_graph_res.1 == cl1), which(seurat$GeneBarcodeMerge %in% clone_candidates)))
  size_cl2 <- length(intersect(which(seurat$css_graph_res.1 == cl2), which(seurat$GeneBarcodeMerge %in% clone_candidates)))
  if (cl1 == cl2)
    return(size_cl1 * size_cl2  / length(which(seurat$GeneBarcodeMerge %in% clone_candidates))^2)
  return(2 * size_cl1 * size_cl2 / length(which(seurat$GeneBarcodeMerge %in% clone_candidates))^2)
}))
expected_cl_linkage_cells <- sum(res_cl_barcode_sim$obs_linkage)/2 * prob_cl_linkage_cells
sd_cl_linkage_cells <- sqrt(sum(res_cl_barcode_sim$obs_linkage)/2 * prob_cl_linkage_cells * (1-prob_cl_linkage_cells))
diag(expected_cl_linkage_cells) <- diag(expected_cl_linkage_cells)*2
diag(sd_cl_linkage_cells) <- diag(sd_cl_linkage_cells)*2
zscores_binom <- (res_cl_barcode_sim$obs_linkage - expected_cl_linkage_cells) / sd_cl_linkage_cells

hcl <- hclust(as.dist(1-(zscores_binom-min(zscores_binom))/max(zscores_binom-min(zscores_binom))), method="ward.D2")
gplots::heatmap.2(sign(zscores_binom)*log(abs(zscores_binom)), Rowv = as.dendrogram(hcl), Colv = as.dendrogram(hcl), dendrogram = "both", scale = "none", key=F, keysize=0.9, trace = "none", col = rev(brownpurple_colscheme(30)), cexRow = 0.8, cexCol = 1.5)

### suppl fig 5d
freq_cl_orgs <- table(seurat$organoid, seurat$css_graph_res.1)
hcl <- hclust(as.dist(1-cor(freq_cl_orgs)), method="ward.D2")
gplots::heatmap.2(cor(freq_cl_orgs), Rowv = as.dendrogram(hcl), Colv = as.dendrogram(hcl), dendrogram = "both", scale = "none", key=F, keysize=0.9, trace = "none", col = rev(brownpurple_colscheme(30)), cexRow = 0.8, cexCol = 1.5)

