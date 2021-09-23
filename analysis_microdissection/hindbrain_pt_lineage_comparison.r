library(Seurat)
library(destiny)
library(WGCNA)
library(glmGamPoi)
library(splines)

seurat_regional <- readRDS("data/microdissected/seurat.rds")

# NPC-neuron scores
DE_neuron_NPC <- read.table("~/Work/brain_organoid_variability/data_integration/tab.DE_NPC_neurons.tsv", sep="\t", stringsAsFactors=F)
high_neuron <- intersect(rownames(seurat_regional), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "neuron" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
high_NPC <- intersect(rownames(seurat_regional), DE_neuron_NPC$feature[which(DE_neuron_NPC$group == "NPC" & DE_neuron_NPC$padj < 0.01 & DE_neuron_NPC$logFC > log(1.2) & DE_neuron_NPC$pct_in > 50 & DE_neuron_NPC$pct_out < 20 & DE_neuron_NPC$auc > 0.6)])
seurat_regional$neuron_vs_NPC <- colMeans(seurat_regional@assays$RNA@data[high_neuron,]) - colMeans(seurat_regional@assays$RNA@data[high_NPC,])

# hindbrain pseudotimes (only cells in CG3-1 and from the hindbrain samples)
idx_CG3_1 <- which(seurat_regional@active.ident %in% as.character(c(4,10,1,12,17,7,5,8)) & seurat_regional$orig.ident %in% c("Region2.1","Region2.2"))
barcodes_candidates <- names(which(table(seurat_regional$GeneBarcodeMerge[idx_CG3_1]) > 100))
scar_fam <- rep(NA, ncol(seurat_regional))
scar_fam[!is.na(seurat_regional$ScarMerge)] <- paste0(seurat_regional$GeneBarcodeMerge[!is.na(seurat_regional$ScarMerge)], "_", seurat_regional$ScarMerge[!is.na(seurat_regional$ScarMerge)])
scars_candidates <- names(which(table(scar_fam[idx_CG3_1]) > 50))
idx_candidate_cells <- intersect(idx_CG3_1, which(seurat_regional$GeneBarcodeMerge %in% barcodes_candidates & scar_fam %in% scars_candidates))

pca_candidate_cells <- Embeddings(seurat_regional, "pca")[idx_candidate_cells,1:20]
dm_candidate_cells <- DiffusionMap(data = pca_candidate_cells, n_pcs = NA)
pt_candidate_cells <- rank(-dm_candidate_cells@eigenvectors[,1]) / length(idx_candidate_cells)

seurat_regional$hindbrain_pt_cand <- NA
seurat_regional$hindbrain_pt_cand[idx_candidate_cells] <- pt_candidate_cells

# visualization: pseudotime distribution for different scar families
barcode_mat <- sapply(barcodes_candidates, function(x) seurat_regional$GeneBarcodeMerge[idx_candidate_cells] == x)
scar_mat <- sapply(scars_candidates, function(x) scar_fam[idx_candidate_cells] == x)

dists_lineages <- (1 - barcode_mat %*% t(barcode_mat)) + (1 - scar_mat %*% t(scar_mat))
dists_pt <- matrix(rep(pt_candidate_cells, each = length(idx_candidate_cells)), nrow = length(idx_candidate_cells)) + matrix(rep(pt_candidate_cells, length(idx_candidate_cells)), nrow = length(idx_candidate_cells))
dists <- as.dist(1 * dists_lineages + dists_pt)
hcl <- hclust(dists, method="complete")
hcl$height[which(round(hcl$height)==4)] <- 4
hcl$height[which(round(hcl$height)==3)] <- 3

scar_fam_ordered <- unique(scar_fam[idx_candidate_cells[hcl$order]])
barcode_fam_ordered <- unique(seurat_regional$GeneBarcodeMerge[idx_candidate_cells[hcl$order]])
x_coord <- 1:length(idx_candidate_cells)
x_coord <- x_coord + setNames(0:(length(scar_fam_ordered)-1) * 50, scar_fam_ordered)[scar_fam[idx_candidate_cells[hcl$order]]]
x_coord <- x_coord + setNames(0:(length(barcode_fam_ordered)-1) * 100, barcode_fam_ordered)[seurat_regional$GeneBarcodeMerge[idx_candidate_cells[hcl$order]]]
x_coord_scarfam <- tapply(x_coord, scar_fam[idx_candidate_cells[hcl$order]], list)
x_coord_bcfam <- tapply(x_coord, seurat_regional$GeneBarcodeMerge[idx_candidate_cells[hcl$order]], list)

plot(x_coord, 1-pt_candidate_cells[hcl$order], ylim=c(0,1.5), frame=F, axes=F, xlab=NA, ylab=NA, type="n")
for(i in 1:length(idx_candidate_cells))
  lines(x = c(x_coord[i],x_coord[i]), y = c(1-pt_candidate_cells[hcl$order[i]], 1.05), lty=3, col="#a0a0a0", lwd=0.3)
points(x_coord, 1-pt_candidate_cells[hcl$order], pch=21, col = "#303030", lwd=0.1, bg=setNames(c(c("#76D7C4","#16A085","#27AE60"), c("#85C1E9","#2471A3"), colorRampPalette(c("#F5B7B1","#E74C3C","#943126"))(4)), scar_fam_ordered)[scar_fam[idx_candidate_cells[hcl$order]]])
for(i in 1:length(scar_fam_ordered)){
  x_ends <- c(min(x_coord_scarfam[[scar_fam_ordered[i]]]), max(x_coord_scarfam[[scar_fam_ordered[i]]]))
  lines(x = x_ends, y = rep(1.05,2), lwd = 0.75, col = "#000000")
  lines(x = rep(mean(x_ends),2), y = c(1.05, mean(c(1.05,1.5))), lwd=0.75, col = "#000000")
}
for(i in 1:length(barcode_fam_ordered)){
  sf <- unique(scar_fam[idx_candidate_cells][which(seurat_regional$GeneBarcodeMerge[idx_candidate_cells]==barcode_fam_ordered[i])])
  x_ends <- c(min(sapply(sf, function(x) mean(x_coord_scarfam[[x]]))), max(sapply(sf, function(x) mean(x_coord_scarfam[[x]]))))
  lines(x = x_ends, y = rep(mean(c(1.05,1.5)),2), lwd = 0.75, col = "#000000")
  lines(x = rep(mean(x_coord_bcfam[[barcode_fam_ordered[i]]]),2), y = c(mean(c(1.05,1.5)), 1.5), lwd=0.75, col = "#000000")
}
lines(x = c(min(sapply(x_coord_bcfam, mean)), max(sapply(x_coord_bcfam, mean))), y = rep(1.5, 2), lwd=0.75, col = "#000000")
abline(h = 1 - cutoff_pt_NPC, lty = 2, col = "#303030", lwd=0.75)

# pt distribution clustering of scar families
pt_quant_scarfam <- sapply(tapply(pt_candidate_cells, scar_fam[idx_candidate_cells], list)[scar_fam_ordered], quantile, seq(0,1,0.05))
pt_diff_quant_scarfam <- t(sapply(2:nrow(pt_quant_scarfam), function(i) pt_quant_scarfam[i,]-pt_quant_scarfam[i-1,]))
hcl_scarfam_pt <- hclust(dist(t(pt_diff_quant_scarfam)), method="ward.D2")
hcl_scarfam_pt$labels <- c("Sf5-3","Sf5-2","Sf5-1","Sf6-2","Sf6-1","Sf4-4","Sf4-3","Sf4-2","Sf4-1")

plotDendroAndColors(hcl_scarfam_pt, c(c("#76D7C4","#16A085","#27AE60"), c("#85C1E9","#2471A3"), colorRampPalette(c("#F5B7B1","#E74C3C","#943126"))(4)))

# pt distances between different scar families, comparing those from the same barcode family and those from different ones
dists_scarfam_pt <- as.matrix(dist(t(pt_diff_quant_scarfam)))
diag(dists_scarfam_pt) <- NA
d1 <- c(as.numeric(dists_scarfam_pt[1:3,1:3]), as.numeric(dists_scarfam_pt[4:5,4:5]), as.numeric(dists_scarfam_pt[6:9,6:9]))
d2 <- dists_scarfam_pt
d2[1:3,1:3] <- NA
d2[4:5,4:5] <- NA
d2[6:9,6:9] <- NA
d2 <- as.numeric(d2)
boxplot(d1, d2, frame=F, names = c("same Bf","diff Bf"), ylab="Distance (Pt)")

# DE between scar families
idx_focused <- idx_candidate_cells[which(pt_candidate_cells < cutoff_pt_NPC)]
detect_rate <- rowMeans(seurat_regional@assays$RNA@data[,idx_focused]>0)
candidate_genes <- rownames(seurat_regional)[which(detect_rate > 0.05)]
expr <- as.matrix(exp(seurat_regional@assays$RNA@data[candidate_genes, idx_focused])-1)
barcodes <- factor(seurat_regional$GeneBarcodeMerge[idx_focused])
scars <- factor(scar_fam[idx_focused])
pt <- pt_candidate_cells[which(pt_candidate_cells < cutoff_pt_NPC)]

fit_scars <- glm_gp(expr, design = ~ scars + pt.1 + pt.2 + pt.3, col_data = data.frame(barcodes = barcodes, scars = scars, pt = ns(pt,3)))
res_scars <- test_de(fit_scars, reduced_design = ~ barcodes + pt.1 + pt.2 + pt.3)
