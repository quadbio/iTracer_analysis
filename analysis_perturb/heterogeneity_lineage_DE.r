library(Matrix)
library(Seurat)
library(simspec)
library(dplyr)
library(presto)
library(destiny)

source("r_functions/DifferentialExpression.r")
source("r_functions/FeaturePlots_baseR.r")
beach_colscheme <- colorRampPalette(c("#cdcdcd","#edf8b1","#7fcdbb","#41b6c4","#1d91c0","#225ea8","#0c2c84"))
bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))

seurat <- readRDS("data/perturb/whole_organoid/seurat.rds")
blacklist <- unique(c(unlist(cc.genes.updated.2019), grep("^MT-", rownames(seurat), value=T), read.table("ext/RPgenes.txt")[,1]))

seurat <- NormalizeData(seurat) %>%
  FindVariableFeatures(seurat, nfeatures = 3000)
VariableFeatures(seurat) <- setdiff(VariableFeatures(seurat), blacklist)
seurat <- ScaleData(seurat, vars.to.regress = c("G2M.Score","S.Score")) %>%
  RunPCA(npcs = 50, verbose=F) %>%
  RunUMAP(dims = 1:20)
UMAPPlot(seurat, group.by="sample") & NoAxes()

# integration with CSS
seurat <- cluster_sim_spectrum(seurat, label_tag = "sample", cluster_resolution = 1) %>%
  run_PCA(reduction="css", npcs = 20, reduction.name="css_pca", reduction.key="CSSPCA_") %>%
  RunUMAP(reduction="css_pca", dims = 1:10, reduction.key = "UMAPCSS_", reduction.name = "umap_css")
DimPlot(seurat, group.by="sample", reduction = "umap_css") & NoAxes()

# clustering (no integration)
seurat <- FindNeighbors(seurat, dims = 1:20) %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 1.5) %>%
  FindClusters(resolution = 2)

# clustering (CSS integrated)
seurat <- FindNeighbors(seurat, reduction="css_pca", dims = 1:10)
seurat[['RNA_css_snn']] <- seurat[['RNA_snn']]
seurat[['RNA_css_nn']] <- seurat[['RNA_nn']]
seurat <- FindClusters(seurat, graph.name = "RNA_css_snn", resolution = 0.1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 0.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 1.5) %>%
  FindClusters(graph.name = "RNA_css_snn", resolution = 2)

FeaturePlot(seurat, c("SOX2","HOPX","DCX","TTR","S100B","POU5F1","DCN","EPCAM","FOXG1","EMX1","DLX2","RSPO3","TCF7L2","OTX2","HOXB2","MKI67"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()
FeaturePlot(seurat, c("BCAN","ERBB3","MPZ","S100B","SOX10","SPP1","RGCC","KANK4"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()
FeaturePlot(seurat, c("GPM6A","PRPH","POU4F1","PHOX2B"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()

cols_cl <- setNames(c(c("#696969","#909090","#cdcdcd"),
                      colorRampPalette(c("#f2d7d5","#922b21"))(4),
                      c("#52BE80","#229954"),
                      "#8E44AD",
                      colorRampPalette(c("#F9E79F","#F4D03F","#D4AC0D","#B7950B"))(5)),
                    c(12,14,11,6,4,0,2,3,5,8,1,9,13,10,7))

idx_scarred <- which(lengths(lapply(strsplit(seurat$ScarMerged, "_"), setdiff, c("0:101M", NA)))>0)

plotFeature(Embeddings(seurat, "umap_css"), seurat$RNA_css_snn_res.0.5, colorPal=cols_cl, cex=0.6, pt_border = T, lwd_border=0.1, do_label=T, label_round = T, cex.label = 1.5)
plotFeature(Embeddings(seurat, "umap_css"), seurat$facs, colorPal=setNames(c("#E74C3C","#5DADE2","#A9DFBF"),c("R","B","G")), cex=0.7, pt_border = T, lwd_border=0.1, do_legend = F)
plot(Embeddings(seurat, "umap_css"), pch=16, col="#cdcdcd", cex=0.8, axes=F, xlab=NA, ylab=NA, bty="n")
points(Embeddings(seurat, "umap_css")[idx_scarred,], cex=0.9, pch=21, col="#303030", lwd=0.2, bg=setNames(c("#E74C3C","#5DADE2","#A9DFBF"),c("R","B","G"))[seurat$facs[idx_scarred]])


# focus on clusters with reasonable scarred proportion (and FOXG1+)
idx_scarred <- intersect(which(seurat$facs %in% c("R","B")), which(lengths(lapply(strsplit(seurat$ScarMerged, "_"), setdiff, c("0:101M", NA)))>0))
scar_freq <- table(seurat$RNA_css_snn_res.0.5[idx_scarred])
scar_prop <- table(seurat$RNA_css_snn_res.0.5[idx_scarred])/table(seurat$RNA_css_snn_res.0.5)

barplot(sort(scar_freq, decreasing=T), ylab="Scarred cell numbers", cex.lab=1.6, las=2, cex.names = 1.3, col = cols_cl[names(scar_freq)[order(scar_freq, decreasing=T)]])
abline(h = 100, lty=2)

#seurat <- subset(seurat, subset = RNA_css_snn_res.0.5 %in% names(which(sort(scar_freq,decreasing=T)>100)))
seurat <- readRDS("data/perturb/telen_subset/seurat.rds") # the subset data is provided and can be directed loaded

seurat <- RunUMAP(seurat, reduction="css_pca", dims = 1:10, reduction.name="umap_css", reduction.key="UMAPCSS_")
FeaturePlot(seurat, c("SOX2","DCX","FOXG1","EMX1","NKX2-1","RSPO3","TCF7L2","OTX2","MKI67"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()
FeaturePlot(seurat, c("DLK1","NES","SOX2","MKI67"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()
FeaturePlot(seurat, c("DLK1","TRIB3","COQ10B","SLC31A1","SOX2","NES","ZIC1","ZIC2","MKI67"), reduction = "umap_css", cols = beach_colscheme(30), order = T) & NoAxes() & NoLegend()

plotFeature(Embeddings(seurat, "umap_css"), seurat$facs, colorPal=setNames(c("#E74C3C","#5DADE2","#A9DFBF"),c("R","B","G")), cex=0.8, pt_border = T, lwd_border=0.2, do_legend = F)
plotFeature(Embeddings(seurat, "umap_css"), seurat$facs, colorPal=prettyrainbow_colscheme, cex=1.5, pt_border = T, lwd_border = 1, do_legend = T, legend_pos = "topleft", legend_cex = 2, emphasize=idx_scarred)
idx_scarred <- which(lengths(lapply(strsplit(seurat$ScarMerged, "_"), setdiff, c("0:101M", NA)))>0)
plotFeature(Embeddings(seurat, "umap_css"), seurat$Phase, colorPal=prettyrainbow_colscheme, cex=1.5, pt_border = T, lwd_border = 1, do_legend = T, legend_pos = "topleft", legend_cex = 2, emphasize=idx_scarred)

## clustering
seurat <- FindNeighbors(seurat, dims = 1:10, reduction = "css_pca") %>%
  FindClusters(resolution = 0.1) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1)

## pseudotime
css_pca_nocc <- apply(Embeddings(seurat, "css_pca")[,1:10], 2, function(x) residuals(lm(x ~ seurat$G2M.Score + seurat$S.Score)))
seurat[['css_pca_nocc']] <- CreateDimReducObject(css_pca_nocc, key = "CSSPCANOCC_", assay = DefaultAssay(seurat))
dm <- DiffusionMap(Embeddings(seurat, "css_pca_nocc")[,1:10], n_pcs = NA, k = 30)
seurat[['diffusionmap_nocc']] <- CreateDimReducObject(dm@eigenvectors, key="DMNOCC_", assay = DefaultAssay(seurat))

seurat$pseudotime <- rank(Embeddings(seurat,"diffusionmap_nocc")[,1])/ncol(seurat) # use ranked DC1 as pseudotime
FeaturePlot(seurat, "pseudotime", cols=bluered_colscheme(30), reduction = "umap_css") & NoAxes() & NoLegend()

## pseudotime distribution of example scar families
tail(sort(table(seurat$scar_clone[intersect(idx_scarred, which(seurat$facs=="R"))])))
tail(sort(table(seurat$scar_clone[intersect(idx_scarred, which(seurat$facs=="B"))])))
idx <- which(seurat$scar_clone %in% c("02_TSC2|B|TAGATCTAAAA|0:43M21N58M",
                                      "02_TSC2|B|TAGATCTAAAA|0:51M9D50M",
                                      "02_TSC2|B|TAGATCTAAAA|0:46M51N55M",
                                      "02_TSC2|B|TGCGCGACAGT|0:52M9D49M",
                                      "02_TSC2|R|AAGTATTAGTT|0:54M6D47M"))
pt_bin_freq <- table(ceiling(seurat$pseudotime * 20)[idx], seurat$scar_clone[idx])

barplot(pt_bin_freq, beside=T, horiz=F, col=rep(c("#2E86C1","#5DADE2","#AED6F1","#1A5276","#E74C3C"),each=20), names.arg=rep("",5), space=c(0,3))
legend("topleft", pch=21, col="#303030", pt.lwd=0.2, pt.bg=c("#2E86C1","#5DADE2","#AED6F1","#1A5276","#E74C3C"), legend = c("B|TAGATCTAAAA|0:43M21N58","B|TAGATCTAAAA|0:51M9D50M","B|TAGATCTAAAA|0:46M51N55M","B|TGCGCGACAGT|0:52M9D49M","R|AAGTATTAGTT|0:54M6D47M"), bty="n", cex=1, pt.cex = 1.5)

plot(Embeddings(seurat, "umap_css"), pch=16, col="#dedede", bty="n", xaxt="n", yaxt="n", xlab=NA, ylab=NA)
points(Embeddings(seurat, "umap_css")[which(seurat$scar_clone == "02_TSC2|B|TAGATCTAAAA|0:43M21N58M"),], pch=21, col="#303030", lwd=0.2, bg="#2E86C1", cex=1.5)
points(Embeddings(seurat, "umap_css")[which(seurat$scar_clone == "02_TSC2|B|TAGATCTAAAA|0:51M9D50M"),], pch=21, col="#303030", lwd=0.2, bg="#5DADE2", cex=1.5)
points(Embeddings(seurat, "umap_css")[which(seurat$scar_clone == "02_TSC2|B|TAGATCTAAAA|0:46M51N55M"),], pch=21, col="#303030", lwd=0.2, bg="#AED6F1", cex=1.5)
points(Embeddings(seurat, "umap_css")[which(seurat$scar_clone == "02_TSC2|B|TGCGCGACAGT|0:52M9D49M"),], pch=21, col="#303030", lwd=0.2, bg="#1A5276", cex=1.5)
points(Embeddings(seurat, "umap_css")[which(seurat$scar_clone == "02_TSC2|R|AAGTATTAGTT|0:54M6D47M"),], pch=21, col="#303030", lwd=0.2, bg="#E74C3C", cex=1.5)
legend("topleft", pch=21, col="#303030", pt.lwd=0.2, pt.bg=c("#2E86C1","#5DADE2","#AED6F1","#1A5276","#E74C3C"), legend = c("B|TAGATCTAAAA|0:43M21N58","B|TAGATCTAAAA|0:51M9D50M","B|TAGATCTAAAA|0:46M51N55M","B|TGCGCGACAGT|0:52M9D49M","R|AAGTATTAGTT|0:54M6D47M"), bty="n", cex=1.5, pt.cex = 2)


# DE analysis
## DE analysis between RFP+ and BFP+
detection_rate <- rowMeans(seurat@assays$RNA@data > 0)
gene_candidates <- names(which(detection_rate > 0.01))
DE_RFP_BFP <- ancova_group_test(expr = seurat@assays$RNA@data[gene_candidates, idx_rb],
                                group = seurat$facs_inferred[idx_rb],
                                covar = data.frame(cl = seurat$RNA_snn_res.1[idx_rb], ncounts = log(seurat$nCount_RNA[idx_rb])),
                                num_threads = 20,
                                return_coef_group = "R")

sig_DE_RFP_BFP <- DE_RFP_BFP[setdiff(rownames(DE_RFP_BFP)[which(p.adjust(DE_RFP_BFP$p_ANOVA, method="bonferroni")<0.01 & abs(DE_RFP_BFP$coef)>0.1)], c(grep("^MT-",rownames(seurat),value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1])),]
sig_DE_RFP_BFP <- sig_DE_RFP_BFP[order(sig_DE_RFP_BFP$p_ANOVA),]

## permutation for DE analysis
idx_rb <- which(seurat$facs %in% c("R","B"))
detection_rate <- rowMeans(seurat@assays$RNA@data > 0)
gene_candidates <- names(which(detection_rate > 0.01))

rand_seeds <- 1:100
permut_grp_labs <- lapply(1:100, function(i){
  set.seed(rand_seeds[i])
  grp <- sample(seurat$facs_inferred[idx_rb])
  return(grp)
})
DE_RFP_BFP_permut <- lapply(1:length(permut_grp_labs), function(i){
  grp_labs <- permut_grp_labs[[i]]
  message(paste0("start permut DE no.",i))
  ancova_group_test(expr = seurat@assays$RNA@data[gene_candidates, idx_rb],
                    group = grp_labs,
                    covar = data.frame(cl = seurat$RNA_snn_res.1[idx_rb], ncounts = log(seurat$nCount_RNA[idx_rb])),
                    num_threads = 20,
                    return_coef_group = "R")
})

## DE analysis between BFP scar clones versus other BFP scar clones
bfp_sf <- names(which(table(seurat$scar_clone[which(seurat$facs == "B")])>50))
idx_bfp <- which(seurat$facs =="B")

DE_BFP_sf <- setNames(lapply(bfp_sf, function(sf){
  grp_labs <- ifelse(seurat$scar_clone[idx_bfp] == sf,"group","bg")
  message(paste0("start BFP clone DE for ",sf))
  ancova_group_test(expr = seurat@assays$RNA@data[gene_candidates, idx_bfp],
                    group = grp_labs,
                    covar = data.frame(cl = seurat$RNA_snn_res.1[idx_bfp], ncounts = log(seurat$nCount_RNA[idx_bfp])),
                    num_threads = 20,
                    return_coef_group = "group")
}), bfp_sf)

DEG_sfs <- lapply(DE_BFP_sf, function(x) setdiff(rownames(x)[which(p.adjust(x$p_ANOVA, method="bonferroni")<0.01 & abs(x$coef)>0.1)], c(grep("^MT-",rownames(seurat),value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1])))
DEG_rfp_bfp <- setdiff(rownames(DE_RFP_BFP)[which(p.adjust(DE_RFP_BFP$p_ANOVA, method="bonferroni")<0.01 & abs(DE_RFP_BFP$coef)>0.1)], c(grep("^MT-",rownames(seurat),value=T), read.table("~/Work/databases/GeneLists/RPgenes.txt")[,1]))
DEG_specific_rfp <- setdiff(DEG_rfp_bfp, unlist(DEG_sfs))
DEG_unspecific_rfp <- intersect(DEG_rfp_bfp, unlist(DEG_sfs))

### volcano plot
par(mar=c(5,5,1,1))
plot(DE_RFP_BFP$coef,
     -log10(ifelse(DE_RFP_BFP$p_ANCOVA < 1E-20, 1E-20, DE_RFP_BFP$p_ANCOVA)),
     bty="n", xlab="Coef(RFP+)", ylab="-Log10(P-anova)", cex.lab=1.5, pch=16, col="#bdbdbd", cex=0.4)
points(sig_DE_RFP_BFP$coef,
       -log10(ifelse(sig_DE_RFP_BFP$p_ANCOVA < 1E-20, 1E-20, sig_DE_RFP_BFP$p_ANCOVA)),
       pch = 21, col="#303030", bg=ifelse(rownames(sig_DE_RFP_BFP) %in% DEG_specific_rfp, "#CB4335", "#F5B7B1"),
       cex= ifelse(rownames(sig_DE_RFP_BFP) %in% DEG_specific_rfp, 1, 0.8))
