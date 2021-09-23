library(Seurat)
library(dplyr)
library(simspec)
library(presto)
library(gplots)

bluered_colscheme <- colorRampPalette(rev(c("#d73027","#f46d43","#fdae61","#fee090","#e0f3f8","#abd9e9","#74add1","#4575b4")))
prettyrainbow_colscheme <- colorRampPalette(c("#8966A9","#3FA4D9","#3EB8B4","#B2C224","#F8CB31","#F6A22B","#EC6325","#DC3838"))


seurat_d15 <- readRDS("data/neuroepithelial/whole_organoid/seurat.rds")
seurat_d15 <- FindVariableFeatures(seurat_d15, nfeatures = 3000) %>% ScaleData() %>% RunPCA(npcs = 50, verbose=F) %>% RunUMAP(dims = 1:20)
seurat_d15 <- CellCycleScoring(seurat_d15,
                               s.features = cc.genes$s.genes,
                               g2m.features = cc.genes$g2m.genes,
                               set.ident = TRUE)
seurat_d15 <- FindNeighbors(seurat_d15) %>% FindClusters(resolution = 0.6)
FeaturePlot(seurat_d15, c("POU5F1","HES1","FOXA2","HAND1","MKI67")) & NoAxes() & NoLegend()


# focus on the neuroepithelial cells
#seurat_d15_n <- subset(seurat_d15, subset=RNA_snn_res.0.6 %in% c(0,2,8,6,7,1,4))
seurat_d15_n <- readRDS("data/neuroepithelial/neuroepithelial/seurat.rds")

seurat_d15_n <- FindVariableFeatures(seurat_d15_n, nfeatures = 3000)
VariableFeatures(seurat_d15_n) <- setdiff(VariableFeatures(seurat_d15_n), unlist(cc.genes))
seurat_d15_n <- ScaleData(seurat_d15_n, vars.to.regress=c("S.Score", "G2M.Score")) %>% RunPCA(npcs = 50) %>% RunUMAP(dims = 1:20)
seurat_d15_n$GeneBarcodeMerged[!is.na(seurat_d15_n$GeneBarcodeMerged)] <- paste0(seurat_d15_n$sample,"_",seurat_d15_n$GeneBarcodeMerged)[!is.na(seurat_d15_n$GeneBarcodeMerged)]

FeaturePlot(seurat_d15_n, c("HES1","DLK1","RSPO3","FOXG1")) & NoAxes() & NoLegend()
FeaturePlot(seurat_d15_n, c("SHH","WNT8B","BMP7","BMP4","WNT3A","WNT1","CNPY1","TRH","FGF8","FGF17","NRG3","NRG1")) & NoAxes() & NoLegend() # organizers
FeaturePlot(seurat_d15_n, c("HOXB2","OTX2","FOXG1","RSPO3")) & NoAxes() & NoLegend() # regions

seurat_d15_n <- FindNeighbors(seurat_d15_n) %>%
  FindClusters(resolution = 0.6) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 1.5) %>%
  FindClusters(resolution = 2)
seurat_d15_n@active.ident <- setNames(seurat_d15_n$RNA_snn_res.0.6, colnames(seurat_d15_n))

# css-based projection to RGs in La Manno et al.
css_model_devmouse_RG <- readRDS(file="data/ext/LaManno_css_model.rds")
seurat_devmouse_RG <- readRDS("data/ext/LaManno_RG_seurat.rds")
meta_niche <- read.table("data/ext/LaManno_meta-niche.tsv")

seurat_d15_n <- css_project(seurat_d15_n, css_model_devmouse_RG)
nn_d15_RG <- RANN::nn2(css_model_devmouse_RG$sim2profiles, Embeddings(seurat_d15_n, "css_proj"), k = 50)
proj_niche_d15 <- apply(nn_d15_RG$nn.idx, 1, function(x)
  names(which.max(table(seurat_devmouse_RG$RNA_snn_res.20[x]))))

seurat_d15_n$is_PSC <- seurat_d15_n$RNA_snn_res.0.6 %in% c("4","9")
seurat_d15_n$proj_forebrain <- ifelse(seurat_d15_n$is_PSC, FALSE, meta_niche[proj_niche_d15,"is_forebrain_snn_res.2"])
seurat_d15_n$proj_midbrain <- ifelse(seurat_d15_n$is_PSC, FALSE, meta_niche[proj_niche_d15,"is_midbrain_snn_res.2"])
seurat_d15_n$proj_hindbrain <- ifelse(seurat_d15_n$is_PSC, FALSE, meta_niche[proj_niche_d15,"is_hindbrain_snn_res.2"])
fate_idx <- ifelse(seurat_d15_n$is_PSC, "PSC", ifelse(seurat_d15_n$proj_hindbrain, "hindbrain", ifelse(seurat_d15_n$proj_midbrain, "midbrain", "forebrain")))
seurat_d15_n$proj_region <- fate_idx

mat_prop <- t(sapply(c("is_PSC","proj_forebrain","proj_midbrain","proj_hindbrain"), function(x) sapply(c(4,7,5,6,2,3,0,1), function(cl) sum(seurat_d15_n@meta.data[,x] & seurat_d15_n$RNA_snn_res.0.6 == cl)/sum(seurat_d15_n$RNA_snn_res.0.6 == cl))))
barplot(mat_prop, beside = T, names.arg = paste0("C",c(4,7,5,6,2,3,0,1)), xlab = "Cluster", ylab = "Proportion", col = prettyrainbow_colscheme(4), border = NA, cex.lab = 1.5)

# further dissect heterogeneity of the forebrain-projected cells
seurat_d15_fb <- subset(seurat_d15_n, subset = proj_region == "forebrain") %>%
  FindVariableFeatures(nfeatures = 2000) %>%
  ScaleData(vars.to.regress=c("S.Score", "G2M.Score")) %>%
  RunPCA(npcs = 20) %>%
  RunUMAP(dims = 1:10) %>%
  FindNeighbors(dims = 1:10) %>%
  FindClusters(resolution = 0.5) %>%
  FindClusters(resolution = 1) %>%
  FindClusters(resolution = 2)
avg_expr_fb_cl <- sapply(levels(seurat_d15_fb$RNA_snn_res.1), function(cl) rowMeans(seurat_d15_fb@assays$RNA@data[,which(seurat_d15_fb$RNA_snn_res.1 == cl)]))

## compare to early mouse RGs
DE_mouse_RG_ct <- wilcoxauc(seurat_devmouse_RG, "annot_snn_res.2")
DEG_mouse_RG_ct <- sort(unique(DE_mouse_RG_ct$feature[which(DE_mouse_RG_ct$logFC>log(1.2) & DE_mouse_RG_ct$auc > 0.6 & DE_mouse_RG_ct$padj < 0.01 & DE_mouse_RG_ct$pct_in > 20 & DE_mouse_RG_ct$pct_in - DE_mouse_RG_ct$pct_out > 20 & DE_mouse_RG_ct$pct_out < 20)]))

avg_expr_mouse_RG <- sapply(sort(unique(seurat_devmouse_RG$annot_snn_res.2)), function(ct) rowMeans(seurat_devmouse_RG@assays$RNA@data[,which(seurat_devmouse_RG$annot_snn_res.2 == ct)]) )
ref_mouse <- avg_expr_mouse_RG[DEG_mouse_RG_ct,]
mouse_orthologs <- read.table("data/ext/mouse_orthologs.txt", sep="\t", header=T, stringsAsFactors=F)
mouse_orthologs <- mouse_orthologs[mouse_orthologs$Mouse.homology.type == "ortholog_one2one",]
rownames(mouse_orthologs) <- make.unique(mouse_orthologs$Mouse.gene.name)
ref_mouse <- ref_mouse[rownames(ref_mouse) %in% rownames(mouse_orthologs),]
rownames(ref_mouse) <- mouse_orthologs[rownames(ref_mouse),"Gene.name"]
sim_fb_cl_RG <- ref_sim_spectrum(avg_expr_fb_cl, ref_mouse, scale = F)

seurat_d15_fb$proj_RG <- setNames(colnames(sim_fb_cl_RG)[apply(sim_fb_cl_RG, 1, which.max)], rownames(sim_fb_cl_RG))[as.character(seurat_d15_fb$RNA_snn_res.1)]

## compare to early mouse organizers
ref_mouse <- avg_expr_mouse_organizers[DEG_mouse_organizers,]
sim_fb_cl_organizers <- ref_sim_spectrum(avg_expr_fb_cl, ref_mouse, scale = F)
sim_fb_organizers <- ref_sim_spectrum(seurat_d15_fb@assays$RNA@data, ref_mouse, scale=T)
seurat_d15_fb$proj_organizer <- setNames(colnames(sim_fb_cl_organizers)[apply(sim_fb_cl_organizers, 1, which.max)],
                                         rownames(sim_fb_cl_organizers))[as.character(seurat_d15_fb$RNA_snn_res.1)]

# final annotation
seurat_d15_n$final_ident <- seurat_d15_n$proj_region
seurat_d15_n$final_ident[which(seurat_d15_n$final_ident == "forebrain")] <- paste0("forebrain: ", seurat_d15_fb$proj_RG[colnames(seurat_d15_n)[which(seurat_d15_n$final_ident == "forebrain")]])

# compare to mouse brain organizers
seurat_organizers <- readRDS("data/ext/LaManno_organizers_seurat.rds")
DE_mouse_organizers <- wilcoxauc(seurat_organizers, "organizer")
DEG_mouse_organizers <- sort(unique(DE_mouse_organizers$feature[which(DE_mouse_organizers$logFC>log(1.2) & DE_mouse_organizers$auc > 0.6 & DE_mouse_organizers$padj < 0.01 & DE_mouse_organizers$pct_in > 20 & DE_mouse_organizers$pct_in - DE_mouse_organizers$pct_out > 20 & DE_mouse_organizers$pct_out < 20)]))
avg_expr_mouse_organizers <- sapply(sort(unique(seurat_organizers$organizer)), function(ct) rowMeans(seurat_organizers@assays$RNA@data[,which(seurat_organizers$organizer == ct)]) )

avg_expr_final_ident <- sapply(c("PSC","hindbrain","midbrain",sort(unique(grep("forebrain",seurat_d15_n$final_ident, value=T)))), function(ct){
  rowMeans(seurat_d15_n@assays$RNA@data[,which(seurat_d15_n$final_ident == ct)])
})
ref_mouse <- avg_expr_mouse_organizers[DEG_mouse_organizers,]
sim2_mouse_organizers <- ref_sim_spectrum(avg_expr_final_ident, ref_mouse, scale=F)[,c(6,8,9,12,4,NA,3,11,10,5,NA,7,1,NA,2)]
heatmap.2(sim2_mouse_organizers-min(sim2_mouse_organizers,na.rm=T), # res: 800*500
          Rowv=NA, Colv=NA, dendrogram="none", scale="row", trace = "none", col=bluered_colscheme(30), key=F, keysize=0.5, margins = c(12,13), cexRow = 0.8)


# barcode and scar distribution across the projected regional identities (coarse)
idx <- which(seurat_d15_n$GeneBarcodeMerged %in% names(which(table(seurat_d15_n$GeneBarcodeMerged)>20)) & !seurat_d15_n$is_PSC)
barplot(table(factor(seurat_d15_n$proj_region)[idx], seurat_d15_n$GeneBarcodeMerged[idx])[c("PSC","forebrain","midbrain","hindbrain"),],
        las=2, border=NA, col=prettyrainbow_colscheme(4)[c(1,4,3,2)])
chisq.test(table(factor(seurat_d15_n$proj_region)[idx], seurat_d15_n$GeneBarcodeMerged[idx]))

idx <- which(seurat_d15_n$GeneBarcodeMerged == "d15_1_ACAGGTAGCTA")
idx <- idx[which(seurat_d15_n$ScarMerged[idx] %in% names(which(table(seurat_d15_n$ScarMerged[idx])>=5)))]
barplot(table(seurat_d15_n$proj_region[idx],seurat_d15_n$ScarMerged[idx])[1:2,])
