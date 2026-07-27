library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(DropletUtils)
library(readr)
library(speckle)
library(MAST)
library(SeuratWrappers)
library(clusterProfiler)
library(org.Mm.eg.db)

setwd("/MSA4A4/")

dir='/MSA4A4/raw' 
samples=list.files( dir )
samples
sceList = lapply(samples,function(pro){ 
  folder=file.path(dir,pro) 
  print(pro)
  print(folder)
  print(list.files(folder))
  sce=CreateSeuratObject(Read10X(folder),project =pro )
  return(sce)
})
names(sceList) <- samples
dim(sceList[[3]])
seurat <- merge(sceList[[1]], c(unlist(sceList[-1])),add.cell.ids=names(sceList))
seurat[["percent.mt"]] <- PercentageFeatureSet(object = seurat, pattern = "^mt-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = '.origident',pt.size = 0)

#QC
filtered_seurat_list <- list()
for (seurat_name in names(sceList)) {
  seurat_object <- sceList[[seurat_name]]
  
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(object = seurat_object, pattern = "^mt-")
  seurat_object <- RunMiQC(
    seurat_object, 
    percent.mt = "percent.mt", 
    nFeature_RNA = "nFeature_RNA", 
    posterior.cutoff = 0.75, 
    model.slot = "flexmix_model", 
    model.type = "spline"
  )
  
  filtered_seurat_list[[seurat_name]] <- seurat_object
}

filtered_seurat <- merge(filtered_seurat_list[[1]], c(unlist(filtered_seurat_list[-1])),add.cell.ids=names(filtered_seurat_list))
filtered_seurat <- filtered_seurat[,filtered_seurat$miQC.keep=='keep']
filtered_seurat <- subset(filtered_seurat, subset = nFeature_RNA > 200)

#DoubletFinder
library(DoubletFinder)
dataL <- SplitObject(filtered_seurat,split.by='orig.ident')
multiplet_rates_10x <- data.frame(
  Multiplet_rate = c(0.004, 0.008, 0.016, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
  Loaded_cells = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
  Recovered_cells = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
)

for (idx in 1:length(dataL)) {
  data <- NormalizeData(dataL[[idx]])
  data <- ScaleData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, verbose = FALSE)
  data <- RunPCA(data, npcs = 30, verbose = FALSE)
  data <- RunUMAP(data, reduction = "pca", dims = 1:30)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:30)
  data <- FindClusters(data, resolution = 0.6)
  sweep.res.list <- paramSweep(data, PCs = 1:30, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- as.numeric(as.vector(bcmvn[which.max(bcmvn$BCmetric), ]$pK))
  
  recovered_cells <- nrow(data@meta.data)
  multiplet_rate <- multiplet_rates_10x$Multiplet_rate[which.min(abs(multiplet_rates_10x$Recovered_cells - recovered_cells))]
  nExp_poi <- round(multiplet_rate * recovered_cells) 
  print(paste("Sample", idx, "Multiplet Rate:", multiplet_rate)) 
  
  annotations <- data@meta.data$ClusteringResults
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  data <- doubletFinder(data, PCs = 1:30, pN = 0.25, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  
  # Save results
  dataL[[idx]]$doubFind_res = data@meta.data %>% select(contains('DF.classifications'))
  dataL[[idx]]$doubFind_score = data@meta.data %>% select(contains('pANN'))
}

filtered_seurat <- merge(dataL[[1]], c(unlist(dataL[-1])))
table(filtered_seurat$doubFind_res)
filtered_seurat <- filtered_seurat[,filtered_seurat$doubFind_res=='Singlet']

#Gene level filter: Only keeping those genes expressed in more than 10 cells
filtered_seurat <- JoinLayers(filtered_seurat)
counts <- GetAssayData(object = filtered_seurat, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)
dim(filtered_seurat)

library(biomaRt)
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes_mouse_pc <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                        filters = "biotype",
                        values = "protein_coding",
                        mart = mart)
pc_genes <- genes_mouse_pc$mgi_symbol
intersect(pc_genes,rownames(filtered_counts))
filtered_counts <- filtered_counts[intersect(pc_genes,rownames(filtered_counts)), ]
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

#Normalization, cluster and RunUMAP
filtered_seurat = SCTransform(filtered_seurat, method = "glmGamPoi", vars.to.regress = c("orig.ident"))
filtered_seurat = RunPCA(filtered_seurat, npcs = 50)
ElbowPlot(filtered_seurat, ndims = 50)
filtered_seurat = FindNeighbors(filtered_seurat, dims = 1:30, verbose = TRUE)
filtered_seurat = FindClusters(filtered_seurat, verbose = TRUE, resolution = 0.6)
filtered_seurat = RunUMAP(filtered_seurat, dims = 1:30, n.neighbors=10,min.dist=1)
DimPlot(filtered_seurat, label = TRUE)

#Marker identification
all.markers <- FindAllMarkers(
  object = filtered_seurat,
  assay = "SCT",
  only.pos = TRUE, 
  min.pct = 0.1,
  thresh.use = 0.25
)

filtered_seurat <- RenameIdents(filtered_seurat,
                                '0'='HOM',
                                '1'='HOM',
                                '2'='HOM',
                                '3'='HOM',
                                '4'='HOM',
                                '5'='HOM',
                                '6'='HOM',
                                '7'='HOM',
                                '8'='DAM2',
                                '9'='IRM',
                                '10'='DAM1',
                                '11'='HOM',
                                '12'='HOM',
                                '13'='HOM',
                                '14'='MHC-II',
                                '15'='Ccl',
                                '16'='HOM',
                                '17'='Proliferating',
                                '18'='Macrophage')


# Cell ratio calculation
Cellratio_1 <- prop.table(table(Idents(microglia), microglia$group), margin = 2) %>%
  as.data.frame() %>%
  set_names(c("Subcluster","group","ratio"))

# DEG MAST
DefaultAssay(filtered_seurat) <- 'RNA'
filtered_seurat <- NormalizeData(filtered_seurat)
Idents(filtered_seurat) = filtered_seurat$orig.ident
Idents(microglia) = microglia$orig.ident

degs_mic_4A_AD <-  FindMarkers(microglia,ident.1 = 'FAD_4A',ident.2 = 'FAD',test.use='MAST',logfc.threshold = 0,latent.vars=c('percent.mt','nFeature_RNA'),min.pct = 0.1)
degs_sub_4A_AD= lapply(unique(filtered_seurat$Subcluster), function(x){
  FindMarkers(filtered_seurat[,filtered_seurat$Subcluster==x],ident.1 = 'FAD_4A',
              ident.2 = 'FAD',test.use='MAST',logfc.threshold = 0, latent.vars=c('percent.mt','nFeature_RNA'),min.pct = 0.1)
})
names(degs_sub_4A_AD)=unique(filtered_seurat$Subcluster)

#gsea pathway analysis
mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
all_syms <- unique(unlist(lapply(degs_sub_4A_AD_ds, rownames)))
map_df <- getBM(
  attributes=c("mgi_symbol","entrezgene_id"),
  filters="mgi_symbol",
  values=all_syms,
  mart=mart
)
sym2entrez <- setNames(map_df$entrezgene_id, map_df$mgi_symbol)

gsea_list <- lapply(names(degs_sub_4A_AD), function(clust_name) {
  message("Running GSEA for subcluster: ", clust_name)
  
  df_sub <- degs_sub_4A_AD[[clust_name]]
  out_df_sorted <- df_sub[order(df_sub$avg_log2FC, decreasing = TRUE), ]
  gene_ranking <- setNames(out_df_sorted$avg_log2FC, out_df_sorted$gene)
  gene_ranking <- gene_ranking[!is.na(gene_ranking)]
  gene_ranking <- sort(gene_ranking, decreasing = TRUE)
  
  gene_entrez_raw <- gene_ranking
  names(gene_entrez_raw) <- unname(sym2entrez[names(gene_ranking)])
  
  gene_entrez_raw <- gene_entrez_raw[!is.na(names(gene_entrez_raw))]
  
  df <- data.frame(
    entrez = as.character(names(gene_entrez_raw)),
    score  = as.numeric(gene_entrez_raw),
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$score), ]
  
  df_agg <- stats::aggregate(score ~ entrez, data = df, FUN = max)
  gene_entrez_vec <- setNames(df_agg$score, df_agg$entrez)
  gene_entrez_vec <- sort(gene_entrez_vec, decreasing = TRUE)
  
  stopifnot(is.numeric(gene_entrez_vec), is.null(dim(gene_entrez_vec)))
  
  gse <- tryCatch({
    clusterProfiler::gseGO(
      geneList     = gene_entrez_vec,
      OrgDb        = org.Mm.eg.db,
      keyType      = "ENTREZID",
      ont          = "All",
      minGSSize    = 10,
      maxGSSize    = 200,
      pvalueCutoff = 1,
      eps          = 0,
      nPermSimple  = 10000,
      seed         = TRUE
    )
  }, error = function(e) {
    message("Error in subcluster ", clust_name, ": ", e$message)
    return(NULL)
  })
  
  return(gse)
})
names(gsea_list) <- names(degs_sub_4A_AD)

gsea_sym_list <- lapply(gsea_list, function(x) {
  setReadable(
    x,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID"
  )
})

