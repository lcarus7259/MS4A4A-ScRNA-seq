library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(ggsci)
library(limma)
library(ggrepel)

# Load the dataset
MS4A.data <- Read10X(data.dir = "C:/Users/Icarus/Desktop/1.1.filtered_feature_bc_matrix/")
MS4A.meta <- read.csv('metadata.csv', row.names = 1)
# 过滤掉不在metadata中的细胞
MS4A.data <- MS4A.data[, colnames(MS4A.data) %in% rownames(MS4A.meta)]
# 创建seurat对象 
MS4A <- CreateSeuratObject(counts = MS4A.data, project = "MS4A4A", min.cells = 3, min.features = 200)
# 为meta.data添加一列样品名称
MS4A@meta.data$Sample <- MS4A.meta$Sample

# QC and selecting cells for further analysis
MS4A[["percent.mt"]] <- PercentageFeatureSet(MS4A, pattern = "^mt-")
VlnPlot(MS4A, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(MS4A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(MS4A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
MS4A <- subset(MS4A, subset = nFeature_RNA > 200 & percent.mt < 5)

# Normalizing the data
MS4A <- NormalizeData(MS4A, normalization.method = "LogNormalize", scale.factor = 10000)

# Identification of highly variable features
MS4A <- FindVariableFeatures(MS4A, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MS4A), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(MS4A)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# Scaling the data
all.genes <- rownames(MS4A)
MS4A <- ScaleData(MS4A, features = all.genes)

# Perform linear dimensional reduction
MS4A <- RunPCA(MS4A, features = VariableFeatures(object = MS4A))
# Examine and visualize PCA results a few different ways
VizDimLoadings(MS4A, dims = 1:2, reduction = "pca")
DimPlot(MS4A, reduction = "pca")
DimPlot(MS4A, reduction = "pca", label = TRUE, pt.size = 0.5, group.by = "Sample")
DimHeatmap(MS4A, dims = 1:15, cells = 500, balanced = TRUE)
# Determine the ‘dimensionality’ of the dataset
MS4A <- JackStraw(MS4A, num.replicate = 100)
MS4A <- ScoreJackStraw(MS4A, dims = 1:20)
JackStrawPlot(MS4A, dims = 1:20)
ElbowPlot(MS4A)

# Determine percent of variation associated with each PC
pct <- MS4A [["pca"]]@stdev / sum( MS4A [["pca"]]@stdev) * 100
# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
# last point where change of % of variation is more than 0.1%.
co2
# Minimum of the two calculation
pcs <- min(co1, co2)
pcs
# Create a dataframe with values
plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))
# Elbow plot to visualize 
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()

# Cluster the cells
MS4A <- FindNeighbors(MS4A, dims = 1:14)
MS4A <- FindClusters(MS4A, resolution = 0.5)
# Run non-linear dimensional reduction
# UMAP
MS4A <- RunUMAP(MS4A, dims = 1:14)
pdf(file = 'umap.pdf', width = 9, height = 6)
DimPlot(MS4A, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6)
DimPlot(MS4A, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, group.by = "Sample")
dev.off()
# 按照Sample画不同的图
pdf(file = 'umap-sample.pdf', width = 23, height = 6)
DimPlot(MS4A, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, split.by = "Sample")
dev.off()
# 统计每个cluster中每个sample的数量
sample_counts <- table(MS4A$seurat_clusters, MS4A@meta.data$Sample)
sample_counts_df <- as.data.frame(sample_counts)
colnames(sample_counts_df) <- c("Cluster", "Sample", "Count")
pdf(file = 'cluster-number.pdf', width = 12, height = 6)
ggplot(sample_counts_df, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Count", title = "Sample Counts by Cluster") +
  theme_bw()+
  theme(panel.grid.major=element_blank (), panel.grid.minor=element_blank ())+
  scale_fill_npg()
dev.off()
# T-SNE
MS4A <- RunTSNE(MS4A, dims = 1:14)
DimPlot(MS4A, reduction = "tsne")
DimPlot(MS4A, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "Sample")
# 按照Sample画不同的图
DimPlot(MS4A, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Sample")

# Finding differentially expressed features (cluster biomarkers)
MS4A.markers <- FindAllMarkers(MS4A, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MS4A.markers %>%
  group_by(cluster) %>%
  top_n(n = 50, wt = avg_log2FC) -> top50
write.csv(top50, file = "markers-1.csv")

# Data visualization for marker feature expression
features <- c("Cx3cr1","P2ry12","Tmem119","C1qa","C1qb","Apoe","Ctsb","Lyz2","Trem2","Lpl","Cst7","Spp1","Ifitm3","Ifit3","Isg15","Ccl3","Ccl4","Cd74","H2-Aa","H2-Ab1","Top2a","Mki67","Mcm2","Atp5g1","Cox6c")
pdf(file = 'expression.pdf', width = 10, height = 6)
DotPlot(MS4A, features = features) + RotatedAxis()
dev.off()

# Assigning cell type identity to clusters
new.cluster.ids <- c("HM2", "HM1", "HM1", "HM2", "IRM", "DAM2", "DAM1", "Unknown", "DAM2", "CPM", "MHCII")
names(new.cluster.ids) <- levels(MS4A)
MS4A <- RenameIdents(MS4A, new.cluster.ids)
pdf(file = 'annotation.pdf', width = 9, height = 6)
DimPlot(MS4A, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6)
dev.off()
pdf(file = 'annotation-sample.pdf', width = 23, height = 6)
DimPlot(MS4A, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 6, split.by = "Sample")
dev.off()

sample_counts <- table(MS4A@active.ident, MS4A@meta.data$Sample)
sample_counts_df <- as.data.frame(sample_counts)
write.csv(sample_counts_df, file = "1.csv")
sample_counts_df_1 <- read.csv("1.csv", row.names = 1)
colnames(sample_counts_df_1) <- c("Cluster", "Sample", "Count")
sample_counts_df_1$Cluster <- factor(sample_counts_df_1$Cluster,levels=c("HM1","HM2","DAM1","DAM2","IRM","MHCII","CPM","Unknown"))
pdf(file = 'cluster-number.pdf', width = 9, height = 6)
ggplot(sample_counts_df_1, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Count", title = "Sample Counts by Cluster") +
  theme_bw()+
  theme(panel.grid.major=element_blank (), panel.grid.minor=element_blank ())+
  scale_fill_npg()
dev.off()

library(reshape2)
library(ggalluvial)

sample_counts <- table(MS4A@meta.data$celltype, MS4A@meta.data$Sample)
sample_propotion <- apply(sample_counts, 2, function(x) x/sum(x)*100)
sample_propotion_df <- as.data.frame(sample_propotion)
sample_propotion_df$cellType <- rownames(sample_propotion_df)
sample_propotion_df <- melt(sample_propotion_df, id="cellType")
colnames(sample_propotion_df) <- c("CellType", "Sample", "Propotion")
sample_propotion_df$Sample <- factor(sample_propotion_df$Sample,levels = c("WT", "FAD_4A", "FAD"))

## 桑基图
ggplot(sample_propotion_df, aes(x =Sample, y= Propotion, fill = CellType,
                                stratum=CellType, alluvium=CellType)) +
  geom_col(width = 0.5, color='black')+
  geom_flow(width=0.5,alpha=0.4, knot.pos=0.5)+ 
  theme_classic() +
  labs(x='Sample',y = 'Ratio') +
  #coord_flip()+
  scale_fill_brewer(palette="Set3")
dev.off()

# 将IRM取出进行差异分析
IRM <- subset(MS4A, idents = "IRM")
Idents(IRM) <- IRM@meta.data$Sample

features <- c("mt-Atp6","mt-Nd4l","mt-Nd4","mt-Nd5","mt-Nd2","mt-Nd3","mt-Atp8")
pdf(file = 'IRM-ATP.pdf', width = 5, height = 5)
DotPlot(IRM, features = features) + RotatedAxis() + labs(x = NULL, y = "") + coord_flip() 
VlnPlot(IRM, features = features, ncol = 7, group.by = "Sample")
dev.off()

IRM.markers <- FindMarkers(IRM, ident.1 = "FAD_4A", ident.2 = "FAD", min.pct = 0.25, logfc.threshold = 0)
IRM.markers <- IRM.markers[order(IRM.markers$avg_log2FC, decreasing = TRUE),]
IRM.markers$expression <- ifelse(IRM.markers$avg_log2FC > 0.25 
                                  & IRM.markers$p_val < 0.05, "Up",
                                  ifelse(IRM.markers$avg_log2FC < -0.25
                                         & IRM.markers$p_val < 0.05, "Down", "None"))
IRM.markers$log10p_val <- -log10(IRM.markers$p_val)
IRM.markers$Gene <- rownames(IRM.markers)
write.csv(IRM.markers, file = 'IRM-marker.csv')

IRM.label <- bind_rows(
  IRM.markers %>%
    filter(expression == 'Up') %>%
    arrange(p_val, desc(abs(avg_log2FC))) %>%
    head(10),
  IRM.markers %>%
    filter(expression == 'Down') %>%
    arrange(p_val, desc(abs(avg_log2FC))) %>%
    head(10)
)

pdf(file = 'IRM.marker.pdf', width = 12, height = 9)
ggplot(IRM.markers, aes(avg_log2FC, log10p_val))+
  geom_point(size = 1, aes(color = expression))+
  geom_hline(yintercept=-log10(0.05), linetype=4)+
  geom_vline(xintercept=c(-1,1), linetype=4)+
  xlab(expression("log"[2]*"FoldChange"))+
  ylab(expression("-log"[10]*"(p value)")) +
  labs(title = "High vs Low", colour = "Significant")+
  # 更换颜色
  scale_color_manual(values = c("Up" = "#FF6C67","Down" = "#519BFF")) +
  # 标题居中
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_label_repel(data = IRM.label,
                   aes(avg_log2FC, log10p_val, label = Gene),
                   size = 2)
dev.off()

# GO term analysis
GO <- read.csv(file = "C:/Users/Icarus/Desktop/IRM-GO.csv")
GOBP <- subset(GO, subset = (Ontology == "BP"))[1:12,]
GO_bar <- ggplot(data = GOBP,
                 aes(x = Description, y = Count, fill = Ontology))+
  geom_bar(stat = "identity", width = 0.9, show.legend = FALSE)+
  coord_flip() + theme_bw()+ 
  scale_x_discrete(labels = function(x) str_wrap(x,width = 70))+
  labs(x = "GO terms", y = "GeneNumber", title = "Barplot of Enriched GO Terms")+
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 10), 
        plot.title = element_text(size = 14,hjust = 0.5,face = "bold"), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 10), 
        plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
ggsave(GO_bar, filename = "GO_Barplot.pdf", width = 10, height = 7)

saveRDS(MS4A, file = "MS4A_tutorial.rds")