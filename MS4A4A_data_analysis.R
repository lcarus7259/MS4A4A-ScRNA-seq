setwd("/home/kxw/Project/230602_MS4A4A")

library(tidyverse)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggpubr)
library(gplots)
library(ggcharts)
library(ggsci)
library(ggrepel)
library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)
library(Nebulosa)

theme <- theme(#panel.grid.major = element_blank(), #主网格线
  #panel.grid.minor = element_blank(), #次网格线
  plot.background=element_rect(fill="white"),
  panel.border = element_rect(fill=NA), #边框
  panel.background = element_rect(fill = 'white'), #背景色
  axis.text = element_text(size = 12,family = "Arial"),
  axis.title = element_text(size = 16,family = "Arial"),
  plot.title = element_text(size = 24,hjust = 0.5,colour = "#000000",face = "bold",family = "Arial"))

mouse.data <- Read10X(data.dir = "result/1.Basic_analysis/1.1.filtered_feature_bc_matrix")
# 将gene转换为大写
mouse.data <- as.matrix(mouse.data)
rownames(mouse.data) <- toupper(rownames(mouse.data))
mouse.meta <- read.csv('metadata.csv', row.names = 1)
# 过滤掉不在metadata中的细胞
mouse.data <- mouse.data[, colnames(mouse.data) %in% rownames(mouse.meta)]

# 创建seurat对象
mouse <- CreateSeuratObject(counts = mouse.data, project = "MS4A4A")
mouse@meta.data$Sample <- mouse.meta$Sample

# QC
mouse[["percent.mt"]] <- PercentageFeatureSet(mouse, pattern = "^MT-")
VlnPlot(mouse, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouse, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

mouse <- subset(mouse, subset = percent.mt < 5)

# 标准化
mouse <- NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)
# 选择变异基因
mouse <- FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
# 标准化
all.genes <- rownames(mouse)
mouse <- ScaleData(mouse, features = all.genes)
# PCA
mouse <- RunPCA(mouse, features = VariableFeatures(object = mouse))
# 画图
VizDimLoadings(mouse, dims = 1:2, reduction = "pca")
DimPlot(mouse, reduction = "pca")
DimPlot(mouse, reduction = "pca", label = TRUE, pt.size = 0.5, group.by = "Sample")

# 挑选主成分
ElbowPlot(mouse)

# 聚类
mouse <- FindNeighbors(mouse, dims = 1:20)
mouse <- FindClusters(mouse, resolution = 0.5)

# UMAP
mouse <- RunUMAP(mouse, dims = 1:20)
# 画图
DimPlot(mouse, reduction = "umap")
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "Sample")
# 按照Sample画不同的图并保存
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Sample")

# 统计每个cluster中每个sample的数量
sample_counts <- table(mouse$seurat_clusters, mouse@meta.data$Sample)

# 绘图
sample_counts_df <- as.data.frame(sample_counts)

# 重命名数据框的列名
colnames(sample_counts_df) <- c("Cluster", "Sample", "Count")

ggplot(sample_counts_df, aes(x = Cluster, y = Count, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Count", title = "Sample Counts by Cluster") +
  theme+
  scale_fill_npg()

# T-SNE
mouse <- RunTSNE(mouse, dims = 1:20)
# 画图
DimPlot(mouse, reduction = "tsne")
DimPlot(mouse, reduction = "tsne", label = TRUE, pt.size = 0.5, group.by = "Sample")
# 按照Sample画不同的图
DimPlot(mouse, reduction = "tsne", label = TRUE, pt.size = 0.5, split.by = "Sample")

VlnPlot(mouse, features = c("MS4A4A"))
FeaturePlot(mouse, features = c("MS4A4A"),split.by = "Sample")

mouse.markers <- FindAllMarkers(mouse, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 保存
write.csv(mouse.markers, file = "markers.csv")

# 求子集
FAD <- subset(mouse,Sample == "FAD")
FAD.4A <- subset(mouse,Sample == "FAD_4A")

mouse.makers.fad <- FindAllMarkers(FAD, min.pct = 0.25, logfc.threshold = 0.25)
mouse.makers.fad.4a <- FindAllMarkers(FAD.4A, min.pct = 0.25, logfc.threshold = 0.25)

# 按照cluster进行计数
sample_deg_df <- data.frame(Cluster = integer(), FAD = integer(), FAD_4A = integer())
fad.deg <- mouse.makers.fad %>%
    group_by(cluster) %>%
    summarise(n = n())
fad.4a.deg <- mouse.makers.fad.4a %>%
    group_by(cluster) %>%
    summarise(n = n())

# 合并
sample_deg_df <- merge(fad.deg, fad.4a.deg, by = "cluster", all = TRUE)
# 重命名列名
colnames(sample_deg_df) <- c("Cluster", "FAD", "FAD_4A")

# 换堆叠方式
sample_deg_df <- cbind(sample_deg_df[1], stack(sample_deg_df[,2:3]))
colnames(sample_deg_df) <- c("Cluster", "Counts", "Sample")

ggplot(sample_deg_df, aes(x = Cluster, y = Counts, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "DEGs", title = "DEGs by Cluster") +
  theme+
  scale_fill_npg()

new.cluster.ids <- c("unknown1","HM","DAM2","DAM2","ARM","DAM1","unknown2","DAM2","CPM","ARM")
names(new.cluster.ids) <- levels(mouse)
mouse <- RenameIdents(mouse, new.cluster.ids)
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5,label.size = 8) + NoLegend()+
  stat_density_2d(mapping=aes(x=UMAP_1, y=UMAP_2),
                  linemitre = 20)
DimPlot(mouse, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Sample")+
  stat_density_2d(mapping=aes(x=UMAP_1, y=UMAP_2),
                  linemitre = 20)



mouse@meta.data$annotation <- Idents(mouse)

features = c("CX3CR1","P2RY12","TMEM119","APOE","LYZ2","CTSB","TREM2","CST7","LPL","CD74","H2-AB1","MCM2","TOP2A")
fplot <- FeaturePlot(mouse, features = features)

ggsave("featurePlot.png",fplot,width = 22, height = 20)

Idents(mouse) <- mouse@meta.data$seurat_clusters
VlnPlot(mouse, features = features,stack = TRUE, flip = TRUE)

# 将DAM1,DAM2分别取出进行差异分析
DAM1 <- subset(mouse, annotation == "DAM1")
DAM2 <- subset(mouse, annotation == "DAM2")

Idents(DAM1) <- DAM1@meta.data$annotation
DimPlot(DAM1, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Sample")+
  stat_density_2d(mapping=aes(x=UMAP_1, y=UMAP_2),
                  linemitre = 20)

Idents(DAM2) <- DAM2@meta.data$annotation
DimPlot(DAM2, reduction = "umap", label = TRUE, pt.size = 0.5, split.by = "Sample")+
  stat_density_2d(mapping=aes(x=UMAP_1, y=UMAP_2),
                  linemitre = 20)

Idents(DAM1) <- DAM1@meta.data$Sample

DAM1.markers <- FindMarkers(DAM1, ident.1 = "FAD_4A", ident.2 = "FAD",logfc.threshold = 0, min.pct = 0.25)
DAM1.markers <- DAM1.markers[order(DAM1.markers$avg_log2FC),]

DAM1.markers$expression <- ifelse(DAM1.markers$avg_log2FC > 0 
                                       & DAM1.markers$p_val < 0.05, "Up",
                                         ifelse(DAM1.markers$avg_log2FC < 0
                                                & DAM1.markers$p_val < 0.05, "Down", "None"))
DAM1.markers$log10p_val <- -log10(DAM1.markers$p_val)
DAM1.markers$Gene <- rownames(DAM1.markers)

DAM1.markers <- DAM1.markers[order(DAM1.markers$avg_log2FC,decreasing = TRUE),]

data <- bind_rows(
  DAM1.markers %>%
    filter(expression == 'Up') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10),
  DAM1.markers %>%
    filter(expression == 'Down') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10)
)

write.csv(DAM1.markers,file = "DAM1_marker_0.csv")

DAM1.degPlot <- ggplot(DAM1.markers, aes(avg_log2FC,log10p_val))+
  geom_point(size = 1,aes(color = expression))+
  #scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  xlab(expression("log"[2]*"FoldChange"))+
  ylab(expression("-log"[10]*"(p value)")) +
  labs(title = "High vs Low",colour = "Significant")+
  # 更换颜色
  scale_color_manual(values = c("Up" = "#FF6C67","Down" = "#519BFF")) +
  # 标题居中
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_label_repel(data = data,
                   aes(avg_log2FC, log10p_val, label = Gene),
                   size = 2)

Idents(DAM2) <- DAM2@meta.data$Sample

DAM2.markers <- FindMarkers(DAM2, ident.1 = "FAD_4A", ident.2 = "FAD",logfc.threshold = 0, min.pct = 0.25)
DAM2.markers <- DAM2.markers[order(DAM2.markers$avg_log2FC),]

DAM2.markers$expression <- ifelse(DAM2.markers$avg_log2FC > 0 
                                       & DAM2.markers$p_val < 0.05, "Up",
                                         ifelse(DAM2.markers$avg_log2FC < 0 
                                                & DAM2.markers$p_val < 0.05, "Down", "None"))
DAM2.markers$log10p_val <- -log10(DAM2.markers$p_val)
DAM2.markers$Gene <- rownames(DAM2.markers)

DAM2.markers <- DAM2.markers[order(DAM2.markers$avg_log2FC,decreasing = TRUE),]

data <- bind_rows(
  DAM2.markers %>%
    filter(expression == 'Up') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10),
  DAM2.markers %>%
    filter(expression == 'Down') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10)
)

write.csv(DAM2.markers,file = "DAM2_marker_0.csv")

DAM2.degPlot <- ggplot(DAM2.markers, aes(avg_log2FC,log10p_val))+
  geom_point(size = 1,aes(color = expression))+
  #scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  xlab(expression("log"[2]*"FoldChange"))+
  ylab(expression("-log"[10]*"(p value)")) +
  labs(title = "High vs Low",colour = "Significant")+
  # 更换颜色
  scale_color_manual(values = c("Up" = "#FF6C67","Down" = "#519BFF")) +
  # 标题居中
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_label_repel(data = data,
                   aes(avg_log2FC, log10p_val, label = Gene),
                   size = 2,max.overlaps=Inf)

cluster2 <- subset(mouse,seurat_clusters == 2)
Idents(cluster2) <- cluster2@meta.data$Sample

cluster2.markers <- FindMarkers(cluster2, ident.1 = "FAD_4A", ident.2 = "FAD",logfc.threshold = 0, min.pct = 0.25)
cluster2.markers <- cluster2.markers[order(cluster2.markers$avg_log2FC),]

cluster2.markers$expression <- ifelse(cluster2.markers$avg_log2FC > 0.1
                                  & cluster2.markers$p_val < 0.05, "Up",
                                  ifelse(cluster2.markers$avg_log2FC < -0.1 
                                         & cluster2.markers$p_val < 0.05, "Down", "None"))
cluster2.markers$log10p_val <- -log10(cluster2.markers$p_val)
cluster2.markers$Gene <- rownames(cluster2.markers)

cluster2.markers <- cluster2.markers[order(cluster2.markers$avg_log2FC,decreasing = TRUE),]

data <- bind_rows(
  cluster2.markers %>%
    filter(expression == 'Up') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10),
  cluster2.markers %>%
    filter(expression == 'Down') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10)
)

write.csv(cluster2.markers,file = "cluster2_marker_0.csv")

cluster2.degPlot <- ggplot(cluster2.markers, aes(avg_log2FC,log10p_val))+
  geom_point(size = 1,aes(color = expression))+
  #scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  xlab(expression("log"[2]*"FoldChange"))+
  ylab(expression("-log"[10]*"(p value)")) +
  labs(title = "High vs Low",colour = "Significant")+
  # 更换颜色
  scale_color_manual(values = c("Up" = "#FF6C67","Down" = "#519BFF")) +
  # 标题居中
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_label_repel(data = data,
                   aes(avg_log2FC, log10p_val, label = Gene),
                   size = 2,max.overlaps=Inf)

cluster3 <- subset(mouse,seurat_clusters == 3)
Idents(cluster3) <- cluster3@meta.data$Sample

cluster3.markers <- FindMarkers(cluster3, ident.1 = "FAD_4A", ident.2 = "FAD",logfc.threshold = 0, min.pct = 0.25)
cluster3.markers <- cluster3.markers[order(cluster3.markers$avg_log2FC),]

cluster3.markers$expression <- ifelse(cluster3.markers$avg_log2FC > 0 
                                  & cluster3.markers$p_val < 0.05, "Up",
                                  ifelse(cluster3.markers$avg_log2FC < 0 
                                         & cluster3.markers$p_val < 0.05, "Down", "None"))
cluster3.markers$log10p_val <- -log10(cluster3.markers$p_val)
cluster3.markers$Gene <- rownames(cluster3.markers)

cluster3.markers <- cluster3.markers[order(cluster3.markers$avg_log2FC,decreasing = TRUE),]

data <- bind_rows(
  cluster3.markers %>%
    filter(expression == 'Up') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10),
  cluster3.markers %>%
    filter(expression == 'Down') %>%
    arrange(p_val,desc(abs(avg_log2FC))) %>%
    head(10)
)

write.csv(cluster3.markers,file = "cluster3_marker_0.csv")

cluster3.degPlot <- ggplot(cluster3.markers, aes(avg_log2FC,log10p_val))+
  geom_point(size = 1,aes(color = expression))+
  #scale_x_continuous(limits = c(-2.5, 2.5))+
  geom_hline(yintercept=-log10(0.05),linetype=4)+
  geom_vline(xintercept=c(-1,1),linetype=4)+
  xlab(expression("log"[2]*"FoldChange"))+
  ylab(expression("-log"[10]*"(p value)")) +
  labs(title = "High vs Low",colour = "Significant")+
  # 更换颜色
  scale_color_manual(values = c("Up" = "#FF6C67","Down" = "#519BFF")) +
  # 标题居中
  theme(plot.title = element_text(hjust = 0.5,size = 16),
        panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA)) +
  geom_label_repel(data = data,
                   aes(avg_log2FC, log10p_val, label = Gene),
                   size = 2,max.overlaps=Inf)

#############################
# miloR
mouse.fad <- subset(mouse,Sample == "FAD")
mouse.fad_4a <- subset(mouse,Sample == "FAD_4A")

Condition.fad <- sample(x = c("A", "B","C"), size = 12390, replace = TRUE)
mouse.fad@meta.data$Group <- Condition.fad

Condition.fad_4a <- sample(x = c("D", "E","F"), size = 9800, replace = TRUE)
mouse.fad_4a@meta.data$Group <- Condition.fad_4a

mouse.ad <- merge(mouse.fad,mouse.fad_4a)
meta <- mouse.ad@meta.data

mouse.ad <- subset(mouse,Sample != "WT")
mouse.ad@meta.data$Group <- meta$Group

mouse.sce <- as.SingleCellExperiment(mouse.ad)
mouse.milo <- Milo(mouse.sce)

mouse.milo <- buildGraph(mouse.milo, k = 20, d = 30,reduced.dim = "PCA")
mouse.milo <- makeNhoods(mouse.milo, prop = 0.1, k = 20, d=30, refined = TRUE,reduced_dims = "PCA")
plotNhoodSizeHist(mouse.milo)
mouse.milo <- countCells(mouse.milo, meta.data = data.frame(colData(mouse.milo)), samples="Group")
head(nhoodCounts(mouse.milo))

mouse.design <- data.frame(colData(mouse.milo))[,c("Group","Sample")]
mouse.design <- distinct(mouse.design)
rownames(mouse.design) <- mouse.design$Group
mouse.design <- mouse.design[colnames(nhoodCounts(mouse.milo)), , drop=FALSE]

mouse.milo <- calcNhoodDistance(mouse.milo, d=30, reduced.dim = "PCA")

da_results <- testNhoods(mouse.milo,design = ~ Sample, design.df = mouse.design)

da_results %>%
  arrange(- SpatialFDR) %>%
  head() 

mouse.milo <- buildNhoodGraph(mouse.milo)

umap_pl <- plotReducedDim(mouse.milo, dimred = "UMAP", colour_by="seurat_clusters", text_by = "seurat_clusters", 
                          text_size = 10, point_size=0.5) +
  guides(fill="none")+
  theme(legend.position = "none")

nh_graph_pl <- plotNhoodGraphDA(mouse.milo, da_results, layout="UMAP",alpha=0.1) 

umap_pl + nh_graph_pl +
  plot_layout(guides="collect")


save(mouse, file = "mouse.RData")
