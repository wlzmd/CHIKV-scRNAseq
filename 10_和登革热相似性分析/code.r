setwd("~/SHIXI/2.基肯孔雅热/10_和登革热相似性分析")
setwd("~/SHIXI/2.基肯孔雅热/10_和登革热相似性分析/SDp_登革热")
list.files("./")

# 获取所有文件
files <- list.files("./", full.names = TRUE)

# 提取样本前缀（去掉 barcodes/features/matrix 部分）
prefixes <- unique(gsub("_(barcodes|features|matrix).*", "", basename(files)))

# 遍历每个样本前缀
for (pref in prefixes) {
  # 创建文件夹
  dir.create(pref, showWarnings = FALSE)
  
  # 找到该前缀对应的文件
  sample_files <- files[grepl(pref, files)]
  
  # 遍历并复制/移动到文件夹里
  for (f in sample_files) {
    # 去掉前缀，只保留后缀部分
    new_name <- sub(paste0(pref, "_"), "", basename(f))
    
    # 构建目标路径
    target <- file.path(pref, new_name)
    
    # 移动文件（如果想复制用 file.copy）
    file.rename(f, target)
  }
}

library(Seurat)

# 列出所有样本文件夹
dirs <- list.files("./", full.names = TRUE)

# 逐个读取文件夹并创建 Seurat 对象
seurat.list <- list()
for (d in dirs) {
  sample_id <- basename(d)
  data <- Read10X(data.dir = d)
  rownames(data) <- make.unique(rownames(data))
  rn <- gsub("_", "-", rownames(data))
  rownames(data) <- make.unique(rn)
  seurat.obj <- CreateSeuratObject(
    counts = data,
    project = sample_id
    )
  seurat.list[[sample_id]] <- seurat.obj
}

# 合并所有样本为一个 Seurat 对象
combined <- merge(
  seurat.list[[1]], 
  y = seurat.list[-1], 
  add.cell.ids = names(seurat.list), 
  project = "GSE220969"
)
combined_seurat <- readRDS("~/SHIXI/2.基肯孔雅热/0.0rawdata/combined_seurat.rds")
chikv<-subset(combined_seurat,subset = group == "CKGYF")

setwd("~/SHIXI/2.基肯孔雅热/10_和登革热相似性分析")
library(Seurat)
library(harmony)

# Step 1: 合并两个 Seurat 对象
merged <- merge(combined, y = chikv, add.cell.ids = c("combined", "chikv"), project = "JointAnalysis")

# Step 2: 标准预处理
merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)

# Step 3: PCA 降维
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)

# Step 4: Harmony 批次效应校正
# 假设你有一个元数据列，比如 "orig.ident" 或 "group" 表示批次/样本来源
merged <- RunHarmony(merged, group.by.vars = "orig.ident")

# Step 5: 用 Harmony 的结果做下游分析
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.4)
merged <-JoinLayers(merged)
markers_merge<-FindAllMarkers(merged,only.pos = T,logfc.threshold = 0.80,min.pct = 0.1)
ids <- c("NK","T","Mono","T","B","T",
         "Mono","Plasma","Plasma","Mono","Na","Platelets",
         "Mono","Na","Mono","DC","B","Plasma",
         "B","pDC","Na","Na","B","Platelets",
         "NK","Na","Na","Na","Mono"
         )
names(ids) <- levels(merged)
merged <- RenameIdents(merged, ids)
Idents(merged)<-factor(Idents(merged),levels = c("NK","T","Mono","DC","pDC","B","Plasma",
                 "Platelets","Na"))

merged_sub<-subset(merged,idents="Na",invert=T)
# 给 merged_sub 添加分组信息
merged_sub$group <- ifelse(grepl("^GSM", merged_sub$orig.ident), 
                           "DENV", 
                           "CKGYF")
DimPlot(merged_sub, reduction = "umap", label = F)+NoLegend()+NoAxes()
dev.off()
merged_sub$celltype <- Idents(merged_sub)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

genes <- c(
  "TNFRSF12A","SEMA6B","NAB2","ASPH","EREG","FCAR","CHKA","TUBB6","THBD",
  "P2RX1","USP18","DDX58","NHSL2","TRMT1L","SERPING1","RSAD2","NEXN",
  "AIM2","DEFB1","CXCL10","MX1","ISG15"
)

DotPlot(merged_sub, features = genes, group.by = "group") +
  coord_flip()+
  scale_color_gradient(low = "grey80", high = "orange") +   # 点的颜色渐变
  scale_size(range = c(2, 8)) +                             # 点大小范围
  theme_minimal(base_size = 12) +                           # 简洁主题
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5,colour = "black",size = 12),      # 横坐标组名倾斜
    axis.text.y = element_text(size = 12,colour = "black"),   # 纵坐标基因加粗
    panel.grid = element_blank(),                           # 去掉网格线
    axis.title.x = element_text(size = 0),
    axis.title.y = element_text(size = 0),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  ) +
  labs(
    x = "Celltype Group", 
    y = "Gene", 
    color = "Avg Expr", 
    size = "Pct Expr"
  )

Idents(combined_seurat)<-combined_seurat$group
Idents(combined_seurat)<-factor(combined_seurat$group,levels = c("Con","CKGYF"))
VlnPlot(combined_seurat,features = c("MX1","ISG15"),raster=FALSE,pt.size = 0)




merged_sub <- AddModuleScore(merged_sub, features=list(genes), name="IFNs")

df_sample_ct <- merged_sub@meta.data %>%
  dplyr::select(orig.ident, group, celltype, IFNs1) %>%
  group_by(orig.ident, group, celltype) %>%
  summarise(IFNs = mean(IFNs1),
            .groups="drop") %>%
  pivot_longer(cols=c(IFNs),
               names_to="ScoreType",
               values_to="Score")

df_sample_ct$Condition <- factor(df_sample_ct$group, levels=c("CKGYF","DENV"))
unique(df_sample_ct$celltype)
df_sample_ct$seurat_clusters <- factor(df_sample_ct$celltype, levels = 
                                         c("NK","T","Mono","DC","pDC","B","Plasma",
                                           "Platelets","Na"))


ggplot(df_sample_ct, aes(x=seurat_clusters, y=Score, fill=Condition)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8))+
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.2, alpha = 0.8) +
  facet_wrap(~ScoreType, ncol=1, scales="free_y") +
  theme_classic() +
  labs(x=NULL, y=" ")+
  labs(title = NULL) + # 去掉标题 
  theme( axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12), # 横坐标文字旋转90°，黑色，字体放大 
         axis.text.y = element_text(color = "black", size = 12), # 纵坐标文字黑色，字体放大
         strip.text = element_blank(),
         panel.spacing = unit(1.3, "cm") # 调整分面之间的间距，数值越大间距越大
  )

ggplot(df_sample_ct, aes(x=seurat_clusters, y=Score, fill=Condition)) +
  geom_boxplot(width = 0.7, outlier.shape = NA,
               position = position_dodge(width = 0.8)) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 0.2, alpha = 0.8
  ) +
  facet_wrap(~ScoreType, ncol=1, scales="free_y") +
  

  stat_compare_means(
    aes(group = Condition),
    method = "t.test",   # 或 "t.test"
    label = "p.signif",       # 显示星号
    hide.ns = T ,          # 不显示 ns
    size = 5  
  ) +
  
  theme_classic() +
  labs(x=NULL, y=" ") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    strip.text = element_blank(),
    panel.spacing = unit(1.3, "cm")
  )





Mono <- subset(merged_sub, subset = celltype == "Mono")
Idents(Mono) <- Mono$group
table(Idents(Mono))

de <- FindMarkers(
  object = Mono,
  ident.1 = "CKGYF",
  ident.2 = "DENV",
  test.use = "wilcox",
  logfc.threshold = 0,
  min.pct = 0.4,
  only.pos = FALSE
)

de <- de %>%
  tibble::rownames_to_column("gene") %>%
  mutate(p_val_adj = ifelse(is.na(p_val_adj), 1, p_val_adj)) %>%
  arrange(p_val_adj, desc(avg_log2FC))

# 过滤显著差异基因（阈值你可改）
de_sig <- de %>% filter(p_val_adj < 0.05)
de_sig <- subset(de_sig, pct.1 > 0.05 & pct.2 > 0.05)




## =========================
## 3) Top10 up/down genes & Heatmap (x = samples orig.idents)
## =========================
top_up <- de_sig %>% arrange(desc(avg_log2FC)) %>% slice_head(n = 10) %>% pull(gene)
top_dn <- de_sig %>% arrange(avg_log2FC)       %>% slice_head(n = 10) %>% pull(gene)
top_genes <- unique(c(top_up, top_dn))

# 以“样本”为单位聚合表达：列=样本（orig.idents）
# 推荐用 slot="data"（log-normalized）
avg_expr_list <- AverageExpression(
  Mono,
  assays = DefaultAssay(Mono),
  features = top_genes,
  group.by = "orig.ident",
  slot = "data",
  verbose = FALSE
)

avg_mat <- avg_expr_list[[DefaultAssay(Mono)]]  # genes x samples

# 行标准化（更适合热图展示）
avg_mat_z <- t(scale(t(as.matrix(avg_mat))))
avg_mat_z[is.na(avg_mat_z)] <- 0

# 样本注释（每个样本对应的 group；若一个样本里有多组会报错/不合理）
sample_anno <- Mono@meta.data %>%
  distinct(orig.ident, group) %>%
  group_by(orig.ident) 
# 对齐列名
sample_anno <- sample_anno[colnames(avg_mat_z), , drop = FALSE]

# 热图：横坐标=样本名（orig.idents）
pheatmap(
  mat = avg_mat_z,
  annotation_col = sample_anno,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = TRUE,
  fontsize_col = 12,
  fontsize_row = 12
  )



## 如果你还想按“上/下调”分块显示（可选）
gene_dir <- data.frame(
  direction = ifelse(rownames(avg_mat_z) %in% top_up, "Up", "Down"),
  row.names = rownames(avg_mat_z)
)

pheatmap(
  mat = avg_mat_z,
  annotation_col = sample_anno,
  annotation_row = gene_dir,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  show_colnames = TRUE,
  fontsize_col = 9,
  fontsize_row = 9,
  main = "Mono Top DE genes (Up/Down)"
)




#Mono
geneupCvsD<-dplyr::filter(de_sig,avg_log2FC>0.58)$gene
de_sig_CvsCon<-FindMarkers(subset(combined_seurat,ident = "Mono"),group.by = "group",
                           ident.1 = "CKGYF",
                           ident.2 = "Con",test.use = "wilcox",
                           logfc.threshold = 0.58,
                           min.pct = 0.4,
                           only.pos = FALSE)

geneset_featuresinCKGYF<-intersect(geneupCvsD,rownames(dplyr::filter(de_sig_CvsCon,avg_log2FC>0.58)))

library(clusterProfiler)
library(org.Hs.eg.db)
ego<-enrichGO(
  gene = geneset_featuresinCKGYF,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head(ego)
ego<-data.frame(ego@result)
ego_sub<-ego[c("GO:0048002",
"GO:0051607",
"GO:0034341",
"GO:0071357","GO:0045071"),]


ego_sub <- ego_sub %>%
  arrange(p.adjust) %>%   # 按 p 值从小到大
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

library(ggplot2)
ego_sub$logP

ggplot(ego_sub,
       aes(x = Count,
           y = Description,
           fill = logP)) +
  geom_col(width = 0.7) +

  scale_fill_gradient(
    low = "#2166AC",   # 蓝色（不显著）
    high = "#B2182B"   # 红色（显著）
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black",size = 12))+
  labs(
    x = "Gene number",
    y = NULL,
    fill = "-Log10(adj.p)"
  )

VlnPlot(Mono,features = c("HLA-DQA1","IFI30","UNC93B1"),pt.size = 0)&
  theme(axis.title.x = element_blank())



cytokine_genes <- c(
  "IL2", "IL7", "CSF3", "CXCL10", "CCL2", "IL15", "TNF", "IL6", "CXCL7", "IL1RN",
  "CSF1", "IFNB", "IL20RA", "IL10", "IL18", "HGF", "CXCL9", "CCL27", "TGFB1", "IL1B",
  "LTA", "CSF2", "IL12", "TNFSF13", "IL11", "CCL12", "CXCL8", "CXCL11", "CCL4", "CXCL1",
  "CXCL2", "CCL3", "CCL18", "CXCL16", "FGA", "CCL5", "IL1A", "CCL20", "OSM", "TNFSF14",
  "FGF10", "CXCL5", "IL9", "IL1C", "TNFSF10", "FGF7", "IL23", "TGFB3", "TNFRSF1B", "CREB1",
  "TXN3A", "IL1W", "XCL2", "XCL1", "CCL25", "CCL16", "CD40LG", "IL20", "FASLG", "F10",
  "SCYL1", "GDF15", "TNFSF8", "GDF11", "VEGFA", "GDF7", "PDGFA", "IL21", "ADGRL1", "ADGRD2",
  "PDGFB", "TNFAIP6", "TAMPAI", "BMP6", "PDGFD", "IL20RB", "GHR", "WISP3", "MIP", "IL4",
  "BMP8B", "TNFSF12", "SCYL2", "TSP1", "SDF1B", "NMBR", "PRP", "FGF11", "PGR-AS1", "FGF2",
  "VEGFC", "IL12A", "BAFF", "AMH", "IL26", "PDGFC", "FGF23", "KLF11", "TNFSF9", "BMP7",
  "IL24", "BMP10", "CXCL6", "CD96", "IL16", "CTD", "FGF9", "PSNL1", "TSC1", "IL32", "IL10R",
  "SPPL", "IL12RB1", "BMP4", "TARBP1", "TGFBR2", "A20MT", "EDA", "TNFRSF11B", "CCL28"
)

inflam_genes   <- c(
  "ABCA1", "ABI3", "ACKR1", "ACKR3", "ADORA2A", "ADORA2B", "AHR", "AIM2", "AKIRIN2", "AP2A2",
  "ATF3", "ATP6V1H", "BCL2A1", "BIRC3", "BST1", "BST2", "BTK", "CARD11", "CASP4", "CCL17",
  "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CCL7", "CCR1", "CCR2", "CCR5", "CD14",
  "CD163", "CD19", "CD3D", "CD3G", "CD4", "CD74", "CENPK", "CHST2", "CHST7", "CHURC1",
  "CIITA", "CSF1", "CSF2", "CSF3", "CXCL10", "CXCL11", "CXCL6", "CXCL9", "CYBB", "CYLD",
  "DCLRE1B", "DUSP1", "EGR1", "EPS8L2", "ERO1A", "F3", "FADD", "FAS", "FBP1", "FGR",
  "FTH1", "G0S2", "GABARAPL1", "GAS1", "GAS6", "GATA3", "GBP1", "GBP1P1", "GBP2", "GBP3",
  "GBP4", "GBP5", "GBP6", "GBP7", "GPR183", "GPR84", "HCK", "HIF1A", "HLA-DRA", "ICAM1",
  "ICOS", "IFI16", "IFI30", "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1",
  "IFITM2", "IFITM3", "IGFBP3", "IL10", "IL10RA", "IL10RB", "IL15", "IL15RA", "IL18", "IL1R1",
  "IL1R2", "IL1RAP", "IL1RN", "IL2RA", "IL4R", "IL6", "IL6R", "IL7R", "IL8", "IRAK3",
  "IRF1", "IRF7", "ITGAX", "ITGB2", "ITGB7", "KCNA3", "KCNK10", "KCNK3", "KLF6", "LAMP3",
  "LAPTM5", "LCN2", "LCP2", "LTA", "LY86", "LYN", "MAP3K8", "MAF", "MAFB", "MARCO",
  "MEFV", "MEP1A", "MEP1B", "MERTK", "MMP14", "MMP9", "MNDA", "MRC1", "MSR1", "NAB2",
  "NCR1", "NFKB1", "NFKB2", "NLRP1", "NLRP3", "NOS2", "NQO1", "NQO2", "OAS1", "OAS2",
  "OAS3", "ORM1", "ORM2", "P2RX7", "P2RY6", "P2RY13", "P2RY14", "PADI4", "PARK7", "PBX1",
  "PDCD1", "PDCD1LG2", "PDE4B", "PDE4D", "PGLYRP1", "PIK3R5", "PLAU", "PLSCR1", "PSEN1", "PTAFR",
  "PTGER2", "PTGER4", "PTGS2", "PTPN6", "PYCARD", "RIPK1", "RIPK3", "RNASE2", "RPL13A", "S100A8",
  "S100A9", "S100B", "SAMD9L", "SASH3", "SBNO2", "SELL", "SEMA4D", "SERPINB1", "SIGLEC1", "SIGLEC14",
  "SIGLEC5", "SIGLEC7", "SIGLEC9", "SLC11A1", "SLC15A3", "SLC2A3", "SLC39A8", "SLC7A2", "SMAD3", "SOCS1",
  "SOCS3", "SPI1", "SRGN", "STAT1", "TACR1", "TACR3", "TAPBP", "TBX21", "TBX2", "TLR1",
  "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TNFAIP3", "TNFRSF1B", "TNFRSF9",
  "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNF", "TYROBP"
)

# 保存为只有一列的 CSV
write.csv(
  data.frame(gene = inflam_genes),   # 转成一列的数据框
  file = "inflam_genes.csv",         # 输出文件名
  row.names = FALSE,                 # 不要行名
  quote = FALSE                      # 不加引号（可选）
)


Mono <- AddModuleScore(
  Mono,
  features = list(cytokine_genes),
  name = "CytokineScore"
)
Mono <- AddModuleScore(
  Mono,
  features = list(inflam_genes),
  name = "InflamScore"
)


head(Mono@meta.data[, c("CytokineScore1", "InflamScore1", "group")])

VlnPlot(Mono,"InflamScore1",pt.size = 0)



Idents(combined_seurat)<-combined_seurat$seurat_clusters
ids <- c("NK","Others","Others","Mono_CD14","Others",
         "Mono_CD14","NK","Others","Mono_CD16","Others",
         "Others","Mono_CD16","Others","Others","Others",
         "Others","Others","Others","Mono_CD14","Others",
         "Others","DC_CD1C","Mono_NEAT1","pDC","Others"
)

names(ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, ids)

combined_seurat <- AddModuleScore(
  combined_seurat,
  features = list(inflam_genes),
  name = "InflamScore"
)



Idents(combined_seurat)<-factor(Idents(combined_seurat),
                       levels=c("Others","NK","Mono_CD14","Mono_CD16",
                                "Mono_NEAT1","DC_CD1C","pDC"
                       ))

combined_seurat$group <- factor(combined_seurat$group,
                                levels=c("Con","CKGYF"
                                ))





p <- VlnPlot(
  combined_seurat,
  features = "InflamScore1",
  pt.size = 0,
  split.by = "group",cols = c("turquoise","orange")
) + NoLegend()

# 2) 取出作图同源数据：分数 + 细胞类型(Identity) + group
df <- FetchData(combined_seurat, vars = c("InflamScore1", "group"))
df$Identity <- Idents(combined_seurat)

# 可选：去掉缺失
df <- df %>% filter(!is.na(group), !is.na(InflamScore1))

# 3) 把检验结果加到原图上：
# group.by = "Identity" 代表“按细胞类型分层，在每个层里比较group”
p + stat_compare_means(
  data = df,
  aes(x = Identity, y = InflamScore1, group = group),
  method = "wilcox.test",
  label = "p.signif",
  group.by = "Identity"
)



combined_seurat_monocon <- subset(
  combined_seurat,
  subset = (Idents(combined_seurat) %in% c("Mono_CD14", "Mono_CD16", "Mono_NEAT1")) & group == "Con"
)

merged_Mono <- merge(Mono, y = combined_seurat_monocon,
                     project = "JointAnalysis")

table(merged_Mono$group)
# Step 2: 标准预处理
merged_Mono <- NormalizeData(merged_Mono)
merged_Mono <- FindVariableFeatures(merged_Mono, selection.method = "vst", nfeatures = 2000)
merged_Mono <- ScaleData(merged_Mono)


merged_Mono <- AddModuleScore(
  merged_Mono,
  features = list(cytokine_genes),
  name = "CytokineScore"
)
merged_Mono <- AddModuleScore(
  merged_Mono,
  features = list(inflam_genes),
  name = "InflamScore"
)



merged_Mono$group <- factor(merged_Mono$group,levels = c("Con","CKGYF","DENV"))

my_comp <- list(
  c("Con", "CKGYF"),
  c("Con", "DENV"),
  c("CKGYF", "DENV")
)

VlnPlot(
  merged_Mono,
  features = "InflamScore1",
  pt.size = 0,
  group.by = "group"
) +NoLegend()

