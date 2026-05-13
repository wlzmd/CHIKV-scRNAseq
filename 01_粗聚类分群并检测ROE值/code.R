# 读取 md5.txt
ref <- read.table("md5.txt", stringsAsFactors = FALSE)
colnames(ref) <- c("md5", "file")
ref$file<-basename(ref$file)

# 列出需要检测的文件（排除 md5.txt 本身）
files <- list.files("./Matrix", recursive = TRUE, full.names = TRUE)
# 计算 MD5
calc <- tools::md5sum(files)
# 转换为数据框
calc_df <- data.frame(file = basename(names(calc)), md5 = as.character(calc), stringsAsFactors = FALSE)
# 合并参考和计算结果
compare <- merge(ref, calc_df, by = "file", all.x = TRUE)

# 增加是否一致的标记
compare$match <- compare$md5.x == compare$md5.y

# 查看比对结果
print(compare)



##读取Seurat对象
setwd("~/SHIXI/3.基肯孔雅热/01_聚类分群并检测ROE值")
library(Seurat)
library(Matrix)
library(ggplot2)

## 已有的 files 向量
files <- sort(list.files("./Matrix", pattern = "\\.tar\\.gz$", full.names = TRUE))

# 辅助函数：读一个 tar.gz，返回一个 Seurat 对象 ---------------------------
read_one_sample <- function(tar_path) {
  # 临时目录
  tmpdir <- tempfile(pattern = "one_sample_")
  dir.create(tmpdir)
  
  # 无论成功还是报错，函数退出时都删除临时目录，释放磁盘空间
  on.exit({
    if (dir.exists(tmpdir)) {
      unlink(tmpdir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  
  # 解压
  untar(tar_path, exdir = tmpdir)
  
  # 样本名 = tar.gz 文件名去掉后缀
  sample_id <- basename(tar_path)
  sample_id <- sub("\\.tar\\.gz$", "", sample_id)
  
  # 10x 矩阵目录：<tmpdir>/<sample_id>/02.count/filter_matrix
  mat_dir <- file.path(tmpdir, sample_id, "02.count", "filter_matrix")
  if (!dir.exists(mat_dir)) {
    stop("没有找到目录: ", mat_dir)
  }
  
  # 读 10x 格式（你指定 gene.column = 1）
  mat <- Read10X(data.dir = mat_dir, gene.column = 1)
  
  # 创建 Seurat 对象
  seu <- CreateSeuratObject(
    counts       = mat,
    project      = sample_id,
    min.cells    = 0,
    min.features = 0
  )
  
  # 给细胞名加上样本前缀，避免不同样本细胞名冲突
  colnames(seu) <- paste(sample_id, colnames(seu), sep = "_")
  
  return(seu)
}

# 批量读入每个样本，生成 seurat list ----------------------------------------
seu_list <- list()
seu_list <- lapply(files, read_one_sample)
names(seu_list) <- sub("\\.tar\\.gz$", "", basename(files))

# 合并为一个大的 Seurat 对象 -----------------------------------------------
combined_seurat <- merge(
    x = seu_list[[1]],
    y = seu_list[-1]
  )

# 最终结果
combined_seurat
combined_seurat <- NormalizeData(combined_seurat) 
combined_seurat[["percent.mt"]] <- PercentageFeatureSet(combined_seurat, pattern = "^MT-")
p<-VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave("combined_seurat_plot.pdf",p)

plot1 <- FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(combined_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

combined_seurat <- subset(combined_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)
combined_seurat<-FindVariableFeatures(combined_seurat) %>%
  ScaleData() %>% 
  RunPCA() 

ElbowPlot(combined_seurat,ndims = 50)

combined_seurat <- IntegrateLayers(object = combined_seurat, 
                                method = CCAIntegration, 
                                orig.reduction = "pca",
                                new.reduction = "integrated.cca",
                                verbose = FALSE)

combined_seurat[["RNA"]] <- JoinLayers(combined_seurat[["RNA"]])


combined_seurat <- FindNeighbors(combined_seurat, reduction = "integrated.cca", dims = 1:30)
combined_seurat <- FindClusters(combined_seurat, resolution = 0.3)
combined_seurat <- RunUMAP(combined_seurat, dims = 1:30, reduction = "integrated.cca")
DimPlot(combined_seurat, reduction = "umap", group.by = c("RNA_snn_res.0.3"),label = T)
markerofallsamples<-FindAllMarkers(combined_seurat,logfc.threshold = 0.58,only.pos = T)
FeaturePlot(combined_seurat, features = "CD8A")
FeaturePlot(combined_seurat, features = "CD3E")
ids <- c("NK","CD8","CD4","Mono","CD4","Mono","NK","B","Mono",
         "CD8","B","Mono","γδ T","Erythroid cells","Plasma","CD8",
         "Platlet","CD4","Mono","CD4","B","DC","Mono",
         "DC","CD8")
names(ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, ids)
Idents(combined_seurat)<-factor(Idents(combined_seurat),levels = c("B","Plasma",
    "NK","CD8","CD4","γδ T","Mono","DC","Erythroid cells","Platlet"))

p <- DimPlot(combined_seurat, label = TRUE, repel = TRUE) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks   = element_blank(),
    axis.line = element_blank(),
    plot.title   = element_text(size = 8, hjust = 0.5),
    legend.text  = element_text(size = 8),
    legend.title = element_text(size = 8),
  )

# 保存为 PDF，分辨率 600 dpi
ggsave("DimPlot.png", plot = p, device = "png", dpi = 600, width = 5, height = 4)


DotPlot(
  combined_seurat,
  features = c("CD79A","CD79B","MS4A1","MZB1","NKG7","KLRF1",
               "CD3D","CD3E","CD3G","CD8A","CD8B","CD40LG",
               "TRDC","TRGC1","CST3","LYZ","CD14","CD1C","ITGAX",
               "HBB","HBA1","HBA2","PPBP","GP9")
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 8),
        axis.text.y = element_text(size=8),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  scale_color_gradient(low = "white", high = "red")

table(combined_seurat$orig.ident)
combined_seurat$group <- gsub("[0-9]+$", "", combined_seurat$orig.ident)

# 构建排序后的水平
ckgyf_levels <- paste0("CKGYF", sort(seq(1,29,1)))
con_levels   <- paste0("Con",   sort(seq(1,10,1)))

# 合并成最终因子水平：CKGYF在前，Con在后
final_levels <- c(ckgyf_levels, con_levels)

# 重新设置因子顺序
combined_seurat$ids <- factor(combined_seurat$orig.ident, levels = final_levels)

cell_counts <- table(combined_seurat$ids)
df <- as.data.frame(cell_counts)
colnames(df) <- c("sample", "n_cells")

sick_samples <- paste0("CKGYF", sort(seq(1,29,1)))
heal_samples <- paste0("Con",   sort(seq(1,10,1)))

df$condition <- NA_character_
df$condition[df$sample %in% sick_samples] <- "CKGYF"
df$condition[df$sample %in% heal_samples] <- "Con"

## 检查有没有没标到组的样本
df[is.na(df$condition), ]   # 这里如果有行，就说明有样本没分组


df <- df %>%
  mutate(log10_n = log10(n_cells),
         condition = factor(condition, levels = c("Con","CKGYF"))) %>%
  group_by(condition) %>%
  arrange(desc(n_cells), .by_group = TRUE) %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  ungroup()

## 4. 画图：横向条形图 + 分面 + 颜色 --------------------------------
ggplot(df, aes(x = sample, y = log10_n, fill = condition)) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_grid(condition ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(CKGYF = "#E41A1C",   # 红
                               Con = "#377EB8"   # 蓝
                               )) + 
  labs(x = NULL, y = expression(Log[10]("cell number"))) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = 12),
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10)
  )

# 保存 Seurat 对象
saveRDS(combined_seurat, file = "./combined_seurat.rds")


Idents(combined_seurat)<-combined_seurat$seurat_clusters
ids <- c("NK","CD8⁺T_effetcor","CD4⁺T_naive/CM","CD14⁺Mono","CD4⁺T_naive/CM","CD14⁺Mono",
         "NK","B","CD16⁺Mono",
         "CD8⁺T_naive/CM","B","CD16⁺Mono","γδ T","Erythroid cells","Plasma","Prop T",
         "Platlet","CD4⁺T_naive/CM","CD14⁺Mono","CD4⁺Treg","B","CD1C⁺DC","LRMDA⁺Mono",
         "LILRA4⁺DC","CD8⁺T_effetcor")
names(ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, ids)
factor(ids)

##检测R O/E值和单因素方差分析
setwd("~/SHIXI/3.基肯孔雅热/01_粗聚类分群并检测ROE值")
seu_re <- readRDS("~/SHIXI/3.基肯孔雅热/0.0rawdata/combined_seurat.rds")
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)

seu_re$celltype<-Idents(seu_re)
meta <- seu_re@meta.data %>%
  dplyr::select(orig.ident, celltype, group) %>%
  mutate(
    orig.ident = as.character(orig.ident),
    celltype   = as.character(celltype),
    group      = as.character(group)
  )

meta$group    <- factor(meta$group, levels = c("Con","CKGYF"))

global_prop <- meta %>%
  dplyr::count(celltype, name = "global_n") %>%
  mutate(global_p = global_n / sum(global_n))

# 每组 observed
obs <- meta %>%
  dplyr::count(group, celltype, name = "obs_n")

# 每组总细胞数
group_total <- meta %>%
  dplyr::count(group, name = "group_n")

# 合并计算 expected 和 R_OE
roe_df <- obs %>%
  left_join(group_total, by = "group") %>%
  left_join(global_prop, by = "celltype") %>%
  mutate(
    exp_n = group_n * global_p,
    R_OE  = obs_n / exp_n,
    log2R = log2(R_OE),
    status = ifelse(R_OE >= 1, "Enrichment", "Depletion")
  )

# 画气泡图（点大小=R_OE，红/蓝=富集/缺失）
p_roe <- ggplot(roe_df, aes(x = celltype, y = group)) +
  geom_point(aes(size = R_OE, color = status), alpha = 0.95) +
  scale_color_manual(values = c(Depletion = "#2c7fb8", Enrichment = "#d7301f")) +
  scale_size(range = c(2, 6),
             breaks = c(0.8, 1, 1.2)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL, size = "R O/E", color = NULL)

p_roe


# heatmap
# 每个样本的 celltype 比例
prop_sample <- meta %>%
  dplyr::count(orig.ident, group, celltype, name = "n") %>%
  group_by(orig.ident) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

anova_res <- prop_sample %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(summary(aov(prop ~ group))[[1]][["Pr(>F)"]][1], error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(q_value = p.adjust(p_value, method = "BH")) %>%
  arrange(q_value)

library(dplyr)

stat_df <- prop_sample %>%
  group_by(celltype) %>%
  summarise(
    p_value = tryCatch(
      summary(aov(prop ~ group))[[1]][["Pr(>F)"]][1],
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    label = ifelse(!is.na(p_value) & p_value < 0.05, "*", "")
  )

plot_df <- prop_sample %>%
  group_by(celltype, group) %>%
  summarise(
    mean_prop = mean(prop),
    sem_prop  = sd(prop) / sqrt(n()),
    .groups = "drop"
  )
plot_df <- plot_df %>%
  left_join(stat_df, by = "celltype")
library(ggplot2)

ggplot(plot_df, aes(x = celltype, y = mean_prop, fill = group)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  geom_errorbar(
    aes(
      ymin = mean_prop - sem_prop,
      ymax = mean_prop + sem_prop
    ),
    width = 0.2,
    position = position_dodge(width = 0.8)
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    x = NULL,
    y = "Cell proportion"
  )

