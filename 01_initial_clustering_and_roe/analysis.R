# Read the reference checksums.
ref <- read.table(here::here("data", "Matrix", "md5.txt"), stringsAsFactors = FALSE)
colnames(ref) <- c("md5", "file")
ref$file<-basename(ref$file)

# List files to verify, excluding md5.txt itself.
files <- list.files(here::here("data", "Matrix"), recursive = TRUE, full.names = TRUE)
files <- files[basename(files) != "md5.txt"]
# Calculate MD5 checksums.
calc <- tools::md5sum(files)
# Convert the results to a data frame.
calc_df <- data.frame(file = basename(names(calc)), md5 = as.character(calc), stringsAsFactors = FALSE)
# Join reference and calculated checksums.
compare <- merge(ref, calc_df, by = "file", all.x = TRUE)

# Flag matching checksums.
compare$match <- compare$md5.x == compare$md5.y

# Inspect the comparison.
print(compare)



## Read and construct the Seurat object.
setwd(here::here("01_initial_clustering_and_roe"))
library(Seurat)
library(Matrix)
library(ggplot2)

## Input archive paths.
files <- sort(list.files(here::here("data", "Matrix"), pattern = "\\.tar\\.gz$", full.names = TRUE))

# Read one tar.gz archive and return a Seurat object.
read_one_sample <- function(tar_path) {
  # Create an isolated temporary directory.
  tmpdir <- tempfile(pattern = "one_sample_")
  dir.create(tmpdir)
  
  # Remove temporary data when the function exits.
  on.exit({
    if (dir.exists(tmpdir)) {
      unlink(tmpdir, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)
  
  # Extract the archive.
  untar(tar_path, exdir = tmpdir)
  
  # Derive the sample identifier from the archive name.
  sample_id <- basename(tar_path)
  sample_id <- sub("\\.tar\\.gz$", "", sample_id)
  
  # Expected 10x matrix directory: <tmpdir>/<sample_id>/02.count/filter_matrix.
  mat_dir <- file.path(tmpdir, sample_id, "02.count", "filter_matrix")
  if (!dir.exists(mat_dir)) {
    stop("Matrix directory not found: ", mat_dir)
  }
  
  # Read the 10x matrix with gene symbols in column 1.
  mat <- Read10X(data.dir = mat_dir, gene.column = 1)
  
  # Create the Seurat object.
  seu <- CreateSeuratObject(
    counts       = mat,
    project      = sample_id,
    min.cells    = 0,
    min.features = 0
  )
  
  # Prefix cell barcodes with the sample identifier.
  colnames(seu) <- paste(sample_id, colnames(seu), sep = "_")
  
  return(seu)
}

# Read all samples into a list of Seurat objects.
seu_list <- list()
seu_list <- lapply(files, read_one_sample)
names(seu_list) <- sub("\\.tar\\.gz$", "", basename(files))

# Merge all samples into one Seurat object.
combined_seurat <- merge(
    x = seu_list[[1]],
    y = seu_list[-1]
  )

# Inspect the merged object.
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

# Save the plot at publication resolution.
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

# Construct ordered factor levels.
ckgyf_levels <- paste0("CKGYF", sort(seq(1,29,1)))
con_levels   <- paste0("Con",   sort(seq(1,10,1)))

# Place CHIKV samples before controls.
final_levels <- c(ckgyf_levels, con_levels)

# Apply the factor order.
combined_seurat$ids <- factor(combined_seurat$orig.ident, levels = final_levels)

cell_counts <- table(combined_seurat$ids)
df <- as.data.frame(cell_counts)
colnames(df) <- c("sample", "n_cells")

sick_samples <- paste0("CKGYF", sort(seq(1,29,1)))
heal_samples <- paste0("Con",   sort(seq(1,10,1)))

df$condition <- NA_character_
df$condition[df$sample %in% sick_samples] <- "CKGYF"
df$condition[df$sample %in% heal_samples] <- "Con"

## Check for samples without an assigned condition.
df[is.na(df$condition), ]


df <- df %>%
  mutate(log10_n = log10(n_cells),
         condition = factor(condition, levels = c("Con","CKGYF"))) %>%
  group_by(condition) %>%
  arrange(desc(n_cells), .by_group = TRUE) %>%
  mutate(sample = factor(sample, levels = sample)) %>%
  ungroup()

## Plot horizontal bars, facets, and condition colours.
ggplot(df, aes(x = sample, y = log10_n, fill = condition)) +
  geom_col(width = 0.8) +
  coord_flip() +
  facet_grid(condition ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(CKGYF = "#E41A1C",
                               Con = "#377EB8"
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

# Save the Seurat object.
saveRDS(combined_seurat, file = here::here("data", "combined_seurat.rds"))


Idents(combined_seurat)<-combined_seurat$seurat_clusters
ids <- c("NK","CD8⁺T_effetcor","CD4⁺T_naive/CM","CD14⁺Mono","CD4⁺T_naive/CM","CD14⁺Mono",
         "NK","B","CD16⁺Mono",
         "CD8⁺T_naive/CM","B","CD16⁺Mono","γδ T","Erythroid cells","Plasma","Prop T",
         "Platlet","CD4⁺T_naive/CM","CD14⁺Mono","CD4⁺Treg","B","CD1C⁺DC","LRMDA⁺Mono",
         "LILRA4⁺DC","CD8⁺T_effetcor")
names(ids) <- levels(combined_seurat)
combined_seurat <- RenameIdents(combined_seurat, ids)
factor(ids)

## Calculate observed/expected ratios and sample-level abundance summaries.
setwd(here::here("01_initial_clustering_and_roe"))
seu_re <- readRDS(here::here("data", "combined_seurat.rds"))
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

# Observed cell counts by group.
obs <- meta %>%
  dplyr::count(group, celltype, name = "obs_n")

# Total cell counts by group.
group_total <- meta %>%
  dplyr::count(group, name = "group_n")

# Calculate expected counts and observed/expected ratios.
roe_df <- obs %>%
  left_join(group_total, by = "group") %>%
  left_join(global_prop, by = "celltype") %>%
  mutate(
    exp_n = group_n * global_p,
    R_OE  = obs_n / exp_n,
    log2R = log2(R_OE),
    status = ifelse(R_OE >= 1, "Enrichment", "Depletion")
  )

# Plot observed/expected ratios; point size is R_OE and colour denotes enrichment.
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
# Calculate cell-type proportions for each sample.
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

