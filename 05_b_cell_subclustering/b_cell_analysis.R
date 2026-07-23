setwd(here::here("05_b_cell_subclustering"))
seu_re <- readRDS(here::here("data", "combined_seurat.rds"))

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(readr)
library(tidyr)

Idents(seu_re)<-seu_re$seurat_clusters
Idents(seu_re)
ids <- c("NK","CD8","CD4","Mono","CD4","Mono","NK","B","Mono",
         "CD8","B","Mono","γδ T","Erythroid cells","Plasma","CD8",
         "Platlet","CD4","Mono","CD4","B","DC","Mono",
         "DC","CD8")
names(ids) <- levels(seu_re)
seu_re <- RenameIdents(seu_re, ids)
B_seu <- subset(seu_re ,idents = c("B","Plasma")) %>%
  NormalizeData() %>%
  FindVariableFeatures()%>%
  ScaleData(vars.to.regress = "percent.mt")
B_seu <- RunPCA(B_seu)
B_seu <- RunHarmony(
  object = B_seu,
  group.by.vars = "orig.ident"
)
B_seu <- B_seu %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 0.5) %>%
  RunUMAP(reduction = "harmony", dims = 1:15)

ids <- c("B_Naive","B_Memory","B_Plasma","B_Naive","B_Memory",
         "B_GC","B_iMemory","B_Memory","B_Naive",
         "B_iMemory","B_Dividing_Plasma","B_Naive",
         "B_Plasma","B_Plasma","B_Naive","B_Naive")
names(ids) <- levels(B_seu)
B_seu <- RenameIdents(B_seu, ids)

Idents(B_seu)<-factor(Idents(B_seu),levels = c("B_Naive","B_GC","B_Memory",
                                               "B_iMemory","B_Dividing_Plasma","B_Plasma"))     
DimPlot(B_seu ,label = T,reduction = "umap",label.size = 0)
saveRDS(B_seu,"B_seu.RDS")


# Marker order used in the figure.
genes <- c("CD79A","MS4A1","MZB1","CD38",
           "IGHD","CD27","NEIL1","MKI67")
p <- FeaturePlot(
  object = B_seu,
  features = genes,
  reduction = "umap",
  ncol = 4,
  pt.size = 0.1,
)
ggsave("FeaturePlot_Bcell_markers.png", p, width = 13, height = 7, dpi = 300)


B_seu$seurat_clusters<-factor(Idents(B_seu),levels = c("B_Naive","B_GC","B_Memory",
     "B_iMemory","B_Dividing_Plasma","B_Plasma")) 

B_seu$group<-factor(B_seu$group,levels = c("Con","CKGYF")) 
p<-ggplot(B_seu@meta.data, aes(x=group, fill= seurat_clusters)) + 
  geom_bar(position = "fill") + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=11),
        axis.title.y = element_text(size=11),
        axis.title.x = element_blank(),
        legend.title = element_blank())  +
  labs(y="Percentage of clusters")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
  )
ggsave("Percentage.png",p,dpi=600,width = 3,height =3.5)


p<-ggplot(B_seu@meta.data, aes(x=ids, fill= seurat_clusters)) + 
  geom_bar(position = "fill") + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5,size=11),
        axis.title.y = element_text(size=11),
        axis.title.x = element_blank(),
        legend.title = element_blank())  +
  labs(y="Percentage of clusters")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
  )
ggsave("Percentage_orig.png",p,dpi=600,width = 20,height =3.5)


Idents(B_seu)<-B_seu$group
deg <- FindMarkers(
  B_seu,
  ident.1 = "CKGYF",
  ident.2 = "Con",
  logfc.threshold = 0
)

deg <- deg %>%
  tibble::rownames_to_column("gene") %>%
  as_tibble()

write.csv(deg,"B_deg.csv")


# Support Seurat outputs that use avg_log2FC or avg_logFC.
fc_col <- dplyr::case_when(
  "avg_log2FC" %in% colnames(deg) ~ "avg_log2FC",
  "avg_logFC"  %in% colnames(deg) ~ "avg_logFC",
  TRUE ~ NA_character_
)

# Standardize the fold-change field name.
deg <- deg %>%
  mutate(
    log2FC = .data[[fc_col]],
    padj = if ("p_val_adj" %in% colnames(.)) p_val_adj else NA_real_,
    pval = if ("p_val" %in% colnames(.)) p_val else NA_real_,
    neglog10padj = -log10(padj + 1e-300),
    sig = case_when(
      !is.na(padj) & padj < 0.05 & log2FC >=  0.58 ~ "Up_in_disease",
      !is.na(padj) & padj < 0.05 & log2FC <= -0.58 ~ "Down_in_disease",
      TRUE ~ "NS"
    )
  )

# Label the strongest significant upregulated and downregulated genes.
# Top upregulated genes.
up_genes <- deg %>%
  filter(sig != "NS", log2FC > 0) %>%
  arrange(desc(log2FC))%>%
  slice_head(n = 10) %>%
  pull(gene)

# Top downregulated genes.
down_genes <- deg %>%
  filter(sig != "NS", log2FC < 0) %>%
  arrange(log2FC) %>%
  slice_head(n = 10) %>%
  pull(gene)

# Combine labels.
lab_genes <- c(up_genes, down_genes)


ggplot(deg, aes(x = log2FC, y = neglog10padj)) +
  geom_point(aes(color = sig), alpha = 0.85, size = 0.6) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  ggrepel::geom_text_repel(
    data = deg %>% filter(gene %in% lab_genes),
    aes(label = gene),
    size = 3,
    max.overlaps = 10
  ) +
  scale_color_manual(values = c(
    "Up_in_disease" = "#d7301f",
    "Down_in_disease" = "#2c7fb8",
    "NS" = "grey70"
  )) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank()) +
  labs(
    title = paste0("B", ": ", "CKGYF", " vs ", "Con"),
    x = "log2 Fold Change",
    y = expression(-log[10]("adj p-value"))
  )

# Use CHIKV-upregulated genes for GO enrichment.
up_genes <- deg %>%
  filter(sig == "Up_in_disease") %>%
  pull(gene) %>%
  unique()

go_res <- NULL
p_go <- NULL
gene_map <- bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez <- unique(gene_map$ENTREZID)
go_res <- enrichGO(
      gene = entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.2,
      readable = TRUE
    )

library(clusterProfiler)
library(enrichplot)
library(dplyr)
library(ggplot2)
library(stringr)

# Continue from the previously calculated go_res object.

df <- as.data.frame(go_res) %>%
  arrange(pvalue) %>%
  slice_head(n = 21)%>%
  mutate(
    neglog10P = -log10(pvalue),
    label = paste0(ID, ": ", Description),
    label = factor(label, levels = rev(label))
  )
df<-df[-4,]

p <- ggplot(df, aes(x = neglog10P, y = label, fill = neglog10P)) +
  geom_col(width = 0.8) +
  scale_fill_gradient(low = "#F2D38A", high = "#8B1A1A") +
  labs(x = expression(-log[10](P)), y = "GO biological process terms") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 12),
    plot.margin = margin(5.5, 250, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

# Place pathway labels to the right of each bar.
p + geom_text(
  aes(x = neglog10P + max(neglog10P) * 0.02, label = label),
  hjust = 0, size = 3.5, color = "black"
) +
  # Hide y-axis labels because pathway labels are drawn directly.
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())




## Plot group comparisons.
Idents(B_seu)<-B_seu$seurat_clusters
B_seu_Plasma <- subset(B_seu ,idents = c("B_Naive"))
B_seu_Plasma <- ScaleData(B_seu_Plasma, features = genes, verbose = FALSE)

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

genes <- c("FCER2","PRDM1","AICDA","IL21R")
df <- FetchData(
  B_seu_Plasma,
  vars = c("group", genes),
  slot = "scale.data"
)

colnames(df)[1] <- "Group"

df_long <- df %>%
  pivot_longer(cols = all_of(genes), names_to = "Gene", values_to = "Expr") %>%
  mutate(
    Group = factor(Group, levels = c("Con","CKGYF"))
  )

# Pairwise comparisons.
my_comparisons <- list(
  c("Con","CKGYF")
)



p <- ggplot(df_long, aes(x = Group, y = Expr, color = Group)) +
  geom_boxplot(width = 0.4, outlier.shape = NA)  +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    p.adjust.method = "BH",
    bracket.size = 0,
    tip.length = 0,
    label = "p.signif",
    label.y.npc = "0",
    hide.ns = FALSE
  ) +
  facet_wrap(~Gene, ncol = 3, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 12)
  ) +
  labs(x = NULL, y = "Expression (log-normalized)")

p


FeaturePlot(B_seu,features = "BCL6",split.by = "group")

VlnPlot(B_seu_Plasma,features = genes,
        split.by = "group",ncol = 4,pt.size = 0.01)


genes <- c("IRF4","CD86","EBI3","DHRS9","KLF10","ZBTB32","CXCR3")
B_seu_Plasma <- subset(B_seu ,idents = c("B_Memory","B_iMemory"))
# The grouping column is expected to be named Group.
# Apply the desired order: Heal, Mild, Seve.
Idents(B_seu_Plasma) <- "group"
avg <- AverageExpression(
  B_seu_Plasma,
  features = genes,
  assays = "RNA",
  slot = "data",      # Log-normalized expression.
  group.by = "group"
)$RNA
library(dplyr)
library(tidyr)
library(ggplot2)

df <- as.data.frame(avg) %>%
  tibble::rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "group", values_to = "Expr") %>%
  mutate(
    group = factor(group, levels = c("Con","CKGYF")),
    Gene  = factor(Gene, levels = genes),
    # Match the plotting range used in the manuscript figure.
    Expr_cap = pmin(pmax(Expr, 0), 0.5)
  )
df2 <- df %>%
  group_by(Gene) %>%
  mutate(Expr_z = as.numeric(scale(Expr))) %>%
  ungroup()

ggplot(df2, aes(x = Gene, y = group, fill = Expr_z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(low = "skyblue", high = "red") +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "z-score")



## Plasma-cell differentiation genes.
Idents(B_seu)<-B_seu$seurat_clusters
genes <- c("PRDM1","XBP1","IRF4")
DotPlot(B_seu,genes,group.by = "group",cols = c("grey","red"))+
  theme_bw(base_size = 16)+
  theme(axis.title = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,linewidth = 0.8),
        axis.text = element_text(colour = "black"))







## Calculate inflammatory scores.

library(Seurat)
library(ggplot2)

# Inflammatory gene set.
inflam_genes_use
# Use normalized data from the RNA assay.
DefaultAssay(B_seu) <- "RNA"
B_seu <- AddModuleScore(
  object = B_seu,
  features = list(inflam_genes_use),
  name = "Inflammation",
  assay = "RNA",
  ctrl = 100
)

# Rename the score column for clarity.
B_seu$InflammationScore <- B_seu$Inflammation1
VlnPlot(
  B_seu,
  features = "InflammationScore",
  group.by = "seurat_clusters",
  pt.size = 0,split.by = "group"
) + ggtitle("Inflammation score by B-cell subclusters")


genes<-c("S100A8","S100A9","S100A12","IL6ST")
VlnPlot(B_seu,features = genes,
        group.by = "group",ncol = 4,pt.size = 0.001)
