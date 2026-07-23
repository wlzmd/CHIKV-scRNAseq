setwd(here::here("06_t_cell_subclustering"))
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
DimPlot(seu_re)

CD8 <-subset(seu_re,idents = c("CD8")) %>%
  FindVariableFeatures()%>%
  RunPCA(npcs = 30) 
CD8 <- CD8 %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:20)
CD8 <- FindClusters(CD8,resolution = 0.3,verbose = F)
DimPlot(CD8,label=T)
CD8Marker<-FindAllMarkers(CD8,logfc.threshold = 0.58,only.pos = T)

FeaturePlot(CD8,"FOXP3")
Idents(CD8)<-CD8$RNA_snn_res.0.3
ids <- c("CD8_Effector-GNLY","CD8_Effector-GZMK","CD8_Naive",
   "CD4_Effector-GNLY","CD8_Pro","NKT_CD160","NKT_Naive",
   "CD8_Effector-GZMK","NKT_CD56","CD8_Pro")
names(ids)<-levels(CD8)
CD8<-RenameIdents(CD8,ids) 



CD4 <- subset(seu_re,idents = c("CD4")) %>%
   FindVariableFeatures()%>%
   RunPCA(npcs = 30) 
CD4 <- CD4 %>%
   FindNeighbors(dims = 1:15) %>%
   RunUMAP(dims = 1:20)
CD4 <- FindClusters(CD4,resolution = 0.3,verbose = F)
library(RColorBrewer)
DimPlot(CD4, label = TRUE, cols = brewer.pal(8, "Set2"))
CD4Marker<-FindAllMarkers(CD4,logfc.threshold = 0.58,only.pos = T)
FeaturePlot(CD4,"GPR183")

Idents(CD4)<-CD4$RNA_snn_res.0.3
ids <- c("CD4_Naive","CD4_Th2","CD4_eMemory",
         "CD4_Naive","CD4_Treg","CD4_Th1","CD4_Memory",
         "CD4_Naive","CD4_Tfh")
names(ids)<-levels(CD4)
CD4<-RenameIdents(CD4,ids)
 

CD8$seurat_clusters <- c("CD8_Effector-GNLY","CD8_Effector-GZMK","CD8_Naive",
  "CD4_Effector-GNLY","CD8_Pro","","",
  "CD8_Effector-GZMK","","CD8_Pro")
DimPlot(CD4,label = T)
DimPlot(CD4, label = F, cols = brewer.pal(8, "Set2"))+NoLegend()+NoAxes()





CD8_sub <- subset(CD8,ident="CD4_Effector-GNLY",invert=T) %>%
  RunUMAP(dims = 1:20)
CD8 <- FindClusters(CD8,resolution = 0.3,verbose = F)
DimPlot(CD8_sub,label=F,split.by = "group")+NoAxes()+NoLegend()
Idents(CD8_sub)<-factor(Idents(CD8_sub),
c("CD8_Naive","CD8_Effector-GNLY","CD8_Effector-GZMK",
"CD8_Pro","NKT_Naive","NKT_CD160",
"NKT_CD56"))
CD8_sub$seurat_clusters<-Idents(CD8_sub)
genes <- c("CD3G","CD3E","CD3D","CD8A","CD8B","SELL","CCR7",
           "GZMK","GNLY","MKI67","NCAM1","CD160")

png("CD8_FeaturePlot.png", width = 1200, height = 1000)
FeaturePlot(CD8_sub, features = genes,ncol = 4)
dev.off()

CD4_sub<-merge(subset(CD8,ident="CD4_Effector-GNLY"),CD4)
CD4_sub<-JoinLayers(CD4_sub)
CD4_sub<-CD4_sub %>%
  FindVariableFeatures()%>%
  RunPCA(npcs = 30) %>%
  RunUMAP(dims = 1:20)
DimPlot(CD4_sub)
Idents(CD4_sub)<-factor(Idents(CD4_sub),
                        c("CD4_Naive","CD4_Effector-GNLY","CD4_eMemory",
                          "CD4_Memory","CD4_Tfh","CD4_Th1","CD4_Th2",
                          "CD4_Treg"))


table(Idents(CD4_sub))
genes <- c("CD3E","CD3D","CD40LG","SELL","CCR7",
           "GPR183","S100A4","GZMA","ICOS","SLAMF1",
           "TBX21","MAF","GNLY","FOXP3")

png("CD4_FeaturePlot.png", width = 1200, height = 1300)
FeaturePlot(CD4_sub, features = genes,ncol = 4)
dev.off()
DimPlot(CD4_sub,label=F,split.by = "group")+NoAxes()+NoLegend()


CD4Marker<-FindAllMarkers(CD4_sub,logfc.threshold = 0.58,only.pos = T)
CD8Marker<-FindAllMarkers(CD8_sub,logfc.threshold = 0.58,only.pos = T)


saveRDS(CD4_sub,"CD4_seu.RDS")
saveRDS(CD8_sub,"CD8_seu.RDS")

write_csv(CD4Marker,"CD4Marker.csv")
write_csv(CD8Marker,"CD8Marker.csv")


##CD4 CD8
CD4_seu$seurat_clusters<-Idents(CD4_seu)
Tseu  <-  merge(CD4_seu, y = CD8_seu, 
                      add.cell.ids = c("CD4", "CD8"), 
                      project = "Tcell")
Tseu$group<-factor(Tseu$group,levels=c("Con","CKGYF"))

colnames(Tseu@meta.data)
library(dplyr)
library(tidyr)

meta <- Tseu@meta.data %>%
  dplyr::select(orig.ident, group, seurat_clusters)

# Calculate each subset as a percentage of T cells within each sample.
df_prop <- meta %>%
  group_by(orig.ident, group, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident,group) %>%
  mutate(percentage = 100 * n / sum(n)) %>%
  ungroup()

library(ggplot2)

p <- ggplot(df_prop, aes(x = group, y = percentage, fill = group)) +
  stat_summary(fun = mean, geom = "bar", width = 0.75, color = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_point(aes(group = orig.ident), position = position_jitter(width = 0.15), size =0.5 ) +
  facet_wrap(~ seurat_clusters, scales = "free_y", ncol = 5) +
  labs(x = NULL, y = "Percentage") +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave("t_cell_subset_proportions.png", plot = p, width = 10, height = 6, dpi = 300)



##
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)
unique(Tseu$seurat_clusters)


celltypes_use <- c("CD4_Effector-GNLY",
                   "CD4_eMemory",
                   "CD8_Effector-GZMK",
                   "CD8_Effector-GNLY",
                   "NKT_CD56",
                   "NKT_CD160")

# Cell-type annotations are expected in Tseu$celltype.
Tsub <- subset(Tseu, subset = seurat_clusters %in% celltypes_use)

library(clusterProfiler)
library(org.Hs.eg.db)

go_cytotox <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0001913",
  keytype = "GOALL",
  columns = c("SYMBOL")
)

cytotox_genes_go <- unique(na.omit(go_cytotox$SYMBOL))
length(cytotox_genes_go)

go_exhaust <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = "GO:0002764",
  keytype = "GOALL",
  columns = c("SYMBOL")
)

exhaust_genes_go <- unique(na.omit(go_exhaust$SYMBOL))
length(exhaust_genes_go)


DefaultAssay(Tsub) <- "RNA"

Tsub <- AddModuleScore(Tsub,
                       features = list(cytotox_genes_go),
                       name = "Cytotoxicity_GO")

Tsub <- AddModuleScore(Tsub,
                       features = list(exhaust_genes_go),
                       name = "Exhaustion_GO")

df <- Tsub@meta.data %>%
  dplyr::select(group,seurat_clusters, Cytotoxicity_GO1, Exhaustion_GO1) %>%
  tidyr::pivot_longer(cols = c(Cytotoxicity_GO1, Exhaustion_GO1),
                      names_to = "Signature",
                      values_to = "Score")

df$celltype <- factor(df$seurat_clusters, levels = celltypes_use)
df$Signature <- recode(df$Signature,
                       Cytotoxicity1 = "Cytotoxicity",
                       Exhaustion1 = "Exhaustion")
# Condition order used in the figure.
df$group <- factor(df$group , levels = c("Con","CKGYF"))

p <- ggplot(df, aes(x = celltype, y = Score, fill = group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(width = 0.6)) +
  facet_wrap(~ Signature, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(x = NULL, y = NULL)

p + stat_compare_means(comparisons = my_comp, method = "wilcox.test", label = "p.format", size = 3) +
  labs(title = NULL) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 10),
         axis.text.y = element_text(color = "black", size = 10),
         strip.text = element_blank()
         )

##DotPlot
celltype_use <- "NKT_CD56"
Tsub <- subset(Tseu, ident = celltype_use)
colnames(Tseu@meta.data)
unique(Tseu$seurat_clusters)
genes_use <- c("CCR7","SELL","LEF1","TCF7",
               "PRF1","IFNG","GNLY","NKG7",
               "GZMB","GZMA","GZMH",
               "KLRK1","KLRD1","KLRC1",
               "CTSW","CST7","FCGR3A",
               "LAG3","TIGIT","PDCD1","CTLA4","HAVCR2","TOX")

p <- DotPlot(
  Tsub,
  features = genes_use,
  group.by = "group",
  dot.scale = 8
) +
  RotatedAxis() +
  scale_color_gradientn(colors = c("navy","white","gold"))

p



##DEGs
Idents(Tseu)<-Tseu$group
Tseu<-JoinLayers(Tseu)
deg <- FindMarkers(
  Tseu,
  ident.1 = "CKGYF",
  ident.2 = "Con",
  logfc.threshold = 0
)

deg <- deg %>%
  tibble::rownames_to_column("gene") %>%
  as_tibble()

write.csv(deg,"T_deg.csv")


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
    max.overlaps = Inf
  ) +
  scale_color_manual(values = c(
    "Up_in_disease" = "#d7301f",
    "Down_in_disease" = "#2c7fb8",
    "NS" = "grey70"
  )) +
  theme_classic(base_size = 12) +
  theme(legend.title = element_blank(),
        legend.position = "None") +
  labs(
    title = paste0("T", ": ", "CKGYF", " vs ", "Con"),
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
  mutate(
    neglog10P = -log10(pvalue),
    label = paste0(ID, ": ", Description),
    label = factor(label, levels = rev(label))
  )

df<-df[c("GO:0034341","GO:0034340","GO:0035456","GO:0035455","GO:0001819",
         "GO:0001906","GO:1903037","GO:0061640","GO:0051607","GO:0001819"),] 
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




##
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)

# Migration: leukocyte migration (GO:0050900).
go_mig <- AnnotationDbi::select(org.Hs.eg.db, keys="GO:0050900", keytype="GOALL", columns="SYMBOL")
mig_genes <- unique(na.omit(go_mig$SYMBOL))

# Apoptosis: apoptotic process (GO:0006915).
go_apo <- AnnotationDbi::select(org.Hs.eg.db, keys="GO:0006915", keytype="GOALL", columns="SYMBOL")
apo_genes <- unique(na.omit(go_apo$SYMBOL))

# Retain genes present in the Seurat object.
mig_genes <- intersect(mig_genes, rownames(Tseu))
apo_genes <- intersect(apo_genes, rownames(Tseu))

DefaultAssay(Tseu) <- "RNA"
Tseu <- AddModuleScore(Tseu, features=list(mig_genes), name="Migration")
Tseu <- AddModuleScore(Tseu, features=list(apo_genes), name="Apoptosis")


df_sample <- Tseu@meta.data %>%
  dplyr::select(orig.ident, group, Migration1, Apoptosis1) %>%
  group_by(orig.ident, group) %>%
  summarise(Migration = mean(Migration1),
            Apoptosis = mean(Apoptosis1),
            .groups="drop") %>%
  pivot_longer(cols=c(Migration, Apoptosis),
               names_to="ScoreType",
               values_to="Score")

df_sample$Condition <- factor(df_sample$group, levels=c("Con","CKGYF"))
df_sample$ScoreType <- factor(df_sample$ScoreType, levels=c("Migration","Apoptosis"))

ggplot(df_sample, aes(x=Condition, y=Score, fill=Condition)) +
  geom_boxplot(width=0.4, outlier.shape = NA) +
  geom_jitter(width=0.15, size=0.8, alpha=0.5) +
  facet_wrap(~ScoreType, ncol=1, scales="free_y") +
  theme_classic() +
  labs(x=NULL, y=" ")+
  labs(title = NULL) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12),
         axis.text.y = element_text(color = "black", size = 12),
         strip.text = element_blank(),
         panel.spacing = unit(1, "cm")
  )


df_sample_ct <- Tseu@meta.data %>%
  dplyr::select(orig.ident, group, seurat_clusters, Migration1, Apoptosis1) %>%
  group_by(orig.ident, group, seurat_clusters) %>%
  summarise(Migration = mean(Migration1),
            Apoptosis = mean(Apoptosis1),
            .groups="drop") %>%
  pivot_longer(cols=c(Migration, Apoptosis),
               names_to="ScoreType",
               values_to="Score")

df_sample_ct$Condition <- factor(df_sample_ct$group, levels=c("Con","CKGYF"))
unique(df_sample_ct$seurat_clusters)
celltype_order <- c("CD4_Naive","CD4_Memory","CD4_eMemory","CD4_Effector-GNLY",
                    "CD4_Th1","CD4_Th2","CD4_Tfh","CD4_Treg",
                    "CD8_Naive","CD8_Effector-GZMK","CD8_Effector-GNLY","CD8_Pro",
                    "NKT_Naive","NKT_CD56","NKT_CD160")
df_sample_ct$seurat_clusters <- factor(df_sample_ct$seurat_clusters, levels = celltype_order)



ggplot(df_sample_ct, aes(x = seurat_clusters, y = Score, fill = Condition)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.2, alpha = 0.8) +
  facet_wrap(~ScoreType, ncol = 1, scales = "free_y") +
  theme_classic()+
  labs(x = NULL, y = NULL)+
  labs(title = NULL) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) 


ggplot(df_sample_ct, aes(x=seurat_clusters, y=Score, fill=Condition)) +
  geom_boxplot(width = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8))+
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
             size = 0.2, alpha = 0.8) +
  facet_wrap(~ScoreType, ncol=1, scales="free_y") +
  theme_classic() +
  labs(x=NULL, y=" ")+
  labs(title = NULL) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1, color = "black", size = 12),
         axis.text.y = element_text(color = "black", size = 12),
         strip.text = element_blank(),
         panel.spacing = unit(1.3, "cm")
  )


##NKT_CD56
CD8Marker <- read_csv(here::here("06_t_cell_subclustering", "CD8Marker.csv"))
markersgenes <-filter(CD8Marker,cluster == "NKT_CD56"& pct.1 > 0.20 )$gene

library(clusterProfiler)
library(org.Hs.eg.db)
ego<-enrichGO(
  gene = markersgenes,
  OrgDb = org.Hs.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

head(ego)
ego<-data.frame(ego@result)
ego_sub<-ego[1:10,]


ego_sub <- ego_sub %>%
  arrange(p.adjust) %>%
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

library(ggplot2)
ego_sub$logP<- -log10(ego_sub$p.adjust)

ggplot(ego_sub,
       aes(x = Count,
           y = Description,
           fill = logP)) +
  geom_col(width = 0.7) +
  
  scale_fill_gradient(
    low = "#2166AC",
    high = "#B2182B"
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(color = "black",size = 12))+
  labs(
    x = "Gene number",
    y = NULL,
    fill = "-Log10(adj.p)"
  )


