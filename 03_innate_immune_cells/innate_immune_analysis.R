rm(list=ls());gc()
setwd(here::here("03_innate_immune_cells"))
seu_re <- readRDS(here::here("data", "combined_seurat.rds"))

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(tidyr)
library(readr)

markerofallsamples <- read_csv(here::here("01_initial_clustering_and_roe", "markerofallsamples.csv"))
Idents(seu_re)<-seu_re$seurat_clusters


ids <- c("NK","Others","Others","Mono_CD14","Others",
         "Mono_CD14","NK","Others","Mono_CD16","Others",
         "Others","Mono_CD16","Others","Others","Others",
         "Others","Others","Others","Mono_CD14","Others",
         "Others","DC_CD1C","Mono_NEAT1","pDC","Others"
)

#CD14,VCAN,S100A8,S100A9,FCN1
#CD1C
#GZMB,LILRA4,TCF4
#FCGR3A,LST1,AIF1
#NEAT1,NAMPT,CSF3R
#Single-cell RNA-seq reveals new types of human blood dendritic cells, monocytes, and progenitors

names(ids) <- levels(seu_re)
seu_re <- RenameIdents(seu_re, ids)

my_colors <- c(
  "NK" = "#1f77b4",
  "Mono_CD14" = "#ff7f0e",
  "Mono_CD16" = "#2ca02c",
  "DC_CD1C" = "#d62728",
  "Mono_NEAT1" = "#9467bd",
  "pDC" = "#8c564b",
  "Others" = "grey"
)

# Plot the annotated innate immune-cell UMAP.
DimPlot(seu_re, label = F, cols = my_colors)
DimPlot(seu_re, label = F, cols = my_colors,split.by = "group")+NoLegend()
Idents(seu_re)<-factor(Idents(seu_re),
                       levels=c("Others","NK","Mono_CD14","Mono_CD16",
                                "Mono_NEAT1","DC_CD1C","pDC"
                                 ))



# Save the plot at publication resolution.
DotPlot(
  subset(seu_re,idents = c("NK","Others"),invert=T),
  features = c(
    "CD14","VCAN","S100A8","S100A9","FCN1",
    "FCGR3A","LST1","AIF1",
    "NEAT1","NAMPT","CSF3R",
    "CD1C",
    "GZMB","LILRA4","TCF4"
    )
) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 10),
        axis.text.y = element_text(size=10),
        axis.title = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8)) +
  scale_color_gradient(low = "white", high = "red")

## Compare innate-cell subset abundance across conditions.
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
meta$celltype    <- factor(meta$celltype, levels = c("Others","NK","Mono_CD14","Mono_CD16",
          "Mono_NEAT1","DC_CD1C","pDC"))


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
  )+
  scale_fill_manual(values = c("Con" = "skyblue", "CKGYF" = "orange"))


# Differential expression and GO enrichment.
library(Seurat)
library(dplyr)
library(ggplot2)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

Idents(seu_re)<-seu_re$seurat_clusters
ids <- c("NK","CD8","CD4","Mono","CD4","Mono","NK","B","Mono",
         "CD8","B","Mono","γδ T","Erythroid cells","Plasma","CD8",
         "Platlet","CD4","Mono","CD4","B","DC","Mono",
         "DC","CD8")
names(ids) <- levels(seu_re)
seu_re <- RenameIdents(seu_re, ids)
seu_re$celltype<-Idents(seu_re)

de_go_volcano <- function(
    seu,
    celltype_col = "celltype",
    celltype_value = "Mono",
    group_col = "group",
    disease_group = "CKGYF",
    control_group = "Con",
    assay = NULL,                 # Use the object's default assay.
    slot = "data",                # "data" or "counts"
    test.use = "wilcox",
    min.pct = 0.1,
    logfc.threshold = 0,
    padj_cutoff = 0.05,
    lfc_cutoff = 0.58,            # Threshold for highlighting and enrichment.
    top_n_label = 10,             # Number of genes labelled on the volcano plot.
    orgdb = org.Hs.eg.db,
    keytype = "SYMBOL",           # Row names are expected to be gene symbols.
    ont = "BP",
    go_show_n = 20
) {

  md <- seu@meta.data
  seu_ct <- subset(seu, subset = celltype == celltype_value)
  Idents(seu_ct) <- seu_ct@meta.data[[group_col]]
  
  # Differential expression: CHIKV versus control.
  deg <- FindMarkers(
    seu_ct,
    ident.1 = disease_group,
    ident.2 = control_group,
    test.use = test.use,
    min.pct = min.pct,
    logfc.threshold = logfc.threshold,
    slot = slot
  )
  
  deg <- deg %>%
    tibble::rownames_to_column("gene") %>%
    as_tibble()
  
  # Support Seurat outputs that use avg_log2FC or avg_logFC.
  fc_col <- dplyr::case_when(
    "avg_log2FC" %in% colnames(deg) ~ "avg_log2FC",
    "avg_logFC"  %in% colnames(deg) ~ "avg_logFC",
    TRUE ~ NA_character_
  )
  if (is.na(fc_col)) stop("Neither avg_log2FC nor avg_logFC was found in the FindMarkers output.")
  
  # Standardize the fold-change field name.
  deg <- deg %>%
    mutate(
      log2FC = .data[[fc_col]],
      padj = if ("p_val_adj" %in% colnames(.)) p_val_adj else NA_real_,
      pval = if ("p_val" %in% colnames(.)) p_val else NA_real_,
      neglog10padj = -log10(padj + 1e-300),
      sig = case_when(
        !is.na(padj) & padj < padj_cutoff & log2FC >=  lfc_cutoff ~ "Up_in_disease",
        !is.na(padj) & padj < padj_cutoff & log2FC <= -lfc_cutoff ~ "Down_in_disease",
        TRUE ~ "NS"
      )
    )
  
  # Label the strongest significant upregulated and downregulated genes.
  # Top upregulated genes.
  up_genes <- deg %>%
    filter(sig != "NS", log2FC > 0) %>%
    arrange(desc(log2FC)) %>%
    slice_head(n = top_n_label) %>%
    pull(gene)
  
  # Top downregulated genes.
  down_genes <- deg %>%
    filter(sig != "NS", log2FC < 0) %>%
    arrange(log2FC) %>%
    slice_head(n = top_n_label) %>%
    pull(gene)
  
  # Combine labels.
  lab_genes <- c(up_genes, down_genes)
  
  
  p_volcano <- ggplot(deg, aes(x = log2FC, y = neglog10padj)) +
    geom_point(aes(color = sig), alpha = 0.85, size = 0.6) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
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
    theme(legend.title = element_blank()) +
    labs(
      title = paste0(celltype_value, ": ", disease_group, " vs ", control_group),
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
  
  if (length(up_genes) >= 10) {
    # SYMBOL -> ENTREZID
    gene_map <- suppressMessages(
      bitr(up_genes, fromType = keytype, toType = "ENTREZID", OrgDb = orgdb)
    )
    
    entrez <- unique(gene_map$ENTREZID)
    
    if (length(entrez) >= 10) {
      go_res <- enrichGO(
        gene = entrez,
        OrgDb = orgdb,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
        p_go <- barplot(go_res, showCategory = go_show_n, title = paste0(celltype_value, " Up genes: GO-", ont))
      } else {
        warning("No significant GO terms were found; there may be too few upregulated or mapped genes.")
      }
    } else {
      warning("Fewer than 10 genes mapped to ENTREZID; skipping GO enrichment.")
    }
  } else {
    warning("Fewer than 10 genes were upregulated in CHIKV; skipping GO enrichment.")
  }
  
  list(
    seu_subset = seu_ct,
    deg_table = deg,
    up_genes = up_genes,
    volcano_plot = p_volcano,
    go_object = go_res,
    go_table = if (!is.null(go_res)) as.data.frame(go_res) else NULL,
    go_plot = p_go
  )
}

res_mono <- de_go_volcano(
  seu = seu_re,
  celltype_col = "celltype",
  celltype_value = "Mono",
  group_col = "group",
  disease_group = "CKGYF",
  control_group = "Con",
  assay = "RNA",     # Change to "integrated" when appropriate.
  slot = "data",
  lfc_cutoff = 0.58,
  padj_cutoff = 0.05,
  ont = "BP",
  go_show_n = 20
)

# Display results.
res_mono$volcano_plot
res_mono$go_plot

# Export differential-expression and GO tables.
head(res_mono$deg_table)
head(res_mono$go_table)


res_nk <- de_go_volcano(
  seu = seu_re,
  celltype_col = "celltype",
  celltype_value = "NK",
  group_col = "group",
  disease_group = "CKGYF",
  control_group = "Con",
  assay = "RNA",     # Change to "integrated" when appropriate.
  slot = "data",
  lfc_cutoff = 0.58,
  padj_cutoff = 0.05,
  ont = "BP",
  go_show_n = 15
)

# Display results.
res_nk$volcano_plot
res_nk$go_plot


write.csv(res_mono$deg_table,"res_mono_deg.csv")
write.csv(res_mono$go_object,"res_mono_gobp.csv")
write.csv(res_nk$deg_table,"res_nk_deg.csv")
write.csv(res_nk$go_object,"res_nk_gobp.csv")

## Apoptosis and migration scores.
library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(tidyr)

# Retrieve human GO Biological Process gene sets.
go_bp <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
  subcategory = "BP"
)

go_apoptosis <- go_bp %>%
  filter(grepl("APOPTOTIC_PROCESS", gs_name)) %>%
  pull(gene_symbol) %>%
  unique()

go_migration <- go_bp %>%
  filter(grepl("MIGRATION|CHEMOTAXIS", gs_name)) %>%
  pull(gene_symbol) %>%
  unique()

seu_re <- AddModuleScore(
  object = seu_re,
  features = list(go_apoptosis),
  name = "Apoptosis"
)

seu_re <- AddModuleScore(
  object = seu_re,
  features = list(go_migration),
  name = "MIGRATION"
)

seu_re$seurat_clusters<-Idents(seu_re)
seu_re$group<-factor(seu_re$group,levels=c("Con","CKGYF"))
plot_df <- seu_re@meta.data %>%
  dplyr::select(
    seurat_clusters,
    group,
    Apoptosis1,
    MIGRATION1
  ) %>%
  pivot_longer(
    cols = c(Apoptosis1, MIGRATION1),
    names_to = "pathway",
    values_to = "score"
  )

plot_df <- plot_df %>%
  mutate(
    pathway = gsub("1$", "", pathway),   # Remove only the trailing "1".
    pathway = recode(
      pathway,
      Apoptosis = "Apoptosis",
      MIGRATION = "Migration"
    )
  )

plot_df<-filter(plot_df,seurat_clusters!="Others")
library(ggpubr)

ggplot(plot_df, aes(x = seurat_clusters, y = score, fill = group)) +
  geom_boxplot(
    outlier.size = 0.3,
    width = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1,color="black",size = 10),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold")
  ) +
  labs(
    x = NULL,
    y = "Pathway score"
  )+
  stat_compare_means(
    aes(group = group),
    method = "wilcox.test",
    label = "p.signif"
  )+
  scale_fill_manual(values = c(
    Con = "#3FA7A3",
    CKGYF   = "#E38D4F"
  ))+NoLegend()

## Redraw the monocyte heat map.
res_mono_deg <- read_csv(here::here("03_innate_immune_cells", "res_mono_deg.csv"))
library(readr)
library(dplyr)
library(ggplot2)
library(ggrepel)

top_up <- res_mono_deg %>%
  filter(sig == "Up_in_disease") %>%
  arrange(desc(log2FC)) %>%
  slice_head(n = 10)

top_down <- res_mono_deg %>%
  filter(sig == "Down_in_disease") %>%
  arrange(log2FC) %>%
  slice_head(n = 10)

top_genes <- bind_rows(top_up, top_down)
ggplot(res_mono_deg, aes(x = log2FC,y = neglog10padj)) +
  geom_point(
    aes(color = sig),
    size = 1.2,
    alpha = 0.7
  ) +
  scale_color_manual(
    values = c(
      "Up_in_disease" = "#D73027",
      "Down_in_disease" = "#4575B4",
      "NS" = "grey80"
    )
  ) +
  geom_vline(xintercept = c(-0.58, 0.58), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60") +
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3,
    box.padding = 0.4,
    max.overlaps = Inf
  ) +
  theme_bw() +
  labs(
    x = "log2 Fold Change",
    y = "-log10(adj. P value)",
    color = "",
    title = "Volcano plot of Mono DEGs"
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )





res_nk_gobp <- res_nk_gobp %>%
  arrange(p.adjust) %>%
  mutate(
    Description = factor(Description, levels = rev(Description))
  )

ggplot(res_nk_gobp[c(1,4,7,13,14,16,19,20),],
       aes(x = Count,
           y = Description,
           fill = -log10(p.adjust))) +
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


