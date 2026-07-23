rm(list=ls());gc()
setwd(here::here("02_ifna_and_inflammation_scores"))
seu_re <- readRDS(here::here("data", "combined_seurat.rds"))

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

unique(go_bp$gs_name[grepl("INTERFERON.*ALPHA|INTERFERON_ALPHA", go_bp$gs_name)])
unique(go_bp$gs_name[grepl("ACUTE.*INFLAM", go_bp$gs_name)])


# 1. IFN-alpha response
genes_ifna <- go_bp %>%
  filter(gs_name == "GOBP_RESPONSE_TO_INTERFERON_ALPHA") %>%
  pull(gene_symbol) %>%
  unique()

# 2. Acute inflammatory response
genes_acute_inflam <- go_bp %>%
  filter(gs_name == "GOBP_ACUTE_INFLAMMATORY_RESPONSE") %>%
  pull(gene_symbol) %>%
  unique()

length(genes_ifna)
length(genes_acute_inflam)

# The Seurat object is expected to have been normalized and scaled.

seu_re <- AddModuleScore(
  object = seu_re,
  features = list(genes_ifna),
  name = "IFNa_Response"
)

seu_re <- AddModuleScore(
  object = seu_re,
  features = list(genes_acute_inflam),
  name = "Acute_Inflammatory_Response"
)

seu_re$seurat_clusters<-Idents(seu_re)
seu_re$group<-factor(seu_re$group,levels=c("Con","CKGYF"))
plot_df <- seu_re@meta.data %>%
  dplyr::select(
    seurat_clusters,
    group,
    IFNa_Response1,
    Acute_Inflammatory_Response1
  ) %>%
  pivot_longer(
    cols = c(IFNa_Response1, Acute_Inflammatory_Response1),
    names_to = "pathway",
    values_to = "score"
  )

plot_df <- plot_df %>%
  mutate(
    pathway = gsub("1$", "", pathway),   # Remove only the trailing "1".
    pathway = recode(
      pathway,
      IFNa_Response = "Response to IFN-α",
      Acute_Inflammatory_Response = "Acute inflammatory response"
    )
  )



ggplot(plot_df, aes(x = seurat_clusters, y = score, fill = group)) +
  geom_boxplot(
    outlier.size = 0.3,
    width = 0.7,
    position = position_dodge(width = 0.8)
  ) +
  facet_wrap(~ pathway, ncol = 1, scales = "free_y") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
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
  ))



