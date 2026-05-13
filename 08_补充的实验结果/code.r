setwd("~/SHIXI/3.基肯孔雅热/08_B细胞的时序分析")
CD8 <- readRDS("~/SHIXI/3.基肯孔雅热/06_T细胞的聚类分群/CD8_seu.RDS")
FeaturePlot(CD8,"ZBTB16")


B_seu <- readRDS("~/SHIXI/3.基肯孔雅热/05_B细胞聚类分群/B_seu.RDS")
plasma <- subset(B_seu, idents = "B_Plasma")
plasma <- subset(B_seu, subset = group == "Con")
ig_genes <- c("IGHM", "IGHG1", "IGHG2", "IGHG3",
              "IGHG4", "IGHA1", "IGHA2")

expr <- GetAssayData(plasma, slot = "data")[ig_genes, ]
dominant_isotype <- apply(expr, 2, function(x) {
  ig_genes[which.max(x)]
})

plasma$Ig_isotype <- dominant_isotype
library(dplyr)

ig_prop <- as.data.frame(table(plasma$Ig_isotype)) %>%
  mutate(Percent = Freq / sum(Freq) * 100)

ig_prop <- ig_prop %>%
  arrange(desc(Var1)) %>%
  mutate(
    ypos = cumsum(Percent) - 0.5 * Percent
  )
library(ggplot2)

library(RColorBrewer)

)

ggplot(ig_prop, aes(x = "", y = Percent, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  geom_text(
    aes(y = ypos, x = 1.8,
        label = paste0(Var1, " ", round(Percent, 2), "%")),
    size = 4
  )




ggplot(ig_prop, aes(x = "", y = Percent, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  
  ## ① 折线（从饼图边缘到文字）
  geom_segment(
    aes(x = 1, xend = 1.7, y = ypos, yend = ypos),
    color = "grey40", linewidth = 0.4
  ) +
  
  ## ② 圈外标签
  geom_text(
    aes(x = 1.8, y = ypos,
        label = paste0(Var1, " ", round(Percent, 2), "%")),
    size = 4, hjust = 0
  )
library(ggrepel)

ggplot(ig_prop, aes(x = "", y = Percent, fill = Var1)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y") +
  theme_void() +
  scale_fill_brewer(palette = "Set3") +
  geom_text_repel(
    aes(label = paste0(Var1, " ", round(Percent, 2), "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    show.legend = FALSE
  )

