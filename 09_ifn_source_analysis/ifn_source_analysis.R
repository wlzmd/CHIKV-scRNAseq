setwd(here::here("09_ifn_source_analysis"))
library(Seurat)
library(ggplot2)

combined_seurat_chikv <- readRDS(here::here("data", "combined_seurat.rds"))

DimPlot(combined_seurat_chikv)
infg<-c("IFNB1","IFNL1","IFNL2","IFNL3","IFNG")
DotPlot(combined_seurat_chikv,infg)+theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,        # Rotate x-axis labels by 90 degrees.
                                   hjust = 1,         # Right-align labels.
                                   face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank()) +
  scale_color_gradient(low = "grey", high = "orange")
  
combined_seurat_chikv$group<-factor(combined_seurat_chikv$group,levels=c("Con","CKGYF"))
VlnPlot(subset(combined_seurat_chikv,idents=c("DC")
               ),"IFNL1",split.by="group")+NoLegend()+
  theme(axis.title = element_blank())

     
