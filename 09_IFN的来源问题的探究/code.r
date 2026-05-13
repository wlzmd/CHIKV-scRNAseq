DimPlot(combined_seurat_chikv)
infg<-c("IFNB1","IFNL1","IFNL2","IFNL3","IFNG")
DotPlot(combined_seurat_chikv,infg)+theme(axis.text.x = element_text(angle = 90))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90,        # 横坐标旋转 90°
                                   hjust = 1,         # 对齐方式
                                   face = "italic"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.line = element_blank()) +  # 字体加粗+斜体
  scale_color_gradient(low = "grey", high = "orange") # 颜色从灰到橘红
  
combined_seurat_chikv$group<-factor(combined_seurat_chikv$group,levels=c("Con","CKGYF"))
VlnPlot(subset(combined_seurat_chikv,idents=c("DC")
               ),"IFNL1",split.by="group")+NoLegend()+
  theme(axis.title = element_blank())

     
     