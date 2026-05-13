rm(list=ls());gc()
setwd("~/SHIXI/3.基肯孔雅热/07_全局cellchat分析")
combined_seurat <- readRDS("~/SHIXI/3.基肯孔雅热/0.0rawdata/combined_seurat.rds")

DimPlot(combined_seurat)

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(patchwork)
  library(ggplot2)
  library(future)
})

combined_seurat$seurat_clusters<-Idents(combined_seurat)

## =========================================================
## 1) Basic settings: please edit these two if needed
## =========================================================
group.by      <- "group"     # 分组列名：疾病/健康
celltype.col  <- "seurat_clusters"  # 细胞类型列名（大类群）

## 你两组的名字（按你对象里的实际值写）
group.disease <- "CKGYF"  # 或 "Disease"
group.healthy <- "Con"   # 或 "Healthy"

## 使用人类数据库（如果是小鼠改为 CellChatDB.mouse）
CellChatDB.use <- CellChatDB.human

## =========================================================
## 2) Helper: build CellChat object for one group
## =========================================================
run_cellchat_one_group <- function(seu, group.name,
                                   group.by = "group",
                                   celltype.col = "seurat_clusters",
                                   db = CellChatDB.human,
                                   assay.use = NULL,  # NULL=DefaultAssay(seu)
                                   min.cells = 10) {
  
  seu.g <- subset(seu, subset = group == group.name)
  
  # 选择表达矩阵（默认用对象当前 DefaultAssay）
  if (is.null(assay.use)) assay.use <- DefaultAssay(seu.g)
  data.input <- GetAssayData(seu.g, assay = assay.use, slot = "data")
  
  # meta + 分群
  meta <- seu.g@meta.data
  stopifnot(celltype.col %in% colnames(meta))
  meta$cellchat_celltype <- meta[[celltype.col]]
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_celltype")
  cellchat@DB <- db
  
  # 可选：过滤数据库类别（这里保留所有：Secreted/ECM/Cell-Cell contact）
  # cellchat@DB <- subsetDB(db, search = "Secreted Signaling")
  
  # 预处理
  cellchat <- subsetData(cellchat)                  # subset to signaling genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # 计算通讯概率
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  # pathway层面
  cellchat <- computeCommunProbPathway(cellchat)
  
  # 聚合网络
  cellchat <- aggregateNet(cellchat)
  
  # 计算中心性（用于 A 图那种 outgoing/incoming role）
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  
  return(cellchat)
}


## =========================================================
## 3) Run CellChat for Disease & Healthy
## =========================================================
cellchat.d <- run_cellchat_one_group(
  seu = combined_seurat,
  group.name = group.disease,
  group.by = group.by,
  celltype.col = celltype.col,
  db = CellChatDB.use
)

cellchat.h <- run_cellchat_one_group(
  seu = combined_seurat,
  group.name = group.healthy,
  group.by = group.by,
  celltype.col = celltype.col,
  db = CellChatDB.use
)


saveRDS(cellchat.d, "cellchat_disease.rds")
saveRDS(cellchat.h, "cellchat_healthy.rds")

## =========================================================
## 4) Compare: merge two CellChat objects
## =========================================================
object.list <- list(GKCYF = cellchat.d, Con = cellchat.h)  # 名字会显示在图例中
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))


## A: signaling role scatter（Outgoing vs Incoming）
pA1 <- netAnalysis_signalingRole_scatter(cellchat.d, title = group.disease)
pA2 <- netAnalysis_signalingRole_scatter(cellchat.h, title = group.healthy)

pA2 + pA1 


## 先看细胞类型顺序（确保名字匹配）
levels(cellchat.d@idents)

focus.cell <- "Mono"   # 改成你对象里的真实名字，如 "Mono" / "CD14+ Mono"

## B-left: outgoing from focus.cell
netVisual_circle(cellchat.d@net$weight,
                 sources.use = focus.cell,label.edge = TRUE,
                 title.name = paste0(group.disease, " - outgoing from ", focus.cell))

netVisual_circle(cellchat.h@net$weight,
                 sources.use = focus.cell,label.edge = TRUE,
                 title.name = paste0(group.healthy, " - outgoing from ", focus.cell))

## B-right: incoming to focus.cell
netVisual_circle(cellchat.d@net$weight,
                 targets.use = focus.cell,label.edge = TRUE,
                 title.name = paste0(group.disease, " - incoming to ", focus.cell))

netVisual_circle(cellchat.h@net$weight,
                 targets.use = focus.cell,label.edge = TRUE,
                 title.name = paste0(group.healthy, " - incoming to ", focus.cell))




## C: bubble plot for LR pairs (Monocytes -> others)
# 你可以限制 receiver 细胞类群（按需改）
receivers <- c("B","CD4","CD8","NK","DC","γδ T","Mono", "Plasma")

pC_ATB <- netVisual_bubble(cellchat.d,
                           sources.use = focus.cell,
                           targets.use = receivers,
                           remove.isolate = TRUE,
                           title.name = paste0(group.disease, ": ", focus.cell, " -> targets"))

pC_HC  <- netVisual_bubble(cellchat.h,
                           sources.use = focus.cell,
                           targets.use = receivers,
                           remove.isolate = TRUE,
                           title.name = paste0(group.healthy, ": ", focus.cell, " -> targets"))

pC_ATB + pC_HC

library(ggplot2)
ggsave("pC_ATB_HC.png", plot = pC_HC + pC_ATB , width = 8, height = 9, dpi = 600)

## D: compare information flow / signaling pathways between conditions
rankNet(cellchat.merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
# 或者
rankNet(cellchat.merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)


## 查看可用 pathways
cellchat.d@netP$pathways %>% head()

netAnalysis_contribution(cellchat.d, signaling = "MHC-I")
netAnalysis_contribution(cellchat.d, signaling = "MHC-II")
netAnalysis_contribution(cellchat.d, signaling = "MIF")
netAnalysis_contribution(cellchat.d, signaling = "CD99")


