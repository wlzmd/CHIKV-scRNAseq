rm(list=ls());gc()
setwd(here::here("07_global_cellchat_analysis"))
combined_seurat <- readRDS(here::here("data", "combined_seurat.rds"))

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
group.by      <- "group"
celltype.col  <- "seurat_clusters"

## Condition labels used in the Seurat metadata.
group.disease <- "CKGYF"
group.healthy <- "Con"

## Use the human CellChat database.
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
  
  # Use the expression matrix from the current default assay.
  if (is.null(assay.use)) assay.use <- DefaultAssay(seu.g)
  data.input <- GetAssayData(seu.g, assay = assay.use, slot = "data")
  
  # Prepare metadata and cell-group labels.
  meta <- seu.g@meta.data
  stopifnot(celltype.col %in% colnames(meta))
  meta$cellchat_celltype <- meta[[celltype.col]]
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cellchat_celltype")
  cellchat@DB <- db
  
  # Retain secreted, ECM, and cell-cell contact interactions.
  # cellchat@DB <- subsetDB(db, search = "Secreted Signaling")
  
  # Preprocess expression and interaction data.
  cellchat <- subsetData(cellchat)                  # subset to signaling genes
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  # Estimate communication probabilities.
  cellchat <- projectData(cellchat, PPI.human)
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)
  
  # Summarize communication at the pathway level.
  cellchat <- computeCommunProbPathway(cellchat)
  
  # Aggregate the communication network.
  cellchat <- aggregateNet(cellchat)
  
  # Calculate outgoing and incoming signalling centrality.
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
object.list <- list(CKGYF = cellchat.d, Con = cellchat.h)
cellchat.merged <- mergeCellChat(object.list, add.names = names(object.list))


## A: signaling role scatter（Outgoing vs Incoming）
pA1 <- netAnalysis_signalingRole_scatter(cellchat.d, title = group.disease)
pA2 <- netAnalysis_signalingRole_scatter(cellchat.h, title = group.healthy)

pA2 + pA1 


## Inspect cell-type order before selecting a focal population.
levels(cellchat.d@idents)

focus.cell <- "Mono"

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
# Optionally restrict receiver populations.
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
# Alternative view.
rankNet(cellchat.merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)


## Inspect available pathways.
cellchat.d@netP$pathways %>% head()

netAnalysis_contribution(cellchat.d, signaling = "MHC-I")
netAnalysis_contribution(cellchat.d, signaling = "MHC-II")
netAnalysis_contribution(cellchat.d, signaling = "MIF")
netAnalysis_contribution(cellchat.d, signaling = "CD99")


