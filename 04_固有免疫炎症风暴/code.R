rm(list=ls());gc()
setwd("~/SHIXI/3.基肯孔雅热/04_固有免疫炎症风暴")
seu_re <- readRDS("~/SHIXI/3.基肯孔雅热/0.0rawdata/combined_seurat.rds")

library(Seurat)
library(dplyr)
library(ggplot2)
library(msigdbr)
library(tidyr)

markerofallsamples <- read_csv("~/SHIXI/3.基肯孔雅热/01_粗聚类分群并检测ROE值/markerofallsamples.csv")
Idents(seu_re)<-seu_re$RNA_snn_res.0.3

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

library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(UCell)

cytokine_genes <- c(
  "IL2", "IL7", "CSF3", "CXCL10", "CCL2", "IL15", "TNF", "IL6", "CXCL7", "IL1RN",
  "CSF1", "IFNB", "IL20RA", "IL10", "IL18", "HGF", "CXCL9", "CCL27", "TGFB1", "IL1B",
  "LTA", "CSF2", "IL12", "TNFSF13", "IL11", "CCL12", "CXCL8", "CXCL11", "CCL4", "CXCL1",
  "CXCL2", "CCL3", "CCL18", "CXCL16", "FGA", "CCL5", "IL1A", "CCL20", "OSM", "TNFSF14",
  "FGF10", "CXCL5", "IL9", "IL1C", "TNFSF10", "FGF7", "IL23", "TGFB3", "TNFRSF1B", "CREB1",
  "TXN3A", "IL1W", "XCL2", "XCL1", "CCL25", "CCL16", "CD40LG", "IL20", "FASLG", "F10",
  "SCYL1", "GDF15", "TNFSF8", "GDF11", "VEGFA", "GDF7", "PDGFA", "IL21", "ADGRL1", "ADGRD2",
  "PDGFB", "TNFAIP6", "TAMPAI", "BMP6", "PDGFD", "IL20RB", "GHR", "WISP3", "MIP", "IL4",
  "BMP8B", "TNFSF12", "SCYL2", "TSP1", "SDF1B", "NMBR", "PRP", "FGF11", "PGR-AS1", "FGF2",
  "VEGFC", "IL12A", "BAFF", "AMH", "IL26", "PDGFC", "FGF23", "KLF11", "TNFSF9", "BMP7",
  "IL24", "BMP10", "CXCL6", "CD96", "IL16", "CTD", "FGF9", "PSNL1", "TSC1", "IL32", "IL10R",
  "SPPL", "IL12RB1", "BMP4", "TARBP1", "TGFBR2", "A20MT", "EDA", "TNFRSF11B", "CCL28"
)

inflam_genes   <- c(
  "ABCA1", "ABI3", "ACKR1", "ACKR3", "ADORA2A", "ADORA2B", "AHR", "AIM2", "AKIRIN2", "AP2A2",
  "ATF3", "ATP6V1H", "BCL2A1", "BIRC3", "BST1", "BST2", "BTK", "CARD11", "CASP4", "CCL17",
  "CCL2", "CCL22", "CCL3", "CCL4", "CCL5", "CCL7", "CCR1", "CCR2", "CCR5", "CD14",
  "CD163", "CD19", "CD3D", "CD3G", "CD4", "CD74", "CENPK", "CHST2", "CHST7", "CHURC1",
  "CIITA", "CSF1", "CSF2", "CSF3", "CXCL10", "CXCL11", "CXCL6", "CXCL9", "CYBB", "CYLD",
  "DCLRE1B", "DUSP1", "EGR1", "EPS8L2", "ERO1A", "F3", "FADD", "FAS", "FBP1", "FGR",
  "FTH1", "G0S2", "GABARAPL1", "GAS1", "GAS6", "GATA3", "GBP1", "GBP1P1", "GBP2", "GBP3",
  "GBP4", "GBP5", "GBP6", "GBP7", "GPR183", "GPR84", "HCK", "HIF1A", "HLA-DRA", "ICAM1",
  "ICOS", "IFI16", "IFI30", "IFI44", "IFI44L", "IFI6", "IFIT1", "IFIT2", "IFIT3", "IFITM1",
  "IFITM2", "IFITM3", "IGFBP3", "IL10", "IL10RA", "IL10RB", "IL15", "IL15RA", "IL18", "IL1R1",
  "IL1R2", "IL1RAP", "IL1RN", "IL2RA", "IL4R", "IL6", "IL6R", "IL7R", "IL8", "IRAK3",
  "IRF1", "IRF7", "ITGAX", "ITGB2", "ITGB7", "KCNA3", "KCNK10", "KCNK3", "KLF6", "LAMP3",
  "LAPTM5", "LCN2", "LCP2", "LTA", "LY86", "LYN", "MAP3K8", "MAF", "MAFB", "MARCO",
  "MEFV", "MEP1A", "MEP1B", "MERTK", "MMP14", "MMP9", "MNDA", "MRC1", "MSR1", "NAB2",
  "NCR1", "NFKB1", "NFKB2", "NLRP1", "NLRP3", "NOS2", "NQO1", "NQO2", "OAS1", "OAS2",
  "OAS3", "ORM1", "ORM2", "P2RX7", "P2RY6", "P2RY13", "P2RY14", "PADI4", "PARK7", "PBX1",
  "PDCD1", "PDCD1LG2", "PDE4B", "PDE4D", "PGLYRP1", "PIK3R5", "PLAU", "PLSCR1", "PSEN1", "PTAFR",
  "PTGER2", "PTGER4", "PTGS2", "PTPN6", "PYCARD", "RIPK1", "RIPK3", "RNASE2", "RPL13A", "S100A8",
  "S100A9", "S100B", "SAMD9L", "SASH3", "SBNO2", "SELL", "SEMA4D", "SERPINB1", "SIGLEC1", "SIGLEC14",
  "SIGLEC5", "SIGLEC7", "SIGLEC9", "SLC11A1", "SLC15A3", "SLC2A3", "SLC39A8", "SLC7A2", "SMAD3", "SOCS1",
  "SOCS3", "SPI1", "SRGN", "STAT1", "TACR1", "TACR3", "TAPBP", "TBX21", "TBX2", "TLR1",
  "TLR2", "TLR3", "TLR4", "TLR5", "TLR6", "TLR7", "TLR8", "TNFAIP3", "TNFRSF1B", "TNFRSF9",
  "TNFSF13", "TNFSF13B", "TNFSF14", "TNFSF15", "TNF", "TYROBP"
)


length(cytokine_genes); head(cytokine_genes)
length(inflam_genes);   head(inflam_genes)

DefaultAssay(seu_re) <- "RNA"  # 如果你用 SCT 就换成 "SCT"

present_genes <- rownames(seu_re)

cytokine_genes_use <- intersect(cytokine_genes, present_genes)
inflam_genes_use   <- intersect(inflam_genes,   present_genes)

length(cytokine_genes_use)
length(inflam_genes_use)

features_list <- list(
  Cytokine      = cytokine_genes_use,
  Inflammatory  = inflam_genes_use
)

# UCell 默认会生成 meta.data 列名：Cytokine_UCell / Inflammatory_UCell
seu_re <- AddModuleScore_UCell(seu_re, features = features_list)

head(seu_re[[]] %>% dplyr::select(Cytokine_UCell, Inflammatory_UCell))

FeaturePlot(
  seu_re,
  features = "Cytokine_UCell",
  split.by = "group",
  cols = c("blue", "red")
)

FeaturePlot(
  seu_re,
  features = "Inflammatory_UCell",
  split.by = "group",
  cols = brewer.pal(n = 9, name = "RdBu")
)



