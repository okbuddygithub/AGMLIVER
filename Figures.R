#FRESH START

pathways.show = "TGFb"

FeaturePlot(AL, c("ADM","CALCRL"), label = T, label.size = 3)
ggsave(filename = paste0(pathways.show, "_ALGP_FP",".png"), last_plot(),dpi = 500 ,bg = "White")
netVisual_bubble(cellchatALGP, signaling = pathways.show, targets.use = c("1_Fib","3_HE","4_HPC"), remove.isolate = T)
netVisual_bubble(cellchatALGP, signaling = pathways.show, remove.isolate = T)
ggsave(filename = paste0(pathways.show, "_ALGP_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")

png(filename =  paste0(pathways.show, "_ALGL_HEAT",".png"), width = 5.04, height = 3.72, units = "in", res = 400)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", font.size = 10,)
dev.off()
netVisual_bubble(cellchatAGMKC, signaling = pathways.show, remove.isolate = T)



FeaturePlot(AL, c("ANGPT1","ANGPT2", "TEK", "ITGA5"), label = T, label.size = 2)

netVisual_bubble(cellchatAL, signaling = pathways.show, remove.isolate = T, targets.use = c("HE","HPC"))
netVisual_bubble(cellchatALGP, signaling = pathways.show, remove.isolate = T, targets.use = c("1_Fib","3_HE","4_HPC"))
netVisual_bubble(cellchatALC, signaling = pathways.show, remove.isolate = T, targets.use = c("HE","HPC"))


netVisual_bubble(cellchatAL, signaling = pathways.show, remove.isolate = T, targets.use = c("HE","HPC"))
netVisual_bubble(cellchatMK, signaling = pathways.show, remove.isolate = T, targets.use = c("HE","HPC"))

netVisual_bubble(cellchat, signaling = pathways.show, remove.isolate = T)
netVisual_bubble(cellchatMKGP, signaling = pathways.show, remove.isolate = T)
pathways.show = "SPP1"

FeaturePlot(AGM, c("ANGPTL1","ITGB1","ITGA1"), label = T, label.size = 2)
ggsave(filename = paste0(pathways.show, "_AGMGP_FP",".png"), last_plot(),dpi = 500 ,bg = "White")
netVisual_bubble(cellchatMKGP, signaling = pathways.show, targets.use = c("3_HE","4_HPC","2_Endo", "7_MK"), remove.isolate = T)
netVisual_bubble(cellchatMKGP, signaling = pathways.show, remove.isolate = T)
netVisual_bubble(cellchatMKG, signaling = pathways.show, remove.isolate = T)
ggsave(filename = paste0(pathways.show, "_MKGP_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")



png(filename =  paste0(pathways.show, "_AGMGP_HEAT",".png"), width = 5.04, height = 3.72, units = "in", res = 400)
netVisual_heatmap(cellchatMKGP, signaling = pathways.show, color.heatmap = "Reds", font.size = 10,)
dev.off()
netVisual_heatmap(cellchatMKG, signaling = pathways.show, color.heatmap = "Reds", font.size = 10,)
netVisual_bubble(cellchatALGP, signaling = pathways.show, remove.isolate = T)


                 

ggsave(filename = paste0(pathways.show, "_ALCtest_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_ALGLose_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")



# AGM Overtime 

p1<-DimPlot(AGMEHT)

cellsremove = "AGGTCATCAGTTTACG-1_4"
AGMEHT <- subset(AGMEHT, cells = cellsremove, invert = TRUE)
gene.show = "HLF"
FeaturePlot(AGMEHT, gene.show, split.by = "orig.ident", label = T,)
ggsave(filename = paste0(gene.show, "EHTTIME",".png"), last_plot(),dpi = 500 ,bg = "White")
VlnPlot(AGMEHT, gene.show, split.by = "orig.ident",)
ggsave(filename = paste0(gene.show, "VLNEHTTIME",".png"), last_plot(),dpi = 500 ,bg = "White")
FeaturePlot(AGM, gene.show, label = T,)
FeaturePlot(AGMEHT, "RAB27B", label = T,)
FeaturePlot(AGMEHT, gene.show)
FeaturePlot(AGMEHT, features = c("RUNX1",gene.show), blend = T, blend.threshold = .25, )
FeaturePlot(AGMEHT, features = c("RUNX1")) +NoLegend() +NoAxes()
#Corrolation 

FeaturePlot(srt, features = c("RUNX1","SNTB1"))

gg1 <- FeaturePlot(AGM, features = c("TNF")) +NoAxes() 
gg2 <- FeaturePlot(AGM, features = c("TNFRSF1A")) + NoLegend() + NoAxes() 
gg3 <- FeaturePlot(AGM, features = c("TNFRSF1B")) + NoLegend() + NoAxes() 
gg4 <- FeaturePlot(AGM, features = c("TPM1")) + NoLegend() + NoAxes() 


(gg3 + gg2) | (gg1)
(gg5 + gg3) | (gg4+ gg2) | (gg1+gg6)
gg1 | gg2

gg2+gg3 | gg1 + gg4


gg1 <- FeaturePlot(AGML, features = c("TNF"), label = T, label.size = 2) +NoAxes() 
gg2 <- FeaturePlot(AGML, features = c("WNT2B"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg3 <- FeaturePlot(AGML, features = c("FZD4"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg4 <- FeaturePlot(AGML, features = c("LRP6"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg5 <- FeaturePlot(AGML, features = c("LRP5"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg6 <- FeaturePlot(AGML, features = c("NOTCH4"), label = T, label.size = 1.5) + NoLegend() + NoAxes()


gg1 <- FeaturePlot(AL, features = c("ANGPT1"), label = T, label.size = 2) +NoAxes() 
gg2 <- FeaturePlot(AL, features = c("ANGPT2"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg3 <- FeaturePlot(AL, features = c("TEK"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg4 <- FeaturePlot(AL, features = c("ITGA5"), label = T, label.size = 2) + NoLegend() + NoAxes() 
#gg5 <- FeaturePlot(AL, features = c("F2RL2"), label = T, label.size = 2) + NoLegend() + NoAxes() 
gg2+gg3 | gg1+gg4


gg2 | gg1

gene = c("TGFB1",  "TGFB2",  "TGFBR1", "TGFBR2", "ACVR1")
VlnPlot(AGML, features = gene) 

FeaturePlot(AL, features = c("ANGPTL4")) 
FeaturePlot(AGML, features = c("TLR4")) 

FeaturePlot(AGM, features = c("WNT2B","FZD6")) + NoLegend() + NoAxes()
#DIfferences 
FeaturePlot(AL, features = c("PROC","PROCR"), combine = T, keep.scale = "feature")
FeaturePlot(AL, features = c("HSP90AA1","PROC","PROCR"), combine = T, keep.scale = "all")
FeaturePlot(AL, features = c("PROC","PROCR"), combine = T, keep.scale = "all")

netVisual_bubble(cellchatALGP, signaling = pathways.show, remove.isolate = T)
netVisual_bubble(cellchatAL, signaling = pathways.show, remove.isolate = T)

netVisual_bubble(cellchatAL, signaling = pathways.show, remove.isolate = T, targets.use = "HE")


plot <- FeaturePlot(object = AGM, features = c("MYB","RUNX1","WNT2B","IGFBP3")) + 
  NoLegend() + NoAxes()

plot[[1]] <- plot[[1]] + theme(legend.position="bottom") + 
  scale_x_continuous(breaks = seq(min(data$x), max(data$x), length=5)) +
  scale_y_continuous(breaks = seq(min(data$y), max(data$y), length=5))

library(CSCORE)
AGMEHT_B = AGMEHT[,AGMEHT$cell.type %in% '4_HPC']

mean_exp = rowMeans(AGMEHT_B@assays$RNA@counts/AGMEHT_B$nCount_RNA)
genes_selected = names(sort.int(mean_exp, decreasing = T))[1:1000]
genes_selected = "RUNX1"

AGMEHT_45 <- AGMEHT_B[, AGMEHT_B$orig.ident == "AGM_6week"]
CSCORE_result <- CSCORE(AGMEHT_B, genes = genes_selected)


CSCORE_coexp <- CSCORE_result$est

# Obtain BH-adjusted p values
CSCORE_p <- CSCORE_result$p_value
p_matrix_BH = matrix(0, length(genes_selected), length(genes_selected))
p_matrix_BH[upper.tri(p_matrix_BH)] = p.adjust(CSCORE_p[upper.tri(CSCORE_p)], method = "BH")
p_matrix_BH <- p_matrix_BH + t(p_matrix_BH)

# Set co-expression entires with BH-adjusted p-values greater than 0.05 to 0
CSCORE_coexp[p_matrix_BH > 0.05] <- 0


adj = WGCNA::adjacency.fromSimilarity(abs(CSCORE_coexp), power = 1)
# Compute the topological overlap matrix
TOM = WGCNA::TOMsimilarity(adj)
dissTOM = 1-TOM
rownames(dissTOM) <- colnames(dissTOM) <- genes_selected
# Run hierarchical clustering as in the WGCNA workflow
hclust_dist = hclust(as.dist(dissTOM), method = "average") 
memb = dynamicTreeCut::cutreeDynamic(dendro = hclust_dist, 
                                     distM = dissTOM, 
                                     deepSplit = 2,
                                     pamRespectsDendro = FALSE,
                                     minClusterSize = 10)
# For more instructions on how to tune the parameters in the WGCNA workflow,
# please refer to https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

names(memb) = genes_selected
memb_tab <- table(memb)
module_list = lapply(sort(unique(memb)), function(i_k) names(which(memb == i_k)))


#THE DOT PLOT

DotPlot(AGML, features =c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4"), group.by = lf)  + RotatedAxis() +  coord_flip()
DotPlot(AGMi, features =c("RUNX1","PTPRC","SPN","HLF","GFI1","MYB","MYCN","C1QA","CD14","LYVE1","LYZ","RNASE2")) + RotatedAxis() +  coord_flip()
DotPlot(AGM, features =c("COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","NKX2-3","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1")) + RotatedAxis() +  coord_flip()
DotPlot(AGML, features =c("EPCAM","AFP","FGB","APOA1","MAL","CALB1"))
DotPlot(AGML, features =c("HBE1","HBZ","GYPA"))
DotPlot(AGML, features =c("MKI67","TOP2A","AURKB")) +theme(
  # Rotate X axis text
  axis.text.x = element_text(angle = 90, hjust = 1),
  # Validate the theme
  validate = TRUE,)

SeuratAxes(RotatedAxis(90))

DotPlot(AGMi, features = c("CD3E","CCR7","IL7R","LTB","CD69","CD7"), idents = "Lymph")

DotPlot(AGM, features =c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4",
"RUNX1","PTPRC","SPN","HLF","GFI1","MYB","MYCN","C1QA","CD14","LYVE1","LYZ","RNASE2",
"COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","NKX2-3","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1",
"EPCAM","AFP","FGB","APOA1","MAL","CALB1",
"HBE1","HBZ","GYPA","MKI67","TOP2A","AURKB"),cols = c("lightgrey", "red")) +theme(
  # Rotate X axis text
  axis.text.x = element_text(angle = 90, hjust = 1),
  # Validate the theme
  validate = TRUE,)




DotPlot(AGM, group.by = "namedclusters",features =c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4",
                          "RUNX1","PTPRC","SPN","HLF","GFI1","MYB","MYCN","C1QA","CD14","LYVE1","LYZ","RNASE2",
                          "COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","NKX2-3","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1",
                          "EPCAM","AFP","FGB","APOA1","MAL","CALB1",
                          "HBE1","HBZ","GYPA","MKI67","TOP2A","AURKB"),cols = c("lightgrey", "red")) + RotatedAxis()+  FontSize(7)





unique_idents <- levels(Idents(AGML))
sort(unique_idents)




listo<- c("Endo1" ,  "Endo2"  , "Endo3" ,  "Endo4"  , "Endo5",  "HE", "HSC", "ERY",   "MK",  "Fibro1",  "Fibro2", 
 "Fibro3",  "Fibro4",  "Fibro5",  "Fibro6",  "Peri", "Kidney1", "Kidney2",
 "Kidney3", "Kidney4", "Liver",   "Lymph",   "Macro1",  "Macro2",  "Macro3",    
    "Mono",    "Muscle1", "Muscle2", "Mel")




DotPlot(AGMi, features =c("TNF","TNFRSF1A","TNFRSF1B")) + RotatedAxis()
DotPlot(AGMi, features =c("SPP1","ITGB1","ITGA5","ITGAV")) + RotatedAxis()
DotPlot(AGMi, features =c("NAMPT","ITGB1","ITGA5")) + RotatedAxis()
DotPlot(AGMi, features =c("TGFB1","TGFB2","TGFBR1","TGFBR2","ACVR1")) + RotatedAxis()
DotPlot(AGMi, features = c("NOTCH1","NOTCH3","NOTCH4","DLK1","JAG1","DLL3","DLL4")) + RotatedAxis()
DotPlot(AGM, features = c("WNT2B","FZD6","LRP6","FZD4","WNT9B")) + RotatedAxis()
DotPlot(AGMi, features =c("CXCL12","CXCR4")) + RotatedAxis()
DotPlot(AGMi, features =c("ANGPTL4","ANGPTL1","CDH5","ITGA1","ITGA5","ITGB1")) + RotatedAxis()

DotPlot(AGMi, features =c("TNF","TNFRSF1A","TNFRSF1B")) + RotatedAxis()
DotPlot(AGMi, features =c("SPP1","ITGB1","ITGA5","ITGAV")) + RotatedAxis()
DotPlot(AGMi, features =c("NAMPT","ITGB1","ITGA5")) + RotatedAxis()
DotPlot(AGMi, features =c("TGFB1","TGFB2","TGFBR1","TGFBR2","ACVR1")) + RotatedAxis()

VlnPlot(AGMi, features = c("TNFRSF1A","TNFRSF1B"),  pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5","HSC"))
VlnPlot(AGM, features = c("NOTCH1","NOTCH3","NOTCH4"),  pt.size = 0)
        
DotPlot(AGMi, features = c("WNT2B","FZD6","LRP6","FZD4","WNT9B")) + RotatedAxis()
DotPlot(AGMi, features =c("CXCL12","CXCR4")) + RotatedAxis()
DotPlot(AGMi, features =c("ANGPTL4","ANGPTL1","CDH5","ITGA1","ITGA5","ITGB1")) + RotatedAxis()

DimPlot(AL, label = T)

DotPlot(AL, features =c("CXCL12","CXCR4")) + RotatedAxis()
ggsave(filename = paste0("CXCL_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features = c("WNT2B","WNT2","FZD6","LRP6","FZD4")) + RotatedAxis()
ggsave(filename = paste0("WNT_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("ANGPTL4","ANGPTL2","ANGPTL1","CDH5","ITGA1","ITGA5","ITGB1")) + RotatedAxis()
ggsave(filename = paste0("ANGPTL_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("PROC","PROCR","F2R","HSP90AA1")) + RotatedAxis()
ggsave(filename = paste0("PROC_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("F2","F2R","PLG","PARD3")) + RotatedAxis()
ggsave(filename = paste0("PAR_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("BMP5","BMP7","BMPR2","ACVR1")) + RotatedAxis()
ggsave(filename = paste0("BMP_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("CALCRL","ADM")) + RotatedAxis()
ggsave(filename = paste0("CALCR_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")
DotPlot(AL, features =c("ANGPT1","ANGPT2","TEK","ITGA5")) + RotatedAxis()
ggsave(filename = paste0("ANGPT_ALG_DOTS",".png"), width = 5, height = 4, units = "in",last_plot(),dpi = 400 ,bg = "White")




Cell_identity_AGMMine = c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4",
                      "RUNX1","PTPRC","SPN",
                      "HBE1","HBZ","GYPA","ITGA2B","SELP","GP9",
                      "C1QA","CD14","LYVE1","LYZ","RNASE2","S100A8","FCGR3A",
                      "HLF","GFI1","MYB","MYCN",
                      "COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1",
                      "EPCAM","AFP","FGB","APOA1","ALB","MAL","CALB1","EDNRB")
DotPlot(AGM, features=Cell_identity_AGM, group.by = "fc" ,cols=c("grey90","red3"), cluster.idents = T, scale = T) +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1)) +theme(axis.text.x = element_text(size = 8)) + theme(axis.text.y = element_text(size = 6))+ NoLegend()

theme(axis.text.y = element_text(size = new_size))


HSC_surface_markers = c("PTPRC","SPN","ITGA2B","CDH5","CD34","THY1","VNN2","ITGAM","EMCN","PROCR","IL3RA","SELP","ACE","PROM1","HLA-DRA")
DotPlot(AGM, features=HSC_surface_markers, cols=c("grey90","red3")) +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1)) + NoLegend()

Hema_cell_identity_AGM = c("RUNX1","SPN","PTPRC",
                           "HOXA9","MLLT3","HLF","MYB","GATA2","MECOM",
                           "CD3E","CCR7","IL7R"," LTB","CD69","CD7",
                           "LYZ","S100A8","FCGR3A","MPO","AZU1","RNASE3",
                           "C1QA","C1QB","C1QC","LYVE1","RNASE1","CD14")
                          
DotPlot(AGM, features=Hema_cell_identity_AGM, cols=c("grey90","red3")) +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))

DotPlot(AL, features=Hema_cell_identity_AGM, cluster.idents = T) + NoLegend()+ RotatedAxis()

DotPlot(AGMi, features=Hema_cell_identity_AGMshort, cols=c("grey90","red3"), idents = c("Macro1","Macro2","Macro3","Mono","Lymph","HSC")) + coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))



pathways.show = "SPP1"
netVisual_heatmap(cellchat, signaling = pathways.show,  remove.isolate = T)
netVisual_bubble(cellchat, signaling = pathways.show,targets.use = "HE", remove.isolate = T,)  + theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(size = 8)) + NoLegend()

netVisual_bubble(cellchat, signaling = pathways.show, targets.use = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5"), remove.isolate = T)  + theme(axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(size = 6))

netVisual_bubble(cellchat, signaling = pathways.show,targets.use = "HE", sources.use = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5","Peri","Gran","Macro1","Macro2","Macro3","MK","HSC") ,remove.isolate = T)  + theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(size = 8)) + NoLegend()
netVisual_bubble(cellchat, signaling = pathways.show,targets.use = "HE", sources.use = c("Peri","Endo4") ,remove.isolate = T)  + theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(size = 8)) + NoLegend()
netVisual_bubble(cellchat, signaling = pathways.show,targets.use = "HE", sources.use = c("HE","Fibro6","Fibro4","Liver","Endo4") ,remove.isolate = T)  + NoLegend() + RotatedAxis() +theme(axis.text.y = element_text(size = 8)) + theme(axis.text.x = element_text(size = 6))
netVisual_bubble


MM<-FindMarkers(AGM, ident.1 = "Mel")


geneX= c("MMP2")
geneG= c("TGFBR1","TGFBR2","ACVR1")
geneG= c("FZD6","FZD4","LRP6")
geneG= c("TNFRSF1A","TNFRSF1B")
geneG= c("NOTCH1","NOTCH3","NOTCH4")
geneG= c("NOTCH1","NOTCH3","NOTCH4","DLK1","DLL4","DLL3","JAG1","JAG2")
geneG= c("ITGA5","ITGB1")     #"ITGA1","CDH5","CLDN5")
geneG= c("ITGA5","ITGAV","ITGA9","ITGB1")
geneG= c("F2R","F2RL2","PARD3")
geneG= c("PROCR","F2R")
geneG= c("PROCR","F2R","HSP90AA1")
geneG= c("CDH5","ITGA5")
geneG= c("ZEB2","RUNX1","SNAI1","TIE1")
geneG= c("TRADD","TRAF2","TRAF3","BCL2","CD40","CCL5","CCL13","BCL2A1","TNF","IRAK1","SLC31A1")
geneG =c("BIRC2","BIRC3", "NFKB1","NFKB2","RELA","REL","RELB")
geneG =c("SMAD2","SMAD3","MMP2")
geneG =c("PROC","F2","ANGPTL4")
geneG= c("ITGA5","ITGB1","ITGA4","ITGA9")
geneG= c("ITGA5","ITGB1","CLDN5","CDH5")
VlnPlot(AGM, geneX, group.by = "fg" ,pt.size = 0) + NoLegend()
VlnPlot(AGM, geneX, pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5")) + NoLegend()

VlnPlot(AGM, geneG, pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"),stack = T, flip = T)+ NoLegend()
VlnPlot(AGM, geneG, pt.size = 0, group.by = "fg", stack = T, flip = T)+ NoLegend()
VlnPlot(AGM, geneG, pt.size = 0, group.by = "fg", stack = T, flip = T, idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5"))+  NoLegend()

VlnPlot(AGM, geneG, pt.size = 0, idents = c("HE","Endo5"), split.by = "orig.ident")
VlnPlot(AGM, geneX, pt.size = 0, idents = c("HE"), split.by = "orig.ident")

DotPlot(AGM, geneG, group.by = "fg") + RotatedAxis() + NoLegend()
DotPlot(AGM, geneG, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5")) + RotatedAxis() 
DotPlot(AGM, geneG, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"), group.by = "fg") + RotatedAxis() +NoLegend()




VlnPlot(AGM, geneX, pt.size = 0)
VlnPlot(AGM, geneX, pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5","HSC"))


VlnPlot(AGM, geneG, pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"),stack = T, flip = T)+ stat_summary(fun.y = mean, geom='point', size = 15, colour = "black", shape = 95,) + NoLegend()

pathways.show = "SPP1"

netVisual_heatmap(cellchatAL, signaling = pathways.show,  remove.isolate = T)
netVisual_bubble(cellchatAL, signaling = pathways.show,targets.use = "HE", remove.isolate = T)  + theme(axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(size = 6)) 

netVisual_bubble(cellchatAL, signaling = pathways.show, targets.use = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5"), remove.isolate = T)  + theme(axis.text.y = element_text(size = 6)) + theme(axis.text.x = element_text(size = 6))


DotPlot(AL, geneG, group.by = "fg") + RotatedAxis() + NoLegend()
DotPlot(AL, geneG, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5","HPSC")) + RotatedAxis() 
DotPlot(AL, geneG, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"), group.by = "fg") + RotatedAxis() +NoLegend()
VlnPlot(AL, geneG, pt.size = 0, group.by = "fg", stack = T, flip = T)+ NoLegend()
VlnPlot(AL, geneG, pt.size = 0, group.by = "fc", stack = T, flip = T)+ NoLegend()


#Chondro
DotPlot(AL, features = c("FDCSP","SPP1","CLU","TM4SF4","ANXA4","S100A2","KRT17","KRT8","CD24","KRT19"))


#endo
DotPlot(AL, features = c("FCN3","DNASE1L3","CCL14","CRHBP","CLEC1B","NTS","RAMP2","S100A16","OIT3","ACP5"))+RotatedAxis()

#ery
DotPlot(AL, features = c("GYPA","HBG2","PRDX2","HEMGN","AHSP","ALAS2","BLVRB","HBG1","HMBS","MYL4"))

#fib
DotPlot(AL, features = c("CXCL14", "TPM1", "PTN", "MDK", "IGFBP3", "RBP1", "BEX3", "KRT18", "COLEC11", "CRABP2", "DCN", "COL1A1"))+ RotatedAxis()

#hepatoblast
DotPlot(AL, features = c("SPINK1", "AFP", "APOA2", "APOA1", "SERPINA1", "APOB", "TTR", "FABP1", "MT1E", "RBP4"))

#HEPTcyte
DotPlot(AL, features = c("ALB", "APOC3", "AHSG", "MT1G", "APOH", "VTN", "TTR", "MT1H", "ALDOB", "FABP1"))

#sin endo
DotPlot(AL, features = c("FGB", "STAB1", "SAA1", "SERPINA3", "ALB", "SAA2", "SERPINA1", "FGG", "STAB2", "HP")) + RotatedAxis()
#peri
DotPlot(AL, features = c("MUSTN1", "MYH11", "CRIP1", "ACTA2", "MYL9", "TAGLN", "TPM2", "SOD3", "ADIRF", "RERGL")) + RotatedAxis()
#fibro better?
DotPlot(AL, features = c("TAGLN", "COL1A1", "ACTA2", "COL3A1", "COL1A2", "RGS5", "MYL9", "BGN", "IGFBP7", "DCN")) + RotatedAxis()




shot<- SingleR(pca.sce, bp.ref, clusters= pca.sce$ident, labels=bp.ref$label.fine)

#ACTA2, COL1A1, TAGLN, COL1A2, COL3A1, RBP1

DotPlot(AL, "ELN")

#CONDENCED 

# Stel 
DotPlot(AL, features = c("DCN","BGN","COL1A1", "CCBE1", "COLEC10"))
#Afib
DotPlot(AL, features = c("ACTA2","TAGLN","MYL9","DES"))
#Endo
DotPlot(AL, features = c("CDH5","ESAM","KDR"))
#Hept
DotPlot(AL, features = c("ALB","TF","KRT18","TTR","HNF4A"))
#Fib
DotPlot(AL, features = c("LUM","FBLN1","FBLN2","PDGFRA"))



DotPlot(AL, features = c("DCN","BGN","COL1A1", "CCBE1", "COLEC10","ACTA2",
                         "TAGLN","MYL9","DES","CDH5","ESAM","KDR", "ALB","TF","KRT18","TTR",
                         "HNF4A","LUM","FBLN1","FBLN2","PDGFRA","SELP","HBZ","LY6H","LYVE1","CD14","C1QA",'LYZ',"RNASE2"), cluster.idents = T, group.by = "fct") + RotatedAxis()+ coord_flip()  + NoLegend() 

DotPlot(AL, features = c("DCN","BGN","COL1A1", "CCBE1", "COLEC10","ACTA2",
                         "TAGLN","MYL9","DES","CDH5","ESAM","KDR", "ALB","TF","KRT18","TTR",
                         "HNF4A","SELP","HBZ","LY6H","LYVE1","CD14","C1QA",'LYZ',"RNASE2","CEBPA","SPI1","FCGR3A","FCGR2","RUNX1","CD34"), cluster.idents = T, group.by = "fg") + RotatedAxis()+ coord_flip()  + NoLegend() 


DotPlot(AL, features = c("DCN","BGN","COL1A1", "CCBE1", "COLEC10","ACTA2",
                         "TAGLN","MYL9","DES","CDH5","ESAM","KDR", "ALB","TF","KRT18","TTR",
                         "HNF4A","SELP","HBZ","LYVE1","CD14","C1QA",'LYZ',"RNASE2"),cluster.idents = T, group.by = "fct") + RotatedAxis()+ coord_flip()  + NoLegend() 

#EHT


VlnPlot(AGM, features = c("SMAD2","SMAD3","RUNX1","JUN","MMP2"))



old_cluster_names <- "Fib1"
new_cluster_names <- "Stel4"
cellchatAL <- updateClusterLabels(cellchatAL, old.cluster.name = old_cluster_names, new.cluster.name = new_cluster_names)



vp_case1 <- function(gene_signature, file_name, test_sign){
  plot_case1 <- function(signature, y_max = NULL){
    VlnPlot(pbmc, features = signature,
            pt.size = 0.1, 
            group.by = "Response", 
            y.max = y_max # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  plot_list <- list()
  y_max_list <- list()
  for (gene in gene_signature) {
    plot_list[[gene]] <- plot_case1(gene)
    y_max_list[[gene]] <- max(plot_list[[gene]]$data[[gene]]) # get the max no. for each gene
    plot_list[[gene]] <- plot_case1(gene, y_max = (y_max_list[[gene]] + 1) )
  }
  cowplot::plot_grid(plotlist = plot_list)
  file_name <- paste0(file_name, "_r.png")
  ggsave(file_name, width = 14, height = 8)
}


gene_name = "ITGA5"
get_gene_value <- function(gene_name, data) {
  if (gene_name %in% rownames(data)) {
    return(data[gene_name, 1])
  } else {
    return("Gene not found")
  }
}
get_log_value <- function(gene_name, data) {
  if (gene_name %in% rownames(data)) {
    return(data[gene_name, 2])
  } else {
    return("Gene not found")
  }
}





get_gene_value("ITGA5", HEMark)
get_log_value("ITGAV", HEMark)

AL <- readRDS("~/ultliver.rds")
cellchatAL <- readRDS("~/LIVERCELLCHATALC.rds")
AGM <- readRDS("~/ultimateAGMf.rds")
AL <- readRDS("~/ultliver.rds")


AGM <- SetIdent(AGM, value = "fg")
HEMark <- FindMarkers(AGM, ident.1 = "HE",ident.2 = "Endo")


geneG= c("ITGA5","ITGB1")
geneG= c("CDH5","CLDN5")
geneG= ""

VlnPlot2(AGM, features = geneG, pt.alpha = F, stat.method = "wilcox.test", box = F, group.by = "fc",  comparisons = list(c(1,6), c(2,6), c(3,6), c(4,6), c(5,6)), hide.ns = FALSE) + stat_summary(fun.y = mean, geom='point', size = 15, colour = "black", shape = 95)
VlnPlot2(AGM, features = geneG, pt.alpha = F, stat.method = "wilcox.test", box = F, group.by = "fc",  comparisons = list(c(1,6), c(2,6), c(3,6), c(4,6), c(5,6)), hide.ns = FALSE)
VlnPlot(AGM, geneG, pt.size = 0, group.by = "fc", idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5"))
VlnPlot(AGM, geneG, pt.size = 0, group.by = "fc", idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5")) + NoLegend() + stat_summary(fun.y = mean, geom='point', size = 15, colour = "black", shape = 95)

VlnPlot2(AGM, features = geneG, pt.alpha = F, cells = cellsfg ,stat.method = "wilcox.test", box = F, group.by = "fg", hide.ns = FALSE)
VlnPlot2(AGM, features = geneG, pt.alpha = F, cells = cells ,stat.method = "wilcox.test", box = F, group.by = "fc",  comparisons = list(c(1,6), c(2,6), c(3,6), c(4,6)), hide.ns = FALSE)
VlnPlot2(AGM, features = geneG, pt.alpha = F, cells = cells ,stat.method = "wilcox.test", box = F, group.by = "fc", hide.ns = FALSE)
cells <- colnames(AGM)[AGM$fc %in% c("HE","Endo1","Endo2","Endo3","Endo4")]
cellsfg <- colnames(AGM)[AGM$fg %in% c("HE","Endo")]


geneX = "ITGB1"
VlnPlot(AGM, geneX, pt.size = 0, group.by = "fc", idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5")) + NoLegend() + stat_compare_means(label = "p.format", method = "t.test", comparisons = my_comparisons) + stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
VlnPlot(AGM, geneX, pt.size = 0, group.by = "fg", idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5")) + NoLegend()+ stat_compare_means()
VlnPlot(AGM, geneG, pt.size = 0, group.by = "fg", idents = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5")) + NoLegend()+ stat_compare_means() +  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
VlnPlot(AGM, geneX, pt.size = 0, group.by = "fc", idents = c("HE","Endo1")) + NoLegend()+ stat_compare_means( method = "t.test") +  stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)
#+ stat_summary(fun.y = median, geom='point', size = 15, colour = "black", shape = 95)



my_comparisonsname = list(c(1,6), c(2,6), c(3,6), c(4,6), c(5,6))
my_comparisonsname = list(c("HE","Endo1"), c("HE","Endo5"))


Cell_identity_AGM = c("CDH5","KDR","TIE1","APLNR","NRP2","GJA5","CXCR4",
                      "RUNX1","PTPRC","SPN",
                      "C1QA","CD14","LYVE1","LYZ","RNASE2",
                      "HLF","GFI1","MYB","MYCN",
                      "COL1A1","PDGFRA","CXCL12","POSTN","PTN","PAX1","SOX9","HAND1","NKX2−3","DCN","ALDH1A2","COL14A1","CRABP1","LUM","PAX3","TBX5","HAND2","REN","GATA3","ACTC1","ACTA2","NPHS2","DSC1","NPHS1",
                      "EPCAM","AFP","FGB","APOA1","MAL","CALB1",
                      "HBE1","HBZ","GYPA","ITGA2B","SELP","GP9")
DotPlot(sample, features=Cell_identity_AGM, cols=c("grey90","red3"), group.by="idents") +
  coord_flip() + theme(axis.text.x=element_text(angle=45, hjust=1))+ theme(axis.text.y = )






AGM <- SetIdent(AGM, value = "fc")
AGM <- RenameIdents(AGM, "Mel"= "Stromal3")


AGM@meta.data$fc <- AGM@active.ident

