# Here is a guide/documentation for the code and the steps I took to 
# make the seurat and cellchat objects as well as any of the graphs.
# Most of this is just modifications from the seurat and cellchat 
# guides on their website or github. Email for contact



#Combined Seurat Object 

#How to make a liver. 
# There were 3 liver objects from weeks 4.5, 5, and 6 respectively.

#The liver objects are from an older version of seurat and thus
#need to be updated 
#All objects are found on this google drive with the inventory 
# of what each object is found on the github link. Additionally the raw data is on 
#GEO. If the drive is down email me and I can send my own copies of the objects

#https://drive.google.com/drive/folders/1bsl4HMPh0ZZb9iAZTXD5lVY56sNj_MCm
#https://github.com/mikkolalab/Human-HSC-Ontogeny/blob/main/Inventory.csv


library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(Seurat)
library(dplyr)
library(data.table)
library(RColorBrewer)
library(SingleR)




load("/Users/CS/Documents/RObjects/EDFig2_5/seurat_object.Rdata")
LIV45 = UpdateSeuratObject(object = sample)

load("/Users/CS/Documents/RObjects/EDFig3_12/seurat_object.Rdata")
LIV5 = UpdateSeuratObject(object = sample)

load("/Users/CS/Documents/RObjects/EDFig3_13/seurat_object.Rdata")
LIV6 = UpdateSeuratObject(object = sample)



metadata <- LIV45@meta.data
metadata$orig.ident <- "LIV45"
LIV45@meta.data <- metadata

metadata <- LIV5@meta.data
metadata$orig.ident <- "LIV5"
LIV5@meta.data <- metadata

metadata <- LIV6@meta.data
metadata$orig.ident <- "LIV6"
LIV6@meta.data <- metadata


metadata <- AGML@meta.data
metadata$orig.ident <- "AGML"
AGML@meta.data <- metadata



#This can crash weaker computers so save before hand. 
LIV45 <- subset(LIV45, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)
LIV5 <- subset(LIV5, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)
LIV6 <- subset(LIV6, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)

LIVXAGM <- merge(LIV45, y = c(LIV5, LIV6, AGML), add.cell.ids = c("LIV45", "LIV5","LIV6", "AGML"))


#We have merged the 3 objects but they will be experiencing strong batch effects
# to negate these we will be integrating the objects. 

rm(LIV45,LIV5,LIV6)

liver.list <- SplitObject(LIVX, split.by = "orig.ident")
reference.list <- liver.list[c("LIV45","LIV5", "LIV6")]
species.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:20)
species.integrated <- IntegrateData(anchorset = species.anchors, dims = 1:20)
DefaultAssay(species.integrated) <- "integrated"
species.integrated <- ScaleData(species.integrated, verbose = FALSE)
species.integrated <- RunPCA(species.integrated, npcs = 20, verbose = FALSE)
species.integrated <- RunUMAP(species.integrated, dims = 1:20)
DimPlot(species.integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(species.integrated, reduction = "umap")



#The liver is all done. Now we need to get the AGM. 

#Insert Thom code for the creation of the AGM object 

#The AGM is certainly fascinating and all, but the total AGM contains 
# 21 thousand cells. Mix that with the liver object at a light 15 thousand cells 
# and you get an absurdly large matrix. One that neither R nor many computer 
#are particularity fond of. As such we have to be picky about which parts of the AGM we use. 
# We decided to use  the cells that make up the hegemonic endothelium and those that 
# are HPCs. The reasoning being that the rest of the tissue's 
# interaction with the liver, while potentially interesting are not a priority for the research.


AGM <- readRDS("/Users/CS/Documents/RObjects/AGM_harmony_merged_filt.rds")


Sample_subset <- subset(AGML, idents = c("HE","HPC","MK"))
AGML <- subset(AGM, idents = c("zero", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14","HPC","16","17","18","HE","ERY","21","MK","23","24","25","26","27","28"))

#if named already #Sample_subset <- subset(AGM, idents = c("HE", "HPC", "MK"))


#We now have a subset of the AGM, and it is time to mix it in with the liver



#I would recommend removing seurat objects that you no longer need.

#rm(LIV45,LIV5,LIV6,AGM,liver.list,species.anchors,sample)

seuratobj1 <- Sample_subset
seuratobj2 <- LIVX

#metadata <- seuratobj1@meta.data
#metadata$orig.ident <- "AGM"
#seuratobj1@meta.data <- metadata

#metadata <- seuratobj2@meta.data
#metadata$orig.ident <- "LIV"
#seuratobj2@meta.data <- metadata



seuratobj1[["percent.mt"]] <- PercentageFeatureSet(seuratobj1, pattern = "^MT-")
AGM500 <- subset(AGM, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 10)
seuratobj1 <- NormalizeData(seuratobj1, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj1 <- FindVariableFeatures(seuratobj1, selection.method = "vst", nfeatures = 2000)


seuratobj2 <- subset(seuratobj2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mito < .1)
seuratobj2 <- NormalizeData(seuratobj2, normalization.method = "LogNormalize", scale.factor = 10000)
seuratobj2 <- FindVariableFeatures(seuratobj2, selection.method = "vst", nfeatures = 2000)



#You also have the ability to combine more than 2 objects. This can be done
#by putting c(x,y,z,...) after y = and also add in the IDs you need to add. 
comboseurat <- merge(seuratobj1, y = seuratobj2, add.cell.ids = c("AGM", "LIV45","LIV5","LIV6"))
table(comboseurat@meta.data$orig.ident)
comboseurat <- FindVariableFeatures(comboseurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(comboseurat)
comboseurat <- ScaleData(comboseurat, features = all.genes)
comboseurat <- RunPCA(comboseurat, features = VariableFeatures(object = comboseurat)) 
#ElbowPlot(pca)
# KEEP THESE SETTINGS 
comboseurat <- FindNeighbors(comboseurat, dims = 1:20)
comboseurat <- FindClusters(comboseurat, resolution = .85)
comboseurat <- RunUMAP(comboseurat, dims = 1:20)
DimPlot(comboseurat, reduction = "umap", label = T)
DimPlot(comboseurat, group.by = "orig.ident")

# We now have our final seurat object. Now its time for SingleR and CellChat


#SingleR to identify groups of cells. SingleR is not compatible with 
#an intergrated object mixed with a non intergrated object. As such we will 
#remove the intergration assay. It results in no changes to the dimplot or cells 
# it is just a holdover from prior seurat objects. 

AL@assays$integrated <- NULL

bp.ref <- BlueprintEncodeData()
hp.ref <- HumanPrimaryCellAtlasData()
bpalt.ref <- celldex::BlueprintEncodeData()


options(ggrepel.max.overlaps = 100)
AL.sce <- as.SingleCellExperiment(AL)
bpref.main <- SingleR(test = AL.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.main)
bpref.ont <- SingleR(test = pca.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.ont)
bpref.fine <- SingleR(test = pca.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.fine)
hpref.main <- SingleR(test = AL.sce,assay.type.test = 1,ref = hp.ref,labels = hp.ref$label.main)
#hpref.fine <- SingleR(test = pca.sce,assay.type.test = 1,ref = hp.ref,labels = hp.ref$label.fine)
comboseurat@meta.data$bp.main   <- bpref.main$labels
comboseurat@meta.data$bp.fine   <- bpref.fine$labels
comboseurat@meta.data$hp.main   <- hpref.main$labels
comboseurat@meta.data$hp.fine   <- hpref.fine$labels
DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "bp.fine") + NoLegend()
DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "hp.main") + NoLegend()


pca.sce <- as.SingleCellExperiment(AGM)
bpref.main <- SingleR(test = pca.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.main)
bpref.fine <- SingleR(test = pca.sce,assay.type.test = 1,ref = bp.ref,labels = bp.ref$label.fine)
hpref.main <- SingleR(test = pca.sce,assay.type.test = 1,ref = hp.ref,labels = hp.ref$label.main)
#AGM@meta.data$bp.main   <- bpref.main$labels
AGM@meta.data$bp.fine   <- bpref.fine$labels
AGM@meta.data$hp.main   <- hpref.main$labels
DimPlot(LIVX, label = T, repel= T, label.size = 2.5, group.by = "hp.main") +NoLegend()
DimPlot(LIVX, label = T, repel= T, label.size = 3, group.by = "bp.main") + NoLegend()
DimPlot(AGM, label = T, repel= T, label.size = 2.5, group.by = "bp.fine") 
#More modular version 
seurat = AGM
reftype = "bp.ont"
DimPlot(AGM, label = T, repel= T, label.size = 3, group.by = reftype) + NoLegend()

AGM <-  SetIdent(AGM, value = "cell.type")
#SingleR graphs can be quite granular so using hoverlocator can make it clear which 
#clusters are which cell types. Additionally, the gene labels will indicate original 
# identity within AGM or Liver. for the AGM cells 1_1 1_2 and so on correspond to 
# week 4.5,5, 5.5, and 6. You can switch between grouping by cluster or ident 
# and showing cluster or ident on the hover locator. You can also show multiple 
#variables by using the code for gg4. There are 2 main annotation databases used for this 
#being the Human primary cell atlas and the blueprint data. The former is better 
#for more general annotations while bp was superior for immune cells. We used both 
# to annotate the cells. As with any cell annotation there is some ambiguity 
# and should be taken with a slight grain of salt. Also the fact that this is 
# embryonic tissue complicates things. 







gg1 <- DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg1, information = FetchData(comboseurat,vars = "seurat_clusters"))

gg2 <- DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "bp.main") + NoLegend()
HoverLocator(gg2, information = FetchData(comboseurat,vars = "bp.main"))


gg3 <- DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg3, information = FetchData(comboseurat,vars = "bp.main"))

gg4 <- DimPlot(comboseurat, label = T, repel= T, label.size = 3, group.by = "seurat_clusters") + NoLegend()
HoverLocator(gg4, information = FetchData(comboseurat,vars = c("bp.main", "seurat_clusters")))

groupings <- "bp.fine"
vars = c("fg", "bp.fine")
gg5 <- DimPlot(AL, label = T, repel= T, label.size = 3, group.by = groupings) + NoLegend()
HoverLocator(gg5, information = FetchData(AL,vars = vars))


hpca2 <- hp.ref
hpca2$label.main <- paste0("HPCA.", hpca2$label.main)

bpe2 <- bp.ref
bpe2$label.fine <- paste0("BPE.", bpe2$label.fine)

shared <- intersect(rownames(hpca2), rownames(bpe2))
combined <- cbind(hpca2[shared,], bpe2[shared,])


plotScoreHeatmap(bpref.fine, clusters = pca.sce$ident)
#plotScoreDistribution(bpref.fine)
plotDeltaDistribution(shotgroup)

shot<- SingleR(pca.sce, bp.ref, clusters= pca.sce$ident, labels=bp.ref$label.main)
shotf<- SingleR(pca.sce, bp.ref, clusters= pca.sce$ident, labels=bp.ref$label.fine)
shotfg<- SingleR(pca.sce, bp.ref, clusters= pca.sce$fg, labels=bp.ref$label.main)
shotffg<- SingleR(pca.sce, bp.ref, clusters= pca.sce$fg, labels=bp.ref$label.fine)


livbulk<- SingleR(AL.sce, hp.ref, clusters= AL.sce$ident, labels=hp.ref$label.main)
livclus<- SingleR(AL.sce, hp.ref, clusters= AL.sce$fc, labels=hp.ref$label.main)

shots<- SingleR(pca.sce, combined, clusters= pca.sce$ident, labels=c(combined$label.main))

shothpcafine<- SingleR(pca.sce, hp.ref, clusters= pca.sce$fc, labels=hp.ref$label.fine)

tab <- table(clusters= levels(pca.sce$ident), label=shot$labels) 
tab <- table(clusters= levels(pca.sce$ident), label=shotf$labels) 
tab <- table(clusters= levels(pca.sce$fg), label=shotfg$labels) 
tab <- table(clusters= levels(pca.sce$fg), label=shotffg$labels) 
tab <- table(clusters= levels(pca.sce$fc), label=shothpca$labels) 
tab <- table(clusters= levels(pca.sce$fc), label=shothpcafine$labels)



tab <- table(clusters= levels(AL.sce$ident), label=livbulk$labels)

pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 


tab <- table(clusters= levels(pca.sce$fct), label=shotgroup$labels) 
pheatmap::pheatmap(log10(tab+10)) # using a larger pseudo-count for smoothing. 

tab <- table(clusters= pca.sce$ident, label=hpref.main$labels) 

tab <- table(clusters= pca.sce$ident, label=bpref.fine$labels) 
tab <- table(clusters= pca.sce$fg, label=bpref.fine$labels) 
pheatmap::pheatmap(log10(tab+10), treeheight_row = 0, treeheight_col = 0) # using a larger pseudo-count for smoothing. 

# Once identified we need to group the cells and prepare for cellchat. 
# while one could simply memorize all the clusters numbers by staring at dimplots 
# for hours on end, it is a bit more practical to label the clusters. 
# Furthermore we will be establishing overall cell classes based on the 
# single R data. 
# There are 8 Classes: HE, HPC, MK, ERY ,Fib, Endo, Hept, and Hemo
# Short for Hemogenic endothelium, Hemopoetic stem cells, Megakaryocytes
# Fibroblasts, Endothelila cells, Hepatocytes, and Blood cells in the liver
#If you are attempting to use this code for your own project, name your object prior 
# to making it a cellchat. You can alter the names, but it is a hassle.

comboseurat <- SetIdent(comboseurat, value = "namedclusters")

comboseurat <- RenameIdents(object = comboseurat, "0" = "2_Endo", "1" = "1_Fib", "2" = "5_Hem", "3" = "5_Hem", "4" = "5_Hem", "5" = "1_Fib", "6" = "2_Endo", "7" = "5_Hem", "8" = "1_Fib", "9" = "1_Fib", "10" = "5_Hem", "11" = "5_Hem", "12" = "5_Hem", "13" = "5_Hem", "14" = "2_Endo", "15" = "8_Hept", "16" = "7_MK", "17" = "3_HE", "18" = "1_Fib", "19" = "4_HPC", "20" = "2_Endo", "21" = "6_Ery", "22" = "8_Hept", "23" = "2_Endo", "24" = "1_Fib", "25" = "5_Hem")
Idents(comboseurat) <- factor(Idents(comboseurat), levels=rev(c("1_Fib", "2_Endo", "3_HE", "4_HPC", "5_Hem","6_Ery","7_MK", "8_Hept")))
comboseurat$'cell.type' <- Idents(comboseurat) 
DimPlot(comboseurat, label = T)  


AL6 <- RenameIdents(object = AL6, "0" = "Endo", "1" = "Fib", "2" = "Hem", "3" = "Hem", "4" = "Hem", "5" = "Fib", "6" = "Endo", "7" = "Hem", "8" = "Fib", "9" = "Fib", "10" = "Hem", "11" = "Hem", "12" = "Hem", "13" = "Hem", "14" = "Endo", "15" = "Hept", "16" = "MK", "17" = "HE", "18" = "Fib", "19" = "HPC", "20" = "Endo", "21" = "Ery", "22" = "Hept", "23" = "Endo", "24" = "Fib", "25" = "Hem")
Idents(AL6) <- factor(Idents(AL6), levels=rev(c("Fib", "Endo", "HE", "HPC", "Hem","Ery","MK", "Hept")))
AL6$'celltype' <- Idents(AL6) 
DimPlot(AL6, label = T)  

data <- subset(x = AL, downsample = 200)


AL <- SetIdent(AL, value = "seurat_clusters")


# The "zero" is essential as cellchat hates the number 0
new.cluster.ids <-  c("zero", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14","Hept1","MK","HE","18","HPC","20","ERY","Hept2","23","24","25")
names(new.cluster.ids) <- levels(comboseurat)
comboseurat <- RenameIdents(comboseurat, new.cluster.ids)
comboseurat$'namedclusters' <- Idents(comboseurat)
DimPlot(comboseurat, label = T)


DimPlot(comboseurat, group.by = "cell.type", label = T)
DimPlot(comboseurat, group.by = "namedclusters", label = T)

#AGM

AGM <- readRDS("/Users/CS/Documents/RObjects/AGM_harmony_merged_filt.rds")
AGM <- FindClusters(AGM, resolution = .85)


new.cluster.ids <-  c("zero", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14","HPC","16","17","18","HE","ERY","21","MK","23","24","25","26","27","28","29")
names(new.cluster.ids) <- levels(AGM)
AGM <- RenameIdents(AGM, new.cluster.ids)
AGM$'namedclusters' <- Idents(AGM)
DimPlot(AGM, label = T)


AGM <- RenameIdents(object = AGM, "zero"="1_Stromal", "1"="2_Endo", "2"="2_Endo", "3"="2_Endo", "4"="1_Stromal", "5"="1_Stromal", "6"="5_Hem", "7"="2_Endo", "8"="5_Hem", "9"="1_Stromal", "10"="2_Endo", "11"="1_Stromal", "12"="1_Stromal", "13"="5_Hem", "14"="5_Hem", "HPC"="4_HPC", "16"="1_Stromal", "17"="1_Stromal", "18"="5_Hem", "HE"="3_HE", "ERY"="6_Ery", "21"="1_Stromal", "MK"="7_MK", "23"="1_Stromal", "24"="1_Stromal", "25"="1_Stromal", "26"="1_Stromal", "27"="1_Stromal", "28"="1_Stromal")
Idents(AGM) <- factor(Idents(AGM), levels=rev(c("1_Stromal", "2_Endo", "3_HE", "4_HPC", "5_Hem","6_Ery","7_MK")))
AGM$'cell.type' <- Idents(AGM) 
DimPlot(AGM, label = T)

AGM <- RenameIdents(object = AGM, "zero"="1_Stromal", "1"="2_Endo", "2"="2_Endo", "3"="2_Endo", "4"="1_Stromal", "5"="1_Stromal", "6"="5_Hem", "7"="2_Endo", "8"="5_Hem", "9"="1_Stromal", "10"="2_Endo", "11"="8_Peri", "12"="1_Stromal", "13"="5_Hem", "14"="5_Hem", "HPC"="4_HPC", "16"="1_Stromal", "17"="1_Stromal", "18"="5_Hem", "HE"="3_HE", "ERY"="6_Ery", "21"="1_Stromal", "MK"="7_MK", "23"="1_Stromal", "24"="1_Stromal", "25"="1_Stromal", "26"="1_Stromal", "27"="1_Stromal", "28"="1_Stromal")
Idents(AGM) <- factor(Idents(AGM), levels=rev(c("1_Stromal", "2_Endo", "3_HE", "4_HPC", "5_Hem","6_Ery","7_MK","8_Peri")))
AGM$'cell.typep' <- Idents(AGM) 
DimPlot(AGM, label = T)

AGM <- RenameIdents(object = AGM, "zero"="1_Stromal", "1"="Endo", "2"="Endo", "3"="Endo", "4"="", "5"="1_Stromal", "6"="5_Hem", "7"="2_Endo", "8"="5_Hem", "9"="1_Stromal", "10"="2_Endo", "11"="8_Peri", "12"="1_Stromal", "13"="5_Hem", "14"="5_Hem", "HPC"="4_HPC", "16"="1_Stromal", "17"="1_Stromal", "18"="5_Hem", "HE"="3_HE", "ERY"="6_Ery", "21"="1_Stromal", "MK"="7_MK", "23"="1_Stromal", "24"="1_Stromal", "25"="1_Stromal", "26"="1_Stromal", "27"="1_Stromal", "28"="1_Stromal")
Idents(AGM) <- factor(Idents(AGM), levels=rev(c("1_Stromal", "2_Endo", "3_HE", "4_HPC", "5_Hem","6_Ery","7_MK","8_Peri")))
AGM$'cell.typep' <- Idents(AGM) 
DimPlot(AGM, label = T)

AGM <- SetIdent(AGM, value = "namedclusters")
AGM <- SetIdent(AGM, value = "cell.type")
DimPlot(AGM, label = T)
#Time to chat about CellChat
# Cellchat is our tool for examining the ligand receptor interactions 
#between the liver and the AGM. It generally runs the same each time, but you 
#can make changes to a few parameters. Database and group.by do have some variability
# the entire cellchat database is quite large and if it needs to look at every type 
# of interaction it may find irrelevent interactions, waste computing power/time, 
# and obsure relevent interactions. If you are looking for long distance signalling 
# then cell cell contact or ECM interactions arent going to be relevent. 
# you can set the program to look at all interactions or only a subset. 
# we choose only secreted signalling. 

# Additionally, you can group cells by any meta label in the object. 
# you can group them by cell type or cluster. We choose to do both and made 2 objects

cellchat <- createCellChat(object = AGM, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human
###########################################################################
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
#If you want it all use the code before instead.
CellChatDB.use <- CellChatDB
#########################################################################
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("CXCL"))
#pathways.show <- c("CXCL")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, "AGMPSSG.rds")

rm(cellchat)
computeCommunProb()
AGM <- SetIdent(AGM, value = "namedclusters")

cellchat <- createCellChat(object = AGM, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.human
###########################################################################
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation") 
#If you want it all use the code before instead.
#CellChatDB.use <- CellChatDB
#########################################################################
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
#df.net <- subsetCommunication(cellchat, signaling = c("CXCL"))
#pathways.show <- c("CXCL")
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)
saveRDS(cellchat, "AGMSSC.rds")


#If graphing isnt your style then dataframes are the way to go with cellchat.
#You can specefy the targets and the recpients. If you are using a named cluster 
#such as the HE, then you need to set it up such that targets.use = c("HE", any other targets you want)
dfaaaa.All <-subsetCommunication(cellchat, 
                             sources.use = c(0-35), targets.use = c(0-35))
df.AGM <-subsetCommunication(cellchatMK, 
                             sources.use = c(0-30), targets.use = c(0-30))
dfAGMk.all <-subsetCommunication(cellchatMKC, 
                              sources.use = c(0-35), targets.use = c(0-35))

df.AGMHE <- subsetCommunication(cellchatMKG, slot.name = "netP", targets.use = "3_HE")
dfg.test<- subsetCommunication(cellchatALGP, slot.name = "netP", targets.use = "3_HE")


df.All <-subsetCommunication(cellchatAGM45G, 
                                 sources.use = c(0-35), targets.use = c(0-35))

pathways.show <- "TGFb"
#This shows overall signalling by pathway. This includes every LR interaction 
#under a single pathway. 

group.cellTypeAL <- c("2_Endo", "1_Fib", "5_Hem", "5_Hem", "5_Hem", "1_Fib", "2_Endo", "5_Hem", "1_Fib", "1_Fib", "5_Hem", "5_Hem", "5_Hem", "5_Hem", "2_Endo", "8_Hept", "7_MK", "3_HE", "1_Fib", "4_HPC", "2_Endo", "6_Ery", "8_Hept","2_Endo","1_Fib", "5_Hem") 
group.cellTypeAL <- factor(group.cellTypeAL, levels = unique(group.cellTypeAL))
names(group.cellTypeAL) <- levels(cellchatAL@idents)


netVisual_chord_cell(cellchatALC, signaling = pathways.show, group = group.cellTypeAL, title.name = paste0(pathways.show, " signaling network"))

#AGM version
group.cellTypeMK <- c("1_Stromal", "2_Endo", "2_Endo", "2_Endo", "1_Stromal", "1_Stromal", "5_Hem", "2_Endo", "5_Hem", "1_Stromal", "2_Endo", "1_Stromal", "1_Stromal", "5_Hem", "5_Hem", "4_HPC", "1_Stromal", "1_Stromal", "5_Hem", "3_HE", "6_Ery", "1_Stromal", "7_MK", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal")
group.cellTypeMK <- factor(group.cellTypeMK, levels = unique(AGM$cell.type))
names(group.cellTypeMK) <- levels(cellchatMKC@idents)
netVisual_chord_cell(cellchatMKC, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))

group.cellTypeMKP <- c("1_Stromal", "2_Endo", "2_Endo", "2_Endo", "1_Stromal", "1_Stromal", "5_Hem", "2_Endo", "5_Hem", "1_Stromal", "2_Endo", "8_Peri", "1_Stromal", "5_Hem", "5_Hem", "4_HPC", "1_Stromal", "1_Stromal", "5_Hem", "3_HE", "6_Ery", "1_Stromal", "7_MK", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal", "1_Stromal")
group.cellTypeMKP <- factor(group.cellTypeMK, levels = unique(AGM$cell.typep))
names(group.cellTypeMKP) <- levels(cellchat@idents)
netVisual_chord_cell(cellchatMKC, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))

group.cellTypeMK <- c("9_Epi", "2_Endo", "2_Endo", "2_Endo", "1_Stromal", "1_Stromal", "5_Hem", "2_Endo", "5_Hem", "1_Stromal", "2_Endo", "1_Stromal", "8_Peri", "5_Hem", "5_Hem", "4_HPC", "1_Stromal", "9_Epi", "5_Hem", "3_HE", "6_Ery", "9_Epi", "7_MK", "1_Stromal", "9_Epi", "1_Stromal", "9_Epi", "1_Stromal", "9_Epi")
group.cellTypeMK <- factor(group.cellTypeMK, levels = unique(AGM6$cell.typestroma))
names(group.cellTypeMK) <- levels(cellchatAGM5@idents)
netVisual_chord_cell(cellchatAGM5, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))

group.cellTypeMK <- c("9_Epi", "2_Endo", "2_Endo", "2_Endo", "1_Stromal", "1_Stromal", "5_Hem", "2_Endo", "5_Hem", "1_Stromal", "2_Endo", "1_Stromal", "8_Peri", "5_Hem", "5_Hem", "4_HPC", "1_Stromal", "9_Epi", "5_Hem", "3_HE", "6_Ery", "9_Epi", "7_MK", "1_Stromal", "9_Epi", "1_Stromal", "9_Epi", "1_Stromal")
group.cellTypeMK <- factor(group.cellTypeMK, levels = unique(FAGML$id))
names(group.cellTypeMK) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))

netVisual_chord_cell(cellchatAGM45, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))
netVisual_chord_cell(cellchatAGM5, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))
netVisual_chord_cell(cellchatAGM6, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))


table(AGM6$cell.typestroma)

pathways.show <- "ANGPTL"
netVisual_chord_cell(cellchat, signaling = pathways.show,)
#If you only want one specefic LR interaction then run the first line and find 
# the LR pair in pairLR and imput it into the second line. 
#In this case we only see CXCR4 and CXCL12 
netAnalysis_contribution(cellchat, signaling = pathways.show, targets.use = "HE")
netVisual_bubble(cellchat, signaling = pathways.show)
netVisual_bubble(cellchat, signaling = pathways.show, targets.use = "HE")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_chord_cell(cellchatAL, signaling = pathways.show, group = group.cellTypeAL, title.name = paste0(pathways.show, " signaling network"))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(cellchat, pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE)
LR.show <- pairLR[1,] # show one ligand-receptor pair
netVisual_individual(cellchat, pathways.show, pairLR.use = LR.show ,layout = "chord")


#Feature Plots. Allows you to see the expression of the genes involved in the 
#interaction. Usually just 2 genes but there can be additional genes.
genes <- strsplit(LR.show, "_")[[1]]
FeaturePlot(AGM, features = c(genes[1], genes[2]))
FeaturePlot(AGM, features = c(genes[1], genes[2], genes[3], "ACVR1"), label = F)


#Overall Content

netAnalysis_signalingRole_network(cellchatALG, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netAnalysis_signalingRole_heatmap(cellchatALC, pattern = "incoming", width = 18, height = 20, font.size = 5, title = "AGM/LIVER Incoming Signals")
netAnalysis_signalingRole_heatmap(cellchatALG, pattern = "incoming", width = 18, height = 20, font.size = 5, title = "AGM/LIVER Incoming Signals")
netAnalysis_signalingRole_heatmap(cellchatALG, pattern = "outgoing", width = 18, height = 20, font.size = 5, title = "AGM/LIVER  Outgoing Signals")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "all", width = 18, height = 30, font.size = 5, title = "AGM/LIVER Groups")

netVisual_heatmap(cellchatALG, measure = "count", slot.name = "netP", color.use = c("Red","Blue","Green","Yellow","Grey","Purple","Orange","Pink"))
netVisual_heatmap(cellchatALG, measure = "weight", slot.name = "netP", color.use = c("Red","Blue","Green","Yellow","Grey","Purple","Orange","Pink"))

netVisual_heatmap(cellchat, measure = "weight", slot.name = "netP", color.use = c("Red","Blue","Green","Yellow","Grey","Purple","Orange","Pink"))
netVisual_heatmap(cellchat, measure = "weight", slot.name = "net", color.heatmap = "Reds")


netVisual_heatmap(cellchat, measure = "weight", slot.name = "netP", color.heatmap = "Reds")
netVisual_heatmap(cellchatMKG, measure = "count", slot.name = "netP", color.heatmap = "Reds")

netVisual_heatmap(cellchatMKC, measure = "weight", slot.name = "netP", color.heatmap = "Reds")
netVisual_heatmap(cellchatALG, measure = "weight", slot.name = "netP", color.heatmap = "Reds")

plotGeneExpression(cellchatMKG, signaling = pathways.show)

netVisual_heatmap(cellchatMKAGMf, measure = "weight", slot.name = "net", color.heatmap = "Reds")


CellChat::computeAveExpr(cellchat.agm, "MMP2")
AverageExpression(AGM, features = "CXCL12")
AggregateExpression(AGM, features = "CXCL12")
CellChat::computeAveExpr(cellchat, "CXCL12")
AverageExpression(comboseurat, features = "CXCL12")
AggregateExpression(comboseurat, features = "CXCL12")

#Ok that is the basic workflow for the graphs 
gene = "MMP2"
FeaturePlot(comboseurat, gene ) + FeaturePlot(AGM, gene)

FeaturePlot(AGM, "RUNX1",label = T)
FeaturePlot(AGM, gene, label = T)


#.agm has all pathways and MK cluster 
#.al is only secreted signalling of the AL object/comboseurat




pathways.show <- "CXCL"


#netAnalysis_contribution(cellchat, signaling = pathways.show)
netVisual_bubble(cellchatMKC, signaling = pathways.show, targets.use = "HE")
netVisual_bubble(cellchatMKC, signaling = pathways.show, targets.use = )
ggsave(filename = paste0(pathways.show, "_AGMHE_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")
"AGM Clusters"
netAnalysis_signalingRole_network(cellchatMKC, signaling = pathways.show, width = 8, height = 4, font.size = 8)
#ggsave(filename = paste0(pathways.show, "_AGMC_ROLE",".png"), last_plot(),dpi = 400 ,bg = "White")
netVisual_heatmap(cellchatMKC, signaling = pathways.show, color.heatmap = "Reds")
#ggsave(filename = paste0(pathways.show, "_AGMC_HEAT",".png"), last_plot(),dpi = 400 ,bg = "White")
netAnalysis_signalingRole_scatter(cellchatMKC, pathways.show)
ggsave(filename = paste0(pathways.show, "_AGMC_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
netVisual_chord_cell(cellchatMKC, signaling = pathways.show, group = group.cellTypeMK, title.name = paste0(pathways.show, " signaling network"))
ggsave(filename = paste0(pathways.show, "_AGMC_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")

pairLR <- extractEnrichedLR(cellchatMKC, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[22,] # show one ligand-receptor pair
netVisual_individual(cellchatMKC, pathways.show, pairLR.use = LR.show, group = group.cellTypeMK ,layout = "chord")

genes <- strsplit(LR.show, "_")[[1]]
FeaturePlot(AGM, features = c(genes[1], genes[2], genes[3], "FZD4"), label = T, label.size = 2)
ggsave(filename = paste0(pathways.show, "_AGMC_FP",".png"), last_plot(),dpi = 500 ,bg = "White")

FeaturePlot(AGM, c("WNT5B","WNT5A"), label = T, label.size = 2)  

"AGM Celltype"
netVisual_bubble(cellchatMKG, signaling = pathways.show, targets.use = "3_HE")
netAnalysis_signalingRole_network(cellchatMKG, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(cellchatMKG, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(cellchatMKG, pathways.show)
netVisual_chord_cell(cellchatMKG, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))

netVisual_bubble(cellchatMKG, signaling = pathways.show)
netAnalysis_contribution(cellchatMKG, signaling = pathways.show)

pairLR <- extractEnrichedLR(cellchatMKG, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[2,] # show one ligand-receptor pair
netVisual_individual(cellchatMKG, pathways.show, pairLR.use = LR.show ,layout = "chord")

pairLRGT <- extractEnrichedLR(cellchatALG, signaling = pathways.show, geneLR.return = TRUE)
DotPlot(AL, features = pairLRGT$geneLR) + RotatedAxis()


#plotGeneExpression(cellchatmkfull, signaling = pathways.show)
plotGeneExpression(cellchat, signaling = pathways.show) + NoLegend()
ggsave(filename = paste0(pathways.show, "_AGMG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")

ggsave(filename = paste0(pathways.show, "_AGMG_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_ROLE",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_HEAT",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_FP",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMG_DOTSEX",".png"),width = 5, height = 3.4, units = "in" ,last_plot(),dpi = 400 ,bg = "White") 


ggsave()


netAnalysis_signalingRole_network(cellchatmkfull, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(cellchatmkfull, signaling = pathways.show, color.heatmap = "Reds", sources.use = c("zero", "1", "2","3","4","5","6","7","8","9","10","11","13","14","HPC","16","17","18","HE","ERY","21","MK","24","25","26","27","28","29"))
netAnalysis_signalingRole_scatter(cellchatmkfull, pathways.show)
netVisual_chord_cell(cellchatmkfull, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))

pathways.show <- "TGFb"

"AGMP GROUP"
netVisual_heatmap(AGMPSSG, signaling = pathways.show, color.heatmap = "Reds")
netVisual_bubble(AGMPSSG, signaling = pathways.show, targets.use = "3_HE")
ggsave(filename = paste0(pathways.show, "_AGMPGHE_BUBBLE",".png"), width = 4, height = 4, units = "in",last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPGHE_BUBBLE",".png"),last_plot(),dpi = 500 ,bg = "White")
plotGeneExpression(AGMPSSG, signaling = pathways.show)
ggsave(filename = paste0(pathways.show, "_AGMPG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")
pairLRGT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE)
DotPlot(AGM, features = pairLRGT$geneLR) + RotatedAxis()
ggsave(filename = paste0(pathways.show, "_AGMPG_DOTSEX",".png"), width = 5, height = 3.4, units = "in",last_plot(),dpi = 400 ,bg = "White")
netAnalysis_signalingRole_scatter(AGMPSSG, pathways.show)
ggsave(filename = paste0(pathways.show, "_AGMPG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")


netAnalysis_signalingRole_network(AGMPSSG, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(AGMPSSG, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(AGMPSSG, pathways.show)
netVisual_chord_cell(AGMPSSG, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
plotGeneExpression(cellchat, signaling = pathways.show)
FeaturePlot(AGM, c("NOTCH1","NOTCH4","JAG1","DLL4"), label = T, label.size = 2)

pairLRGT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE)
DotPlot(AGM, features = pairLRGT$geneLR) + RotatedAxis()


ggsave(filename = paste0(pathways.show, "_AGMPG_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_ROLE",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_HEAT",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_FP",".png"), last_plot(),dpi = 800 , width = 5, height = 5, units = "in",  bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPG_DOTS",".png"), last_plot(),dpi = 400 ,bg = "White")


pairLR <- extractEnrichedLR(cellchatMKG, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[c(1,2,3,4,5,6,7,8,17,18,19,20,21,22,23,24),]

pathways.show = "WNT"

netVisual_heatmap(cellchatMKG, signaling = pathways.show, color.heatmap = "Reds",)
netVisual_bubble(cellchatMKG, signaling = pathways.show, targets.use = "3_HE" )
netVisual_bubble(cellchatMKG, targets.use = "3_HE", pairLR.use = new_data )
ggsave(filename = paste0(pathways.show, "_AGMPGHE_BUBBLE",".png"), width = 4, height = 4, units = "in",last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPGHE_BUBBLE",".png"),last_plot(),dpi = 500 ,bg = "White")
plotGeneExpression(AGMPSSG, signaling = pathways.show)
ggsave(filename = paste0(pathways.show, "_AGMPG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")
pairLRGT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = TRUE)
DotPlot(AGM, features = pairLRGT$geneLR) + RotatedAxis()
ggsave(filename = paste0(pathways.show, "_AGMPG_DOTSEX",".png"), width = 5, height = 3.4, units = "in",last_plot(),dpi = 400 ,bg = "White")
netAnalysis_signalingRole_scatter(AGMPSSG, pathways.show)
ggsave(filename = paste0(pathways.show, "_AGMPG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
FeaturePlot(AGM, c("WNT2B","NOTCH4","JAG1","DLL4"), label = T, label.size = 2)

netAnalysis_signalingRole_network(AGMPSSG, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(AGMPSSG, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(AGMPSSG, pathways.show)
netVisual_chord_cell(AGMPSSG, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
plotGeneExpression(cellchat, signaling = pathways.show)
FeaturePlot(AGM, c("NOTCH1","NOTCH4","JAG1","DLL4"), label = T, label.size = 2)


pathways.show = "TNF"
png(filename =  paste0(pathways.show, "_AGMPSG_HEAT",".png"), width = 5.04, height = 3.72, units = "in", res = 400)
netVisual_heatmap(cellchatMKG, signaling = pathways.show, color.heatmap = "Reds", font.size = 10,)
dev.off()
netVisual_bubble(cellchatMKG, signaling = pathways.show, targets.use = "3_HE", remove.isolate = T)
ggsave(filename = paste0(pathways.show, "_AGMPSG_BUBBLE",".png"),last_plot(),dpi = 500 ,bg = "White")
netVisual_bubble(cellchatMKG, targets.use = "3_HE", pairLR.use = new_data , remove.isolate = T)
FeaturePlot(AGM, c("WNT2B","NOTCH4","JAG1","DLL4"), label = T, label.size = 2)


pairLR <- extractEnrichedLR(cellchatMKG, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[5,] # show one ligand-receptor pair

genes <- strsplit(LR.show, "_")[[1]]
FeaturePlot(AGM, features = c(genes[1], genes[2]), label = T, label.size = 2)
FeaturePlot(AGM, features = c(genes[1], genes[2], genes[3], "FZD4"), label = T, label.size = 1.5)
ggsave(filename = paste0(pathways.show, "_AGMPSG_FP",".png"), last_plot(),dpi = 500 ,bg = "White")

AGM <- SetIdent(AGM, value = "cell.typestroma")
AGM<- RenameIdents( AGM, c("9_Epi" = "Epi",    "8_Peri"  = "Peri",    "7_MK"   = "MK",  "6_Ery" = "Ery",    "5_Hem"  = "Hem",   "4_HPC" = "HPC",      "3_HE" =   "HE",  "2_Endo" = "Endo", "1_Stromal" = "Stroma" ))


ggsave(filename = paste0(pathways.show, "_AGMPSG_BUBBLE",".png"), last_plot(),dpi = 500 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_ROLE",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_HEAT",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_FP",".png"), last_plot(),dpi = 800 , width = 5, height = 5, units = "in",  bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_AGMPSG_DOTS",".png"), last_plot(),dpi = 400 ,bg = "White")


pathways.show <- "PARs"

"Liver Clusters"
netVisual_bubble(cellchatALC, signaling = pathways.show, targets.use = "HE")
netAnalysis_signalingRole_network(cellchatALC, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(cellchatALC, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(cellchatALC, pathways.show)
netVisual_chord_cell(cellchatALC, signaling = pathways.show, group = group.cellTypeAL, title.name = paste0(pathways.show, " signaling network"))
plotGeneExpression(cellchatALC, signaling = pathways.show)
ggsave(filename = paste0(pathways.show, "_ALC_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")
"Liver Celltype"
netVisual_bubble(cellchatALG, signaling = pathways.show)
netAnalysis_signalingRole_network(cellchatALG, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(cellchatALG, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(cellchatALG, pathways.show)
ggsave(filename = paste0(pathways.show, "_ALG_BUBBLE",".png"), last_plot(),dpi = 400 ,bg = "White")
netVisual_chord_cell(cellchatALG, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
plotGeneExpression(cellchatALG, signaling = pathways.show)
ggsave(filename = paste0(pathways.show, "_ALG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")

netVisual_bubble(cellchatALGP, signaling = pathways.show, targets.use = "3_HE", remove.isolate = T)
netAnalysis_signalingRole_network(cellchatALGP, signaling = pathways.show, width = 10, height = 4, font.size = 8)
netVisual_heatmap(cellchatALGP, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_signalingRole_scatter(cellchatALGP, pathways.show)
ggsave(filename = paste0(pathways.show, "_ALGP_BUBBLE",".png"), last_plot(),dpi = 400 ,bg = "White")
netVisual_chord_cell(cellchatALGP, signaling = pathways.show, title.name = paste0(pathways.show, " signaling network"))
plotGeneExpression(cellchatALGP, signaling = pathways.show)
ggsave(filename = paste0(pathways.show, "_ALGP_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")








#AL <- SetIdent(AL, value = "cell.type")
pairLRGT <- extractEnrichedLR(cellchatALC, signaling = pathways.show, geneLR.return = TRUE)
pairLRGT <- extractEnrichedLR(cellchatALG, signaling = pathways.show, geneLR.return = TRUE)
DotPlot(AL, features = pairLRGT$geneLR) + RotatedAxis()

FeaturePlot(AL, features = pairLRGT$geneLR)

ggsave(filename = paste0(pathways.show, "_ALG_DOTSEX",".png"), last_plot(),dpi = 400 ,bg = "White")


#Winners 
#Fibroblasts 

"BMP4, VEGFA, VEGFB, WNT2, RA"



DotPlot(AL, features = pairLRGT$geneLR)

pathways.show <- c("BMP")


AL <- SetIdent(AL, value = "cell.type")

FeaturePlot(AL,c("WNT2","WNT2B","FZD6"), label = T, repel = T)
ggsave(filename = paste0(pathways.show, "_ALG_CHORD",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_ALG_HEAT",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_ALG_GEX",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_ALG_SCATTER",".png"), last_plot(),dpi = 400 ,bg = "White")
ggsave(filename = paste0(pathways.show, "_ALG_FP",".png"), last_plot(),dpi = 400 ,bg = "White")


FeaturePlot(AGM,c("IGF1","IGF2"), label = T)

FeaturePlot(AGM, features = c(genes[1], genes[2]))
VlnPlot(AGM, features = c(genes[1], genes[2], "TNFRSF1B"))
DotPlot(AGM, features = c(genes[1], genes[2], "TNFRSF1B"))
DotPlot(AGM, features = c("NOTCH1","NOTCH2","NOTCH3","NOTCH4")) + FontSize(6)
DotPlot(AGM, features = c("DLK1","DLL1","DLL3","DLL4","JAG1","JAG2")) +FontSize(6)

DotPlot(AGM, features = c("BMP4","CXCL12","WNT3A","DLK1")) + FontSize(6)

ggsave(filename = paste0(pathways.show, "_AL_DOTS",".png"), last_plot(),dpi = 400 ,bg = "White")


DotPlot(AL, features = c("WNT2","FZD6","FZD4","LRP6")) + FontSize(6)

DotPlot(AL, features = c("BMP4","CXCL12","WNT2","VEGFA","VEGFB", "IGF2", "TGFB2")) + FontSize(6)
DotPlot(AL, features = c("IGFBP3","WNT2","VEGFA","VEGFB", "WNT2B")) + FontSize(6)



DotPlot(AGM, features = c("BMP5","CXCL12","IGFBP3","IGFBP2","IGF2","IGF1")) + FontSize(6) + RotatedAxis()

DotPlot(AGM, features = c("RUNX1","CSPG4","ACTA2","MCAM","RGS5","PDGFA","PDGFB","CD248")) + FontSize(6) + RotatedAxis()
ggsave(filename = paste0("IGFDOTS_DOTS",".png"), last_plot(),dpi = 400 ,bg = "White")




#The Linear Flow of Time 
AGML <- subset(AGM, idents = c("zero", "1", "2","3","4","5","6","7","8","9","10","11","12","13","14","HPC","16","17","18","HE","ERY","21","MK","23","24","25","26","27","28"))


AGM45 <- subset(AGML, subset = orig.ident == "AGM_4_5week")
AGM5 <- subset(AGML, subset = orig.ident == "AGM_5week")
AGM55 <- subset(AGML, subset = orig.ident == "AGM_5_5week")
AGM6 <- subset(AGML, subset = orig.ident == "AGM_6week")

AGM5 <- subset(AGML, subset = orig.ident == "AGM_5week" | orig.ident == "AGM_5_5week")

my.seuObj@meta.data$sample.New[which(my.seuObj@meta.data$sample == "Control")] <- "A"
AGML@meta.data$orig.ident[which(AGML@meta.data$orig.ident == "AGM_5_5week")] <- "AGM_5week"

FeaturePlot(AGML, "SPP1", split.by = "orig.ident")
#Merged 



merged.group <- factor(AGMSSC@meta[["namedclusters"]], 
                       levels = c("zero", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "HPC", "16", "17", "18", "HE", "ERY", "21", "MK", "23", "24", "25", "26", "27", "28", "29"),
                       labels = c("zero", "1", "Merged", "3", "4", "5", "Merged", "7", "8", "9", "10", "11", "12", "13", "14", "HPC", "16", "17", "18", "HE", "ERY", "21", "MK", "23", "24", "25", "26", "27", "28", "29"))

AGMSSCMerge<- mergeInteractions(AGMSSC, merged.group)


#MARKERS OVER TIME 
#Raw #produces the same result 
markersHPC <- FindMarkers(AGM, ident.1 = "AGM_4_5week", ident.2 = "AGM_5week", group.by = 'orig.ident', subset.ident = "HE")

#broken down 

AGMEHT <- subset(AGM, idents = c("HE","HPC","MK"))
DimPlot(AGMEHT, split.by = "orig.ident")
FeaturePlot(AGMEHT, "CXCR4", split.by = "orig.ident")


AGMEHT$clusters <- Idents(AGMEHT)
AGMEHT<- RenameIdents(AGMEHT, "MK" = "MK", "HPC" = "EHT", "HE" = "EHT")
AGMEHT$combo <- Idents(AGMEHT)
AGMEHT<- SetIdent(AGMEHT, value = "clusters")

markersEHT <- FindMarkers(AGMEHT, ident.1 = "AGM_4_5week", ident.2 = "AGM_5week", group.by = 'orig.ident', subset.ident = "HE")


AGMEHT<- SetIdent(AGMEHTalt, value = "clusters" )
AGMEHT<- SetIdent(AGMEHTalt, value = "combo" )
VlnPlot(AGMEHT, "RUNX1", split.by = "orig.ident")
FeaturePlot(AGMEHT, "MECOM", split.by = "orig.ident")
DimPlot(AGMEHT, split.by = "orig.ident")

################################

secret_pathways <- df.All[df.All$`annotation` == "Secreted Signaling", "pathway_name"]

# If 'annotation type' has spaces, use 'annotation.type' instead
# secret_pathways <- df[df$`annotation.type` == "secreted signalling", "signalling pathway"]

secpaths <- unique(secret_pathways)



AGM <- SetIdent(AGM, value = "namedclusters")

library(scMayoMap)
library(Census)
library(ggrepel)

plan("multiprocess", workers = 4)

seurat.markers <- FindAllMarkers(AGM)
scMayoMap.obj <- scMayoMap(data = AGM.markers, database=scMayoMapDatabase)

plt <- scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
scMayoMap.plot(scMayoMap.object = scMayoMap.obj)
View(scMayoMap.obj[["res"]])


Mayo.markers <- FindAllMarkers(AGM6, method = 'MAST')
scMayoMap.obj <- scMayoMap(data = Mayo.markers, database=scMayoMapDatabase)
scMayoMap.plot(scMayoMap.obj)

gns <- scMayoMap.obj$markers$genes[scMayoMap.obj$markers$cluster==19 & scMayoMap.obj$markers$celltype=='kidney:Endothelial cell']
gns <- strsplit(gns, ',')[[1]]
DotPlot(AGM, features = gns)



#Just what is the AGM object? 
#Before we talk about intercellular communications it is best to actually figure out
#what those cells are. Now cell annotation is more of an art then a science and the fact
#we are dealing with embryonic tissues only complicates things as their marker genes
#and RNA profiles are going to be different compared to adults, but we play the cards
#we have not the cards we want. 

#Methods. SingleR, Decouplr, scMayoMap, census. No one tool is perfect except maybe scannoX 
#but the idea is to go for general tissue types and a few oddballs 

res = census_main(LIVX, organ = 'Liver')

des = census_main(AL)

# Plot Census annotations
gg1 <- ggplot(res$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = .1, shape = 16) + 
  geom_label_repel(data = res$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 2, color = 'black', 
                   label.padding = 0.05, 
                   max.overlaps = 200) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(res$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = '5', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))



ggplot(des$pred, aes(umap1, umap2, color = celltype)) + 
  geom_point(size = 0.1, shape = 16) + 
  geom_label_repel(data = des$pred %>% 
                     group_by(celltype) %>% 
                     summarize(u1 = mean(umap1), u2 = mean(umap2)), 
                   force = 100, 
                   aes(u1,u2, label=celltype), 
                   size = 2, color = 'black', 
                   label.padding = 0.1, 
                   max.overlaps = 20) +
  scale_color_manual(values=colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(des$pred$celltype)))) +
  theme_classic() + 
  theme(legend.position = '5', 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(color = 'black', size = 7))



ggplot(des$pred, aes(umap1, umap2, color = celltype)) +
  geom_point(size = 0.1, shape = 16) +  # Adjust geom based on your data
  geom_label_repel(data = des$pred %>%  # No changes to label generation
                     group_by(celltype) %>%
                     summarize(u1 = mean(umap1), u2 = mean(umap2)),
                   force = 100,
                   aes(u1,u2, label=celltype),
                   size = 2, color = 'black',
                   label.padding = 0.1,
                   max.overlaps = 20) +
  scale_color_manual(values = colorRampPalette(brewer.pal(10, "Spectral"))(length(unique(des$pred$celltype))),
                     # Legend customizations
                     breaks = unique(des$pred$celltype),  # Explicitly set legend breaks
                     labels = unique(des$pred$celltype),  # Set legend labels
  ) +
  theme_classic()

#ScX
library(scX)
cseo <- createSCEobject(xx = AGM, 
                        partitionVars = "fg")

cseo <- createSCEobject(xx = AGM, partitionVars = "namedclusters", metadata = AGM@meta.data ,metadataVars = "orig.ident", calcRedDim = FALSE, markerList = markersX)

launch_scX(cseo)




#markers

#Fcoex #doesnt work

library(fcoex)
library(FCBF)

exprs <- data.frame(GetAssayData(AGM))
target <- Idents(AGM)
fc <- new_fcoex(data.frame(exprs),target)
fc <- discretize(fc)
fc <- find_cbf_modules(fc,n_genes = 50, verbose = TRUE, is_parallel = TRUE)
fc <- get_nets(fc)
mod_names(fc)
network_plots <- show_net(fc)
network_plots[["MYB"]]


AGM <- SetIdent(AGM, value = "fc")

AGM <- RenameIdents(object = AGM, "zero"="Muscle", "1"="Endo", "2"="Endo", "3"="Endo", "4"="Fibro", "5"="Fibro", "6"="Macro", "7"="Endo", "8"="Lymph", "9"="Fibro", "10"="Endo", "11"="Peri", "12"="Fibro", "13"="Gran", "14"="Macro", "HPC"="HPSC", "16"="Fibro", "17"="Muscle", "18"="Macro", "HE"="HE", "ERY"="ERY", "21"="Kidney", "MK"="MK", "23"="Fibro", "24"="Liver", "25"="Kidney", "26"="Kidney", "27"="Fibro", "28"="Mel")
AGM@meta.data$fg <- AGM@active.ident




netVisual_bubble(cellchatAL, signaling = pathways.show, targets.use = c("HE"))
netVisual_heatmap(cellchatAL, signaling = pathways.show,remove.isolate = T)



pathways.show = "WNT"
netVisual_bubble(cellchatAGMSSG, signaling = pathways.show, targets.use = c("2_Endo","3_HE"), remove.isolate = T)
netVisual_heatmap(cellchatAGMSSG, signaling = pathways.show)

netVisual_bubble(cellchatAGMC, signaling = pathways.show, targets.use = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"))
netVisual_bubble(cellchat, signaling = pathways.show, targets.use = c("HE", "Endo2"))
netVisual_heatmap(cellchat, signaling = pathways.show,remove.isolate = T)

netVisual_bubble(cellchatAGMG, signaling = pathways.show, remove.isolate = T)

netVisual_heatmap(cellchatAGMC, signaling = pathways.show)
netVisual_heatmap(cellchatAGMG, signaling = pathways.show)
netVisual_heatmap(cellchat, signaling = pathways.show, remove.isolate = T)
netVisual_heatmap(cellchat, signaling = pathways.show, remove.isolate = T, targets.use = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5","HPSC"))


netVisual_bubble(cellchatAL, signaling = pathways.show)
netVisual_heatmap(cellchatAL, signaling = pathways.show)
VlnPlot(AGM, "BMP4", pt.size = 0)

VlnPlot(AGMi, "BMP4", pt.size = 0)
VlnPlot(AGM, "PROCR", idents = "HE", split.by = "orig.ident")
FeaturePlot(AGM, "BMP4")

geneX= "CXCR4"
geneG= c("PARD3","F2R")
VlnPlot(AGM, geneX, pt.size = 0)
VlnPlot(AGMi, geneX, pt.size = 0)
VlnPlot(AGMi, geneX, pt.size = 0, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5","HSC"))
FeaturePlot(AGM, geneX)


VlnPlot(AGM, geneX, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"), pt.size = 0)
VlnPlot(AGMi, geneG, idents = c("HE", "Endo1","Endo2","Endo3","Endo4","Endo5"), pt.size = 0)
VlnPlot(AGM, geneG,  pt.size = 0)

AverageExpression(AGM, features = "TNFRSF1B")
CellChat::computeAveExpr(cellchat, features = "TNFRSF1B")
 
#WAOW
netAnalysis_signalingRole_network(cellchat, signaling = c("CXCL","ANGPTL"))
netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL"))
netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "NOTCH","WNT","TNF","TGFb"))
netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "NOTCH","WNT","TNF","TGFb","ANGPTL","SPP1"), pattern = "incoming")
netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "NOTCH","WNT","TNF","TGFb","ANGPTL","SPP1","BMP"), pattern = "outgoing")
VlnPlot(AGM_harmony_merged_filt, "F2R", pt.size = 0, idents = c(1,2,5,9,12,19), split.by = "seurat_clusters", group.by = "orig.ident")


netAnalysis_signalingRole_heatmap(cellchat, "ANGPTL", "incoming" )


netAnalysis_signalingRole_scatter(cellchatAGMC, pathways.show)


plotGeneExpression(cellchatAGMC, signaling = pathways.show)

#netVisual_hierarchy1(cellchat@net$weight, vertex.receiver = (1,4))
netVisual_aggregrate()
netVisual_aggregate(cellchat, signaling = "WNT", layout = "hierarchy", vertex.receiver = c(2,3,4,8,11,19))

netVisual_heatmap(cellchat, signaling = pathways.show,  remove.isolate = T, targets.use = c("HE","Endo1","Endo2","Endo3","Endo4","Endo5"), sources.use = c("Endo1" ,  "Endo2",   "Endo3",   "Endo4"  , "Endo5" ,  "ERY"  ,   "Fibro1"  ,"Fibro2",  "Fibro3",  "Fibro4", 
"Fibro5",  "Fibro6",  "Fibro7",  "Gran",    "HE",      "HPSC",    "Kidney1", "Kidney2", "Kidney3", "Liver",  
 "Lymph",   "Macro1",  "Macro2",  "Macro3",  "Mel",       "Muscle1", "Muscle2", "Peri"))  




df <- data.frame(
  gene = rep(paste0("Gene", 1:10), 3),
  cluster = rep(c("Cluster1", "Cluster2", "Cluster3"), each = 10),
  expression = runif(30)
)


cell_type_annotations <- data.frame(
  y = c(2, 5, 8),
  label = c("Cell Type A", "Cell Type B", "Cell Type C")
)

ggplot(df, aes(x = cluster, y = gene, size = expression)) +
  geom_point() +
  geom_text(data = cell_type_annotations, aes(x = Inf, y = y, label = label), hjust = 0, vjust = 0.5)












