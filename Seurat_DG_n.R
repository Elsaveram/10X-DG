# install.packages('Seurat')

library(Seurat)
library(dplyr)
library(Matrix)
library(stringr)

# This assumes your script is in a folder called "src" at the top level of all the results.
project_folder <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), ".." )
analysis_folder <- file.path(project_folder,"Results/Sample234")

#Load Seurat object if it exist
load( file = file.path(project_folder, 'Negative_Low_High_filtering_h_0417.Robj') )
load( file = file.path(project_folder, 'Negative_Low_High_filtering_h_AllMarkers_res1_4_0418.Robj') )

# Load the dataset
pbmc.data <- Read10X("/Volumes/Elsa Work/10X/Results/Sample234/outs/filtered_gene_bc_matrices_mex/mm10ECEV4806")

dim(pbmc.data)

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200, project = "10X_DG")

# The number of genes and UMIs (nGene and nUMI) are automatically calculated
# for every object by Seurat.  For non-UMI data, nUMI represents the sum of
# the non-normalized values within a cell We calculate the percentage of
# mitochondrial genes here and store it in percent.mito using AddMetaData.
# We use object@raw.data since this represents non-transformed and
# non-log-normalized counts The % of UMI mapping to MT-genes is a common
# scRNA-seq QC metric.  NOTE: You must have the Matrix package loaded to
# calculate the percent.mito values.
mito.genes <- grep(pattern = "^mt-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to
# stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("Actb","Gapdh", "Ubb"), nCol = 3)
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
# GenePlot is typically used to visualize gene-gene relationships, but can
# be used for anything calculated by the object, i.e. columns in
# object@meta.data, PC scores etc.  Since there is a rare subset of cells
# with an outlier level of high mitochondrial percentage and also low UMI
# content, we filter these as well
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")
table(pbmc@meta.data$percent.mito < 0.05 & pbmc@meta.data$nGene<2500)
table(pbmc@meta.data$percent.mito < 0.2 & pbmc@meta.data$nGene<4000)

# We filter out cells that have unique gene counts over 2,500 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate' -Inf and Inf should be used if you don't want a lower or upper
# threshold.
#For GFP positive and negative
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene","nUMI"), 
                    low.thresholds = c(600, 800), high.thresholds = c(Inf, 20000))

#For all the samples
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene","nUMI", "percent.mito"), 
                    low.thresholds = c(1500, 4, -Inf), high.thresholds = c(4000, 20000, 0.1))

# Normalize the data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


pbmc <- FilterCells(object = pbmc, subset.names = c("Actb","Gapdh", "Ubb"), 
                    low.thresholds = c(0.5, 0.5, 0.5), high.thresholds = c(Inf, Inf, Inf))


VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito","Actb","Gapdh", "Ubb"), nCol = 3)

# Find variable genes
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

#For GFP positive and negative
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI"))

# Scale the data
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

# Run and Predict PCA
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

# Examine and visualize PCA results a few different ways
PrintPCA(object = pbmc, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
pbmc <- ProjectPCA(object = pbmc, do.print = FALSE)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)

# Generate a heatmap
PCHeatmap(object = pbmc, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
PCHeatmap(object = pbmc, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)

# Find the number of PCs (clusters) to be used in our analysis
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)
JackStrawPlot(object = pbmc, PCs = 1:100)
PCElbowPlot(object = pbmc)

#Clustering

# Add sample data to our Seurat object
sample.numbers <- as.numeric(str_extract(colnames(pbmc@data), "[^-]+$"))
names(sample.numbers) <- colnames(pbmc@data)
pbmc <- AddMetaData(pbmc, sample.numbers, "sample.number.column")

# Map the sample numbers into a human readable form
current.sample.ids <- c(1, 2, 3)
new.sample.ids <- c("High", "Negative", "Low")
sample.number.names <- plyr::mapvalues(pbmc@meta.data$sample.number.column, from = current.sample.ids, to = new.sample.ids)
names(sample.number.names) <- colnames(pbmc@data)
pbmc <- AddMetaData(pbmc, sample.number.names, "sample.number.names")

#Name the clusters
pbmc <- SetAllIdent(pbmc, "res.0.3")
current.clusters <- 0:15
new.clusters <- c("NSCs", "OPCs", "NSCs", "mature OLs", "early NPCs", 
               "commited OPC", "pericytes", "mature granular neurons", "late NPCs", "micoglia", 
               "activated OPC", "interneurons", "OL progenitor?", "neuroblast?", "endothelial cells", 
               "pericytes" )
for (i in current.clusters) {
  pbmc <- RenameIdent(object = pbmc, old.ident.name = i, 
                                 new.ident.name = new.clusters[i + 1])
}

cluster.names <- plyr::mapvalues(pbmc@meta.data$res.0.3, from = current.clusters, to = new.clusters)
names(cluster.names) <- colnames(pbmc@data)
pbmc <- AddMetaData(pbmc, cluster.names, "Clusters")

TSNEPlot(pbmc, do.label = T, pt.size = 1)

save (pbmc, file = file.path(project_folder, 'Samples234_Artegianifilter.Robj') )

#Name the clusters For Positive And Negative
pbmc <- SetAllIdent(pbmc, "res.1.4")
current.clusters <- 0:22
new.clusters <- c("astrocytes", "astrocytes", "astrocytes", "astrocytes", "astrocytes", 
                  "mature OL", "microglia", "pericytes", "immature granule", "OPC", 
                  "astrocytes", "RGL-activated", "interneurons", "mature OL", "immature granule", 
                  "committed OPC","RGL-quiescent", "immature granule", "mature granule", "endothelial", "IPC", "pericytes", "macrophages")
for (i in current.clusters) {
  pbmc <- RenameIdent(object = pbmc, old.ident.name = i, 
                      new.ident.name = new.clusters[i + 1])
}

cluster.names <- plyr::mapvalues(pbmc@meta.data$res.1.4, from = current.clusters, to = new.clusters)
names(cluster.names) <- colnames(pbmc@data)
pbmc <- AddMetaData(pbmc, cluster.names, "Clusters")
TSNEPlot(pbmc, do.label = T, pt.size = 0.5, no.axes=TRUE)

save (pbmc, file = file.path(project_folder, 'Negative_Low_High_filtering_h_AllMarkers.Robj') )

#Name the clusters For Negative, Low, High
pbmc <- SetAllIdent(pbmc, "res.1.4")
current.clusters <- 0:25
new.clusters <- c("astrocytes", "astrocytes", "OPC", "astrocytes", "astrocytes", 
                  "mature OL", "RG-like active", "microglia", "immature neurons", "pericytes", 
                  "neuroblast", "astrocytes", "mature OL", "young OL", "immature granule", 
                  "activated OPC","mature granule", "NPC", "interneurons", "young OL", "early OPC", "RG-like quiescent", "endothelial", "macrophages", "VLMC", "OPC")
for (i in current.clusters) {
  pbmc <- RenameIdent(object = pbmc, old.ident.name = i, 
                      new.ident.name = new.clusters[i + 1])
}

cluster.names <- plyr::mapvalues(pbmc@meta.data$res.1.4, from = current.clusters, to = new.clusters)
names(cluster.names) <- colnames(pbmc@data)
pbmc <- AddMetaData(pbmc, cluster.names, "Clusters")
TSNEPlot(pbmc, do.label = T, pt.size = 0.5, no.axes=TRUE)

save (pbmc, file = file.path(project_folder, 'Negative_Low_High_filtering_h_0417.Robj') )
# Why do this now?
# VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI"), nCol = 2)

head(pbmc@meta.data)

# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value (see docs for
# full details)
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:15, 
                     resolution = 1.4, print.output = 0, save.SNN = TRUE)
PrintFindClustersParams(object = pbmc)
pbmc <- RunTSNE(object = pbmc, dims.use = 1:15, do.fast = TRUE)
TSNEPlot(object = pbmc, do.label=TRUE, no.axes=TRUE, pt.size = 0.5 )


#To set new identity to samples

pbmc <- SetAllIdent(object = pbmc, id = "sample.number.names")
TSNEPlot(object = pbmc, colors.use=c("dark green", "green", "blue"), no.axes=TRUE,  pt.size = 0.5, do.label=FALSE )
table(pbmc@meta.data$sample.number.column)
pbmc <- SetAllIdent(pbmc, "res.0.3")

#To set new identity to samples for GFP Neg and Pos

pbmc <- SetAllIdent(object = pbmc, id = "sample.number.names")
TSNEPlot(object = pbmc, colors.use=c("blue", "dark green"), no.axes=TRUE,  pt.size = 0.5, do.label=FALSE )
table(pbmc@meta.data$sample.number.column)


# Count of cells by sample and cluster (with 0.6 resolution) Doesn't work
pbmc <- SetAllIdent(pbmc, "Clusters")
Transcriptpercluster<- with(pbmc@meta.data, table(sample.number.column, Clusters))
write.csv(Transcriptpercluster, file = file.path(project_folder, 'Negative_Low_High_filtering_h_TranscriptsClusters0418.csv') )

#Violn plot for different samples
VlnPlot(object = pbmc, features.plot = c("sample.number.column"), point.size.use = .25, group.by ="res.0.6" )


#For dark background in the tSNE plot dark.theme= TRUE

#Finding differentially expressed genes (cluster biomarkers)
pbmc <- SetAllIdent(pbmc, "Clusters")
pbmc <- SetAllIdent(pbmc, "res.1.4")
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)

save (pbmc.markers, file = file.path(project_folder, 'Negative_Low_High_filtering_h_AllMarkers_res1_4_0418.Robj'))


#Vln plot ny sample type (doesn't work)
DotPlot(object = pbmc, features.plot = c("sample.number.column"), point.size.use = .5)


#Sungle tSNE graoh

FeaturePlot(object = pbmc, features.plot = c("ECEV4806transgene"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne",  no.axes=TRUE,  pt.size = 0.5)

#pairwise comparisons of gene expression RGL.and.IPC.early <- c("RG-like quiescent", "RG-like active", "NPC", "neuroblast")
pbmc <- SetAllIdent(pbmc, "Clusters")
RG_like_quiescent <- SubsetData(pbmc, ident.use = "RG-like quiescent", subset.raw = T)
RG_like_active <- SubsetData(pbmc, ident.use = "RG-like active", subset.raw = T)
avg.RG_like_quiescent <- log1p(AverageExpression(RG_like_quiescent, show.progress = FALSE))
avg.RG_like_quiescent$gene <- rownames(avg.RG_like_quiescent)
avg.RG_like_active <- log1p(AverageExpression(RG-like_active, show.progress = FALSE))
avg.RG_like_quiescent$gene <- rownames(avg.RG_like_active)

genes.to.label1 = c("Ascl1" , "Dbi","Kcne1l", "Nudt4", "Mia", "Sox4")

p1 <- ggplot(aes(avg.RG_like_quiescent, avg.RG_like_active)) + geom_point()
p1 <- LabelUR(p1, genes = genes.to.label1, avg.t.cells, 
              adj.u.t = 0.3, adj.u.s = 0.23)

#other.lineage.clusters <- c("immature neurons", "interneurons", "pericytes", "microglia", "endothelial", "macrophages", "VLMC")
#Ben's code for Violin plot and Feature plot for loop
genes_of_interest <- list(
 
  #'NSC lineage'              = c("Sox9", "Id4", "Clu", "Hmgn2", "Ccnd2", "Neurod1", "Cd24a", "Syt5", "ECEV4806transgene" ),
  'gene_select'              = c("ECEV4806transgene","Aqp4", "Gfap", "Aldh1l1", "Sox9", "Lpar1", "Sox2", "Nes", "Prom1", "Top2a", "Cdk1", "Ccnd2", "Dcx", "Sox11", "Ascl1" , "Dbi","Kcne1l", "Nudt4", "Mia", "Sox4", "Ube2c", "Islr2")
  'NSC lineage'              = c("Sox9", "Aldh1l1", "Lpar1", "Nes", "Cdk1", "Eomes", "Neurod1" , "Rbfox3","ECEV4806transgene" ),
  'Transgene'                = c("ECEV4806transgene"),
  'Stem Cells'               = c("Aldoc", "Apoe", "Id4", "Hopx", "Sox9", "Gfap", "Slc1a3", "Sox2", "Fabp7", "Clu"),
  'NSC lineage'              = c("Sox9", "Id4", "Clu", "Hmgn2", "Ccnd2", "Neurod1", "Cd24a", "Syt5", "ECEV4806transgene" ),
  'Oligodendrocyte lineage'  = c( "Sox10", "Olig2", "Mki67", "Pdgfra", "Cspg4", "Slc1a3", "Mbp", "Mal", "Bmp4"), 
  'Other cell types'         = c( "Tagln", 'Vwf', "Csf1r", "Cx3cr1", "Reln", "Ndnf", "S100b", "Fzd2"), 
  'Early Neural progenitors' = c("Eomes","Ccnd2","Hmgn2"),
  'Late Neural progenitors'  = c('Neurod1', 'Dcx', "Sox11"),
  'Neural progenitors'       = c("Eomes","Ccnd2", 'Neurod1', 'Dcx', "Sox11", "Hmgn2"),
  'Young neurons'            = c("Dcx","Calb2",'Cd24a'),
  'Mature granular neurons'  = c("Calb1", "Gria1"),
  'Mature neurons'           = c("Rbfox3","Syt1","Syt5","Myt1l"),
  'OPC'                      = c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"),  
  'Commited OP'              = c("Bmp4", "Fyn", "Gpr17"),
  'Mature Oligodendrocytes'  = c("Plp1", "Mal", "Mog", "Mbp"),  
  'Microglia'                = c("Csf1r", "Cx3cr1"),
  'Macrophages'              = c('Ptprc'),
  'Endothelial cells'        = c('Vwf'),
  'Astrocytes'               = c("S100b", "Fzd2"),
  'Proliferation'            = c("Ccnd2", "Mki67", "Pcna", "Mcm2"),
  'Interneurons'             = c("Reln", "Ndnf"),
  'Pericytes'                = c("Tbx18", "Vtn", "Tagln", "Des"),
  'Universally expressed'    = c("Actb", "Gapdh", "Ubb"), 
  'CD4 T Cells'              = c("Il7r"), 
  'RBC'                      = c("Hba-a1", "Hba-a2"), 
  'Transgene_h'                = c("ECEV4806transgene"),
  'Astrocytes_h'               = c("Gfap", "Hes5", "Sox9", "Aqp4"),
  'RGL_h'                      = c("Lpar1", "Nes", "Prom1", "Hes5", "Tfap2c", "Rhcg", "Wnt8b"),
  'nIPC_h'                     = c("Cdk1","Aurkb", 'Top2a', "Neurog2", "Eomes", "Neurod4", "Tfap2c"), 
  'NB1_h'                      = c( "Eomes", "Tac2", "Calb2", "Igfbpl1", "Neurod4"), 
  'NB2_h'                      = c("Gal","Sox11","Dcx", "Igfbpl1"),
  'Granule inmature_h'         = c("Fxyd7"),
  'Granule mature_h'           = c("Plk5", "Ntng1"), 
  'Cell Cycle genes_h'         = c("Cdk1","Aurkb", 'Top2a', 'Mki67'))    



for ( focus_genes in names(genes_of_interest) ) {
  print(paste('Generating graphs for', focus_genes ))
  genes <- genes_of_interest[[focus_genes]]
  number_of_colums <- max(round(length(genes)/3),1)
  vln_plot_width <- 1000*number_of_colums
  print(paste('Auto calculating columns and plot width to be', number_of_colums, vln_plot_width ))
  
  # Generate the cluster plot
  png( file.path(analysis_folder,"Plots", paste(focus_genes, "Cluster.png")), res=100, height=1000, width=1600)
  FeaturePlot(object = pbmc, features.plot = genes, cols.use = c("grey", "blue"), reduction.use = "tsne")
  dev.off()
  
  png( file.path(analysis_folder,"Plots", paste(focus_genes, "Heatmap.png")), res=100, height=1000, width=2200)
  FeatureHeatmap(pbmc, features.plot = genes, group.by = "sample.number.names", pt.size = 0.25,  key.position = "top", max.exp = 3)
  dev.off()
  
  for ( identity in c( "sample.number.names", "Clusters")) {
    pbmc <- SetAllIdent(object = pbmc, id = identity)
    
    #generate the dot plot
    png( file.path(analysis_folder,"Plots", paste(focus_genes, identity, "Dot.png")), res=100, height=800, width=800)
    DotPlot(pbmc,genes.plot = genes, x.lab.rot = T, dot.scale = 8, plot.legend=T)
    dev.off()
    
    # Generate the violin plot
    png( file.path(analysis_folder,"Plots", paste(focus_genes, identity, "Violn.png")), res=100, height=1000, width=vln_plot_width)
    print(VlnPlot(object = pbmc, do.sort=T, x.lab.rot = T, point.size.use = .5, features.plot = genes, nCol = number_of_colums))
    dev.off()
  }
}


# Get the top 10 genes per cluster 
pbmc <- SetAllIdent(pbmc, "res.1.4")
pbmc <- SetAllIdent(pbmc, "Clusters")
top30 <- pbmc.markers %>% group_by(cluster) %>% top_n(30, avg_logFC)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
top3 <- pbmc.markers %>% group_by(cluster) %>% top_n(3, avg_logFC)
top2 <- pbmc.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)

# Write the top10 genes by cluster to a CSV file
write.csv( top30, file = file.path(project_folder, 'Negative_Low_High_filtering_h_top30markers_Clusters_0418.csv') )
write.csv( top10, file = file.path(project_folder, 'top10markers.csv') )
write.csv( pbmc.markers, file = file.path(project_folder, 'allmarkers.csv') )

# setting slim.col.label to TRUE will print just the cluster IDS instead of
# every cell name
#print heatmap per cluster


oligodendrocytes.lineage.clusters <- c("early OPC", "OPC", "activated OPC", "young OL","mature OL")
neuronal.cluster<-c("immature neurons", "immature granule", "mature granule", "interneurons")
other.lineage.clusters <- c("pericytes", "microglia", "endothelial", "macrophages", "VLMC")
RGL.and.IPC <- c( "astrocytes", "RGL-quiescent", "RGL-activated", "IPC", "neuroblast")
RGL.and.IPC.early <- c("RG-like quiescent", "RG-like active", "NPC", "neuroblast")
Astrocytes <- c("astrocytes","RG-like active","RG-like quiescent")
focus.clusters <- other.lineage.clusters

new.clusters <- 
top_genes <- pbmc.markers %>% group_by(cluster) %>% filter(cluster %in% focus.clusters) %>% mutate(order = match(cluster, focus.clusters)) %>% top_n(10,  avg_logFC) %>% arrange(order)
focus.cells <- WhichCells( pbmc, ident = focus.clusters, max.cells.per.ident = 100)
 
DoHeatmap(object = pbmc, cells.use = focus.cells, genes.use = top_genes$gene, slim.col.label = TRUE, group.order = focus.clusters, col.low="#330099",col.mid = "#000000", col.high = "#CC0000")

?DoHeatmap

DoHeatmap(object = pbmc, genes.use = top30$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = pbmc, genes.use = top5$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = pbmc, genes.use = top3$gene, slim.col.label = TRUE, remove.key = TRUE)
DoHeatmap(object = pbmc, genes.use = top2$gene, slim.col.label = TRUE, remove.key = TRUE)

pbmc.markers(head)

#Get the 35 top genes per sample

pbmc <- SetAllIdent(pbmc, "sample.number.names")
sample.markers <- FindAllMarkers(object = pbmc, only.pos = TRUE, min.pct = 0.25, 
                               thresh.use = 0.25)
top35sample <- pbmc.markers %>% group_by("sample.number.names") %>% top_n(35, avg_logFC)
write.csv( top35sample, file = file.path(project_folder, 'top35sample.csv') )
DoHeatmap(object = pbmc, genes.use = top35sample$gene, slim.col.label = TRUE, remove.key = TRUE)




# Generate 2D plots for gene expression between clsuters

#"RG-like quiescent", "RG-like active"

cluster_1 <- SubsetData(pbmc, ident.use = "RG-like quiescent", subset.raw = T)
avg.cluster_1 <- log1p(AverageExpression(cluster_1))
avg.cluster_1$gene <- rownames(avg.cluster_1)

cluster_2 <- SubsetData(pbmc, ident.use = "RG-like active", subset.raw = T)
avg.cluster_2 <- log1p(AverageExpression(cluster_2))
avg.cluster_2$gene <- rownames(avg.cluster_2)




merged_clusters = merge(avg.cluster_1,avg.cluster_2,by="gene")

merged_clusters['x'] = merged_clusters["RG-like quiescent"]
merged_clusters['y'] = merged_clusters["RG-like active"]


highlight.gene <- c("Igfbp5")
mycolours <- c("highlight" = "red", "normal" = "grey50")

merged_clusters$highlight <- ifelse(merged_clusters$gene == highlight.gene, "highlight", "normal")


ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x) +geom_point(size = 1, aes(colour = highlight))
ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x) 

#"Astrocytes", "RG-like active"

cluster_1 <- SubsetData(pbmc, ident.use = "astrocytes", subset.raw = T)
avg.cluster_1 <- log1p(AverageExpression(cluster_1))
avg.cluster_1$gene <- rownames(avg.cluster_1)

cluster_2 <- SubsetData(pbmc, ident.use = "RG-like active", subset.raw = T)
avg.cluster_2 <- log1p(AverageExpression(cluster_2))
avg.cluster_2$gene <- rownames(avg.cluster_2)


merged_clusters = merge(avg.cluster_1,avg.cluster_2,by="gene")

merged_clusters['x'] = merged_clusters["astrocytes"]
merged_clusters['y'] = merged_clusters["RG-like active"]

highlight.gene <- "Igfbp5"
mycolours <- c("highlight" = "red", "normal" = "grey50")

merged_clusters$highlight <- ifelse(merged_clusters$gene == highlight.gene, "highlight", "normal")
ggplot(merged_clusters, aes(x,y)) + geom_point() + geom_point(size = 1, aes(colour = highlight))

ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x) 

#"Astrocytes", "RG-like quiescent"

cluster_1 <- SubsetData(pbmc, ident.use = "astrocytes", subset.raw = T)
avg.cluster_1 <- log1p(AverageExpression(cluster_1))
avg.cluster_1$gene <- rownames(avg.cluster_1)

cluster_2 <- SubsetData(pbmc, ident.use = "RG-like quiescent", subset.raw = T)
avg.cluster_2 <- log1p(AverageExpression(cluster_2))
avg.cluster_2$gene <- rownames(avg.cluster_2)


merged_clusters = merge(avg.cluster_1,avg.cluster_2,by="gene")

merged_clusters['x'] = merged_clusters["astrocytes"]
merged_clusters['y'] = merged_clusters["RG-like quiescent"]

highlight.gene <- "Hsd11b1"
mycolours <- c("highlight" = "red", "normal" = "grey50")

merged_clusters$highlight <- ifelse(merged_clusters$gene == highlight.gene, "highlight", "normal")
ggplot(merged_clusters, aes(x,y)) + geom_point() + geom_point(size = 1, aes(colour = highlight))
ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x, se=F) 


#"RG-like active", "NPC"

cluster_1 <- SubsetData(pbmc, ident.use = "RG-like active", subset.raw = T)
avg.cluster_1 <- log1p(AverageExpression(cluster_1))
avg.cluster_1$gene <- rownames(avg.cluster_1)

cluster_2 <- SubsetData(pbmc, ident.use = "NPC", subset.raw = T)
avg.cluster_2 <- log1p(AverageExpression(cluster_2))
avg.cluster_2$gene <- rownames(avg.cluster_2)


merged_clusters = merge(avg.cluster_1,avg.cluster_2,by="gene")

merged_clusters['x'] = merged_clusters["RG-like active"]
merged_clusters['y'] = merged_clusters["NPC"]

highlight.gene <- "Cks2"
mycolours <- c("highlight" = "red", "normal" = "grey50")

merged_clusters$highlight <- ifelse(merged_clusters$gene == highlight.gene, "highlight", "normal")
ggplot(merged_clusters, aes(x,y)) + geom_point() + geom_point(size = 1, aes(colour = highlight))
ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x, se=F) 

#"NPC" vss neuroblast

cluster_1 <- SubsetData(pbmc, ident.use = "NPC", subset.raw = T)
avg.cluster_1 <- log1p(AverageExpression(cluster_1))
avg.cluster_1$gene <- rownames(avg.cluster_1)

cluster_2 <- SubsetData(pbmc, ident.use = "neuroblast", subset.raw = T)
avg.cluster_2 <- log1p(AverageExpression(cluster_2))
avg.cluster_2$gene <- rownames(avg.cluster_2)


merged_clusters = merge(avg.cluster_1,avg.cluster_2,by="gene")

merged_clusters['x'] = merged_clusters["NPC"]
merged_clusters['y'] = merged_clusters["neuroblast"]

highlight.gene <- "Islr2"
mycolours <- c("highlight" = "red", "normal" = "grey50")

merged_clusters$highlight <- ifelse(merged_clusters$gene == highlight.gene, "highlight", "normal")
ggplot(merged_clusters, aes(x,y)) + geom_point() + geom_point(size = 1, aes(colour = highlight))

ggplot(merged_clusters, aes(x,y)) + geom_point()+stat_smooth(method='lm', formula=y ~ x, se=F) 

