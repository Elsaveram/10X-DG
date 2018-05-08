
#source("http://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("monocle")
#library(reshape2)
#library(stringr)
#install.packages("devtools")
#devtools::install_github("cole-trapnell-lab/monocle-release@develop")
#biocLite(c("DDRTree", "pheatmap"))

library(monocle)
library(Seurat)

# This assumes your script is in a folder called "src" at the top level of all the results.
project_folder <- file.path(dirname(rstudioapi::getActiveDocumentContext()$path), ".." )
analysis_folder <- file.path(project_folder,"Results/Sample234")

#Load Seurat or Monocle object if it exist
load( file = file.path(project_folder, 'HSMM_Negative_Low_High_filtering_h_0418.Robj') )
load( file = file.path(project_folder, 'PositiveAndNegative.Robj') )

#Import Seurat object
HSMM <- importCDS(pbmc, import_all = T)

# We must calcualate these 
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)

HSMM <- detectGenes(HSMM, min_expr = 0.1)
head(fData(HSMM))
head(pData(HSMM))

# c( 0, 1, 3, 4, 18, 19, 8, 10, 12)
RGL.and.IPC <- c( "RG-like quiescent", "RG-like active", "NPC", "neuroblast")
neurogenic.lineage.clusters <- c( "NSCs", "neuroblast?", "early NPCs", "late NPCs", "mature granular neurons")
oligodendrocytes.lineage.clusters <- c("OL progenitor?", "OPCs", "activated OPC", "commited OPC","mature OLs")
other.lineage.clusters <- c("pericytes", "micoglia", "interneurons", "endothelial cells")

valid_cells <- row.names(subset(pData(HSMM), is.element(Clusters, c(RGL.and.IPC))))
HSMM <- HSMM[,valid_cells]

expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

HSMM <- detectGenes(HSMM, min_expr = 0.1)

# Unsupervised clustering
disp_table <- dispersionTable(HSMM)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
HSMM <- setOrderingFilter(HSMM, unsup_clustering_genes$gene_id)
plot_ordering_genes(HSMM)

# HSMM@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(HSMM, return_all = F) # norm_method='log'

# Cluster plot using TSNE
HSMM <- reduceDimension(HSMM, max_components = 2, num_dim = 20,
                        reduction_method = 'tSNE', verbose = T)

HSMM <- clusterCells(HSMM, num_clusters=4)
plot_cell_clusters(HSMM, color_by="Clusters")
plot_cell_clusters(HSMM)

# Constructing cell trajectories
# Trajectory step 1: choose genes that define a cell's progress
RGL.and.IPC <- c("Lpar1", "Nes", "Prom1", "Hes5", "Tfap2c", "Rhcg", "Wnt8b", "Cdk1","Aurkb", 'Top2a', "Neurog2", "Eomes", "Neurod4", "Tfap2c", "Aldoc", "Apoe", "Id4", "Hopx", "Sox9", "Gfap", "Slc1a3", "Sox2", "Fabp7","Ccnd2", 'Neurod1', "Sox11","Dcx")
neurogenic.lineage.genes <- c("Aldoc", "Apoe", "Id4", "Hopx", "Sox9", "Gfap", "Slc1a3", "Sox2", "Fabp7", "Eomes","Ccnd2", 'Neurod1', "Sox11","Dcx","Calb2", "Calb1", 'Cd24a', "Gria1", "Rbfox3","Syt1","Syt5","Myt1l" )
oligodendrocytes.lineage.genes <- c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4", "Bmp4", "Fyn", "Gpr17","Plp1", "Mal", "Mog", "Mbp")

ordering_genes <- RGL.and.IPC
HSMM <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)


#Trajectory step 2: reduce data dimensionality
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree', num_centers=3)

#Trajectory step 3: order cells along the trajectory
HSMM <- orderCells(HSMM, reverse = NULL)
plot_cell_trajectory(HSMM, color_by="Clusters")
plot_cell_trajectory(HSMM, color_by="State")
plot_cell_trajectory(HSMM, color_by = "Pseudotime")

# #####
#Ordering based on genes that differ between cluster
# #####
HSMM <- detectGenes(HSMM, min_expr = 0.1)
fData(HSMM)$use_for_ordering <-
  fData(HSMM)$num_cells_expressed > 0.05 * ncol(HSMM)

plot_pc_variance_explained(HSMM, return_all = F)
HSMM <- reduceDimension(HSMM,
                            max_components = 2,
                            norm_method = 'log',
                            num_dim = 3,
                            reduction_method = 'tSNE',
                            verbose = T)

HSMM <- clusterCells(HSMM, verbose = F)
plot_cell_clusters(HSMM, color_by = 'as.factor(res.0.6)')
plot_cell_clusters(HSMM, color_by = 'as.factor(State)')
plot_rho_delta(HSMM, rho_threshold = 2, delta_threshold = 4 )
HSMM<- clusterCells(HSMM,
                         rho_threshold = 60,
                         delta_threshold = 10,
                         skip_rho_sigma = T,
                         verbose = F)
plot_cell_clusters(HSMM, color_by = 'as.factor(res.0.6)')
plot_cell_clusters(HSMM, color_by = 'as.factor(State)')




# #####
# Ordering cells using known marker genes
# #####

# Save our objects as they take a VERY long time to generate (with 1 core)
save (HSMM, file = file.path(project_folder, 'HSMM_Negative_Low_High_filtering_h_0418.Robj') )

load( file = file.path(project_folder, 'HSMM_Neurogenic_Lineage_0224.Robj') )

#1. Run "expressed genes" and NSC genes_id when uploading a Monocle object:
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

#RGL_IPC


Nes_id <- row.names(subset(fData(HSMM), gene_short_name == "Nes"))
Lpar1_id <- row.names(subset(fData(HSMM), gene_short_name == "Lpar1"))
Cdk1_id <- row.names(subset(fData(HSMM), gene_short_name == "Cdk1"))
Sox4_id <- row.names(subset(fData(HSMM), gene_short_name == "Sox4"))

cth <- newCellTypeHierarchy()

cth <- addCellType(cth,
                   "RG-like quiescent",
                   classify_func = function(x) { x[Nes_id,] >= 0.5 })
cth <- addCellType(cth,
                   "RG-like active",
                   classify_func = function(x) { x[Lpar1_id,] >=0.5 })
cth <- addCellType(cth,
                   "NPC",
                   classify_func = function(x) { x[Cdk1_id,] >=1 })
cth <- addCellType(cth,
                   "neuroblast",
                   classify_func = function(x) { x[Sox4_id,] >=1 })
#2. NSC lineage


Hmgn2_id <- row.names(subset(fData(HSMM), gene_short_name == "Hmgn2"))
Sox11_id <- row.names(subset(fData(HSMM), gene_short_name == "Sox11"))
Syt1_id <- row.names(subset(fData(HSMM), gene_short_name == "Syt1"))

cth <- newCellTypeHierarchy()

cth <- addCellType(cth,
                   "NSC",
                   classify_func = function(x) { x[Clu_id,] >= 4 })
cth <- addCellType(cth,
                   "early NPC",
                   classify_func = function(x) { x[Hmgn2_id,] >= 2 })
cth <- addCellType(cth,
                   "late NPC",
                   classify_func = function(x) { x[Sox11_id,] >= 2.5 })
cth <- addCellType(cth,
                   "Neuron",
                   classify_func = function(x) { x[Syt1_id,] >= 1.5 })

# Save the file as markerDiffTable is slow.......
#save (HSMM, file = file.path(project_folder, 'HSMM_OL_Lineage_0225.Robj') )
load( file = file.path(project_folder, 'HSMM_OL_Lineage_0225.Robj') )

#Run "expressed genes" and genes_id when uploading a Monocle object:
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 10))

#########
#OPC lineage
#########
Pdgfra_id <- row.names(subset(fData(HSMM), gene_short_name == "Pdgfra"))
Fyn_id <-row.names(subset(fData(HSMM), gene_short_name == "Fyn"))
Plp1_id <- row.names(subset(fData(HSMM), gene_short_name == "Plp1"))
cth <- newCellTypeHierarchy()

cth <- addCellType(cth,
                   "OPCs",
                   classify_func = function(x) { x[Pdgfra_id,] >= 1.5 })
cth <- addCellType(cth,
                   "commited OPC",
                   classify_func = function(x) { x[Fyn_id,] >= 2 })
cth <- addCellType(cth,
                   "mature OLs",
                   classify_func = function(x) { x[Plp1_id,] >= 4 })



#3. Classify

HSMM <- classifyCells(HSMM, cth)



#4. Change cores based on the CPU and system memory availible. ~4GB meme per core for ~5k cells and 100% CPU
marker_diff <- markerDiffTable(HSMM[expressed_genes,], cth, cores = 5)



#5. Add in 1,000 more genes to help Moncole order based on the supervised genes above.
semisup_clustering_genes <- row.names(marker_diff)[order(marker_diff$qval)][1:1000]
HSMM <- setOrderingFilter(HSMM, semisup_clustering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2, method = 'DDRTree', norm_method = 'log')
HSMM <- orderCells(HSMM)

# Make beautiful plots and save them!
plot_cell_trajectory(HSMM, color_by = "Pseudotime") + 
  theme(legend.position = "right") +
  ggsave(file.path(analysis_folder,"Plots", "RGL_Pseudotime.png"), height=5, width=9)
plot_cell_trajectory(HSMM, color_by = "CellType") + 
  theme(legend.position = "right") +
  ggsave(file.path(analysis_folder,"Plots", "RGL_CellType.png"), height=5, width=9)
plot_cell_trajectory(HSMM, color_by = "Clusters") + 
  theme(legend.position = "right") +
  ggsave(file.path(analysis_folder,"Plots", "RGL_Clusters.png"), height=5, width=9)

# Run plots for branched pseudotime
#HSMM_filtered <- HSMM[expressed_genes,]

#my_genes <- row.names(subset(fData(HSMM_filtered),
#                             gene_short_name %in% c("Olig1", "Olig2", "Sox10")))

#cds_subset <- HSMM_filtered[my_genes,]
#plot_genes_branched_pseudotime(cds_subset,
#                               branch_point = 0,
#                               color_by = "Clusters",
#                               ncol = 1)


#6. Run function and plot graphs. 

generate.pseudo.plot <- function( genes ) {
  #Finding Genes that Change as a Function of Pseudotime
  to_be_tested <- row.names(subset(fData(HSMM), gene_short_name %in% genes ))
  cds_subset <- HSMM[to_be_tested,]
  diff_test_res <- differentialGeneTest(cds_subset, fullModelFormulaStr = "~sm.ns(Pseudotime)")
  diff_test_res[,c("gene_short_name", "pval", "qval")]
  
  file.name <- file.path(analysis_folder,"Plots", paste(genes[1], "Pseudotime.png"))
  print(paste("Saving plot to",file.name))

  plot_genes_in_pseudotime(cds_subset, color_by = "Clusters") + 
    theme(legend.position = "top") + ggsave(file.name, height=9, width=9)
}

#OL genes:

#generate.pseudo.plot(c("ECEV4806transgene"))
#generate.pseudo.plot(c("Olig1", "Olig2", "Sox10", "Pdgfra", "Cspg4"))
#generate.pseudo.plot(c("Bmp4", "Fyn", "Gpr17"))
#generate.pseudo.plot(c("Plp1", "Mal", "Mog", "Mbp"))

#Summary_genes OL:
generate.pseudo.plot(c("Lpar1", "Nes", "Prom1"))
generate.pseudo.plot(c("Cdk1","Aurkb", 'Top2a' ))
generate.pseudo.plot(c("Dcx","Calb2", "Sox11"))

#NSC genes:
#generate.pseudo.plot(c("ECEV4806transgene", "Aldoc", "Apoe"))
#generate.pseudo.plot(c("Id4", "Hopx", "Sox9", "Gfap"))
#generate.pseudo.plot(c("Slc1a3", "Sox2", "Fabp7", "Clu"))
#generate.pseudo.plot(c("Eomes","Ccnd2","Hmgn2"))
#generate.pseudo.plot(c('Neurod1', 'Dcx', "Sox11"))
#generate.pseudo.plot(c("Calb2",'Cd24a'))
#generate.pseudo.plot(c("Calb1", "Gria1"))
#generate.pseudo.plot(c("Rbfox3","Syt1","Syt5","Myt1l"))

#Summary_genes NSC:
generate.pseudo.plot(c("Sox2", "Clu","Slc1a3", "Gfap"))
generate.pseudo.plot(c('Ccnd2', "Neurod1",'Dcx', 'Gria1'))

#New_gene_candidates:"Mia", "Dynlrb2", "Ascl1", "Kcne1l", "Igfbp5", "Ube2c", 'Cks2', "Islr2"
generate.pseudo.plot(c("Mia", "Kcne1l", "Igfbp5"))
generate.pseudo.plot(c("Ube2c", 'Cks2', "Islr2"))

#Clustering Genes by Pseudotemporal Expression Pattern
#Run top_genes from Seurat FIRST
# Use this line to get the top 25 genes used to generate the Seurat heatmaps for the specificied clusters.
# (top_genes %>% filter(cluster %in% ("RGL-activated") ))$gene


top_15 <- (top_genes %>% filter(cluster %in% c("RG-like quiescent", "RG-like active", "NPC", "neuroblast") ))$gene

marker_genes <- row.names(subset(fData(HSMM),
                                 gene_short_name %in% top_15 ))
diff_test_res <- differentialGeneTest(HSMM[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(HSMM[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
write.csv(, file = file.path(project_folder, 'Negative_Low_High_filtering_h_pseudotime_markers_30.csv') )
?plot_pseudotime_heatmap
