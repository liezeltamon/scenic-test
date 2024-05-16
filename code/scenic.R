library(qs)
library(scCustomize)
library(SCENIC)
library(Seurat)
library(SingleCellExperiment)

src_path = file.path("data", "scrnaseq", "merged_firstANDnsco_seu_integrated.qs")
clustering_path = file.path("data", "scrnaseq", "exp_sample_idZZsex.harmony_merged_firstANDnsco_noscSHC_cluster_annotation.rds")

cisdbdata_dir = file.path("data-raw", "cisTarget_databases")
out_dir = file.path("results", "scenic")
  
subset_timepoint = "15"
subset_genotype = "WT"
subset_sample_id = "BIHi005-A"
sample_n_cells = 200
sample_seed = 290

org = "hgnc"
dataset_title = paste0(subset_genotype, "_tp", subset_timepoint, "_sample_id", subset_sample_id)

# Load and prepare data

## Load

seu_merged <- qread(src_path)
cluster_df <- readRDS(clustering_path)
seu_merged <- AddMetaData(seu_merged, cluster_df)

is_subset <- seu_merged$timepoint %in% subset_timepoint & 
  seu_merged$genotype %in% subset_genotype &
  seu_merged$sample_id %in% subset_sample_id
sc_obj <- subset(seu_merged, cells = which(is_subset))

## Prepare 

# **Can provide normalised counts as input but just use raw counts for gene filtering step**
expr_mat <- LayerData(sc_obj, layer = "counts", assay = "RNA")
expr_mat <- as.matrix(expr_mat)

# To make dataset smaller
expr_mat <- expr_mat[VariableFeatures(seu_merged), ]
set.seed(sample_seed)
expr_mat <- expr_mat[,sample(1:ncol(expr_mat), size = sample_n_cells, replace = FALSE)]
dim(expr_mat)

#DimPlot_scCustom(sc_obj, group.by = paste0("leiden_", seq(0.1, 1, 0.1))) # To eyeball which cluster to use
cellinfo_df <- as.data.frame(sc_obj[["leiden_0.1"]])
saveRDS(cellinfo_df, file="int/cellInfo.Rds")

# Running scenic - refer to http://htmlpreview.github.io/?https://github.com/aertslab/SCENIC/blob/master/inst/doc/SCENIC_Running.html, https://www.youtube.com/watch?v=L9r3KP5w1yw

## Intialise scenic

#dbs <- defaultDbNames[[org]] # Gives hg19 dataset names
dbs <- c(`500bp` = "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
         `10kb` = "hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather")
scenic_options <- initializeScenic(org = org, dbDir = cisdbdata_dir, nCores = 10, dbs = dbs, datasetTitle = dataset_title)

motifAnnotations_hgnc <- motifAnnotations # Then rerun initializeScenic()

saveRDS(scenic_options, file.path(out_dir, "scenic_options.rds"))

## Identify co-expression network

### 1. Gene filter / selection

genes_kept <- geneFiltering(expr_mat, scenic_options,
                            # Defaults:
                            minCountsPerGene = 3 * 0.01 * ncol(expr_mat),
                            minSamples = ncol(expr_mat) * 0.01)
# Check whether interesting genes kept
#c("FOXG1", "GLI3", "NEUROD6") %in% genes_kept
expr_mat_filtered <- expr_mat[genes_kept, ]
dim(expr_mat_filtered)

### 2. Correlation on input expression matrix

runCorrelation(expr_mat_filtered, scenic_options)

### 3. Run GENIE3 to infer potential TF targets

# Optional: Add log if data is not log or normalised already
expr_mat_filtered_log <- log2(expr_mat_filtered + 1) 
runGenie3(expr_mat_filtered_log, scenic_options)

## Build and score the GRN

# From tutorial (link above)
# exprMat_log <- log2(exprMat+1)
# scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings

### 1. Get co-expression modules

scenic_options <- runSCENIC_1_coexNetwork2modules(scenic_options)

### 2. Get regulons
scenic_options <- runSCENIC_2_createRegulons(scenic_options, coexMethod = c("top5perTarget")) # Toy run settings

## Identify cell states

### 3. Score GRN (regulons) in the cells with AUCell

expr_mat_log <- log2(expr_mat + 1)
scenic_options <- runSCENIC_3_scoreCells(scenic_options, expr_mat_log)

### Note - there is a an option to binarize activity - see tutorial

saveRDS(scenic_options, file.path(out_dir, "scenic_options.rds"))

### 4. (Optional) Binarize regulon activity


aucellApp <- plotTsne_AUCellApp(scenic_options, expr_mat_log)
savedSelections <- shiny::runApp(aucellApp)

##

n_pcs <- 5
scenic_options@settings$seed <- 290
file_names <- tsneAUC(scenic_options, aucType = "AUC", nPcs = n_pcs, perpl = c(5, 15, 50))
file_names

par(mfrow = c(2,3))

cellInfo <- cellinfo_df
plotTsne_compareSettings("int/tSNE_AUC_05pcs_05perpl.Rds", scenic_options, showLegend = FALSE, varName = "leiden_0.1")

getAUC(scenic_options)
