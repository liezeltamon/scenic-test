renv::load()

BiocManager::version()

# Required packages

BiocManager::install(c("AUCell", "RcisTarget"))
BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost

# Optional packages

# Also install optional packages when packages have been upgraded to BioC 3.18
## Optional (but highly recommended):
# To score the network on cells (i.e. run AUCell):
BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# For various visualizations and perform t-SNEs:
BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# To support paralell execution (not available in Windows):
BiocManager::install(c("doMC", "doRNG"))
# To export/visualize in http://scope.aertslab.org
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)

# Install SCENIC
devtools::install_github("aertslab/SCENIC")

# Install Seurat and SCE
devtools::install_version("Matrix", "1.6.5")
renv::install("Seurat")
BiocManager::install("SingleCellExperiment")

# cisTargets database files

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather.sha1sum.txt",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/old/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather.sha1sum.txt")

dir.create(file.path("data-raw", "cisTarget_databases"))
setwd(file.path("data-raw", "cisTarget_databases"))

for (featherURL in dbFiles) {
  download.file(featherURL, destfile = basename(featherURL)) # saved in current dir
}

# Reset working directory to project directory
setwd(rprojroot::find_rstudio_root_file())
