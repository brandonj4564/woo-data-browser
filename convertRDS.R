library(Seurat)
library(SeuratDisk)

wd <- getwd()
filepath = paste(wd, "/hpf_12_qc_dr_cl_dge_annot5.rds", sep="")

# Load the Seurat object from the RDS file
seurat_obj <- readRDS(filepath)

# Save the Seurat object in H5Seurat format
SaveH5Seurat(seurat_obj, filename = "12hpf_data.h5Seurat", overwrite = TRUE)

# Convert H5Seurat to H5AD format
Convert("12hpf_data.h5Seurat", dest = "h5ad", overwrite = TRUE)
