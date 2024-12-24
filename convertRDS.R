library(Seurat)
library(SeuratDisk)

wd <- getwd()
filepath <- paste(wd, "/integrated_data.rds", sep = "")

# Load the Seurat object from the RDS file
seurat_obj <- readRDS(filepath)

# Save the Seurat object in H5Seurat format
SaveH5Seurat(seurat_obj, filename = "integrated_data.h5Seurat", overwrite = TRUE)

# Convert H5Seurat to H5AD format
Convert("integrated_data.h5Seurat", dest = "h5ad", overwrite = TRUE)
