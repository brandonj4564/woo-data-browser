from pathlib import Path
import scanpy as sc
import pandas as pd

app_dir = Path(__file__).parent

hpf8_adata = sc.read_h5ad(app_dir / "data/8hpf_data_float16.h5ad")
hpf8_obs_metadata = hpf8_adata.obs_keys()
hpf8_genes = list(hpf8_adata.var_names)

hpf10_adata = sc.read_h5ad(app_dir / "data/10hpf_data_float16.h5ad")
hpf10_obs_metadata = hpf10_adata.obs_keys()
hpf10_genes = list(hpf10_adata.var_names)

hpf12_adata = sc.read_h5ad(app_dir / "data/12hpf_data_float16.h5ad")
hpf12_obs_metadata = hpf12_adata.obs_keys()
hpf12_genes = list(hpf12_adata.var_names)

adata_dict = {
    "8 hours past fertilization": hpf8_adata,
    "10 hours past fertilization": hpf10_adata,
    "12 hours past fertilization": hpf12_adata,
}

obs_metadata_dict = {
    "8 hours past fertilization": hpf8_obs_metadata,
    "10 hours past fertilization": hpf10_obs_metadata,
    "12 hours past fertilization": hpf12_obs_metadata,
}

genes_dict = {
    "8 hours past fertilization": hpf8_genes,
    "10 hours past fertilization": hpf10_genes,
    "12 hours past fertilization": hpf12_genes,
}

# This section changes the obs keys of all datasets to be in the same order to have the same color scheme
# Define datasets and obs_key
obs_key = "celltypes"  # Replace with your actual obs_key

# Collect all unique categories across datasets
all_categories = set()
for adata in adata_dict.values():
    all_categories.update(adata.obs[obs_key].unique())

# Define the common order for the categories
common_order = sorted(all_categories)  # Adjust the sorting if needed

# Apply the common order to all datasets
for adata in adata_dict.values():
    adata.obs[obs_key] = pd.Categorical(
        adata.obs[obs_key], categories=common_order, ordered=True
    )
