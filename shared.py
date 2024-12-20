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

shared_categories = set.intersection(
    *(set(adata.obs[obs_key].unique()) for adata in adata_dict.values())
)
print("Shared categories:", shared_categories)

# Ensure shared categories are consistently ordered across datasets
for adata in adata_dict.values():
    shared_order = sorted(shared_categories)
    unique_categories = [
        cat for cat in adata.obs[obs_key].unique() if cat not in shared_categories
    ]
    full_order = shared_order + unique_categories

    # Reorder categories for the obs_key
    adata.obs[obs_key] = pd.Categorical(
        adata.obs[obs_key],
        categories=full_order,
        ordered=True,
    )
