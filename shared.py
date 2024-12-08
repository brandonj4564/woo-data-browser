from pathlib import Path
import scanpy as sc
import pandas as pd

app_dir = Path(__file__).parent

hpf5_adata = sc.read_h5ad(app_dir / "data/5hpf_data.h5ad")
hpf5_obs_metadata = hpf5_adata.obs_keys()
hpf5_genes = list(hpf5_adata.var_names)

hpf8_adata = sc.read_h5ad(app_dir / "data/8hpf_data.h5ad")
hpf8_obs_metadata = hpf8_adata.obs_keys()
hpf8_genes = list(hpf8_adata.var_names)

hpf10_adata = sc.read_h5ad(app_dir / "data/10hpf_data.h5ad")
hpf10_obs_metadata = hpf10_adata.obs_keys()
hpf10_genes = list(hpf10_adata.var_names)

hpf12_adata = sc.read_h5ad(app_dir / "data/12hpf_data.h5ad")
hpf12_obs_metadata = hpf12_adata.obs_keys()
hpf12_genes = list(hpf12_adata.var_names)

adata_dict = {
    "5 hours past fertilization": hpf5_adata,
    "8 hours past fertilization": hpf8_adata,
    "10 hours past fertilization": hpf10_adata,
    "12 hours past fertilization": hpf12_adata
}

obs_metadata_dict = {
    "5 hours past fertilization": hpf5_obs_metadata,
    "8 hours past fertilization": hpf8_obs_metadata,
    "10 hours past fertilization": hpf10_obs_metadata,
    "12 hours past fertilization": hpf12_obs_metadata
}

genes_dict = {
    "5 hours past fertilization": hpf5_genes,
    "8 hours past fertilization": hpf8_genes,
    "10 hours past fertilization": hpf10_genes,
    "12 hours past fertilization": hpf12_genes
}