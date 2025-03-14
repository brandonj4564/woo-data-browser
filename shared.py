from pathlib import Path
import scanpy as sc
import pandas as pd

app_dir = Path(__file__).parent


# Renames certain columns to be more understandable
def renameColumns(adata):
    if adata.obs["Cell Types"].dtype.name == "category":
        adata.obs["Cell Types"] = adata.obs["Cell Types"].cat.rename_categories(
            {"Neural Plate": "Neural Ectoderm", "Adaxial": "Adaxial Mesoderm"}
        )
    else:
        # Fallback for non-categorical columns
        adata.obs["Cell Types"] = adata.obs["Cell Types"].replace(
            "Neural Plate", "Neural Ectoderm"
        )
        adata.obs["Cell Types"] = adata.obs["Cell Types"].replace(
            "Adaxial", "Adaxial Mesoderm"
        )


# Assigns each category in Cell Types a specific color
def createCellTypeColors(adata, colors):
    # Ensure the column is categorical
    if adata.obs["Cell Types"].dtype.name != "category":
        adata.obs["Cell Types"] = adata.obs["Cell Types"].astype("category")

    # Reorder categories to match the color map keys
    adata.obs["Cell Types"] = adata.obs["Cell Types"].cat.set_categories(colors.keys())

    # Assign colors to the groups
    adata.uns["Cell Types_colors"] = [
        colors[group] for group in adata.obs["Cell Types"].cat.categories
    ]


hpf8_adata = sc.read_h5ad(app_dir / "data/8hpf_data.h5ad")
hpf8_obs_metadata = hpf8_adata.obs_keys()
hpf8_genes = list(hpf8_adata.var_names)

hpf10_adata = sc.read_h5ad(app_dir / "data/10hpf_data.h5ad")
hpf10_obs_metadata = hpf10_adata.obs_keys()
hpf10_genes = list(hpf10_adata.var_names)

hpf12_adata = sc.read_h5ad(app_dir / "data/12hpf_data.h5ad")
hpf12_obs_metadata = hpf12_adata.obs_keys()
hpf12_genes = list(hpf12_adata.var_names)

integrated_adata = sc.read_h5ad(app_dir / "data/integrated_data.h5ad")
integrated_obs_metadata = integrated_adata.obs_keys()
integrated_genes = list(integrated_adata.var_names)

renameColumns(hpf8_adata)
renameColumns(hpf10_adata)
renameColumns(hpf12_adata)
renameColumns(integrated_adata)

adata_dict = {
    "8 hours past fertilization": hpf8_adata,
    "10 hours past fertilization": hpf10_adata,
    "12 hours past fertilization": hpf12_adata,
    "Integrated data": integrated_adata,
}

obs_metadata_dict = {
    "8 hours past fertilization": hpf8_obs_metadata,
    "10 hours past fertilization": hpf10_obs_metadata,
    "12 hours past fertilization": hpf12_obs_metadata,
    "Integrated data": integrated_obs_metadata,
}

genes_dict = {
    "8 hours past fertilization": hpf8_genes,
    "10 hours past fertilization": hpf10_genes,
    "12 hours past fertilization": hpf12_genes,
    "Integrated data": integrated_genes,
}

# Define a set color scheme for all categories in Cell Types
group_colors = {
    "Endoderm": "#7CC46D",
    "Notochord": "#F298BA",
    "Prechordal Plate": "#F69534",
    "Axial Mesoderm": "#F8CCCE",
    "Adaxial Mesoderm": "#C29AC7",
    "Neural Ectoderm": "#C9CCB4",
    "Non-neural Ectoderm": "#D0A47A",
    "DFC": "#8ABBE5",
    "LPM": "#F49A98",
    "PSM": "#E2CCE3",
    "YSL": "#EDE246",
    "EVL": "#B2B0CA",
    "Hypochord": "#0073B3",
}

createCellTypeColors(hpf8_adata, group_colors)
createCellTypeColors(hpf10_adata, group_colors)
createCellTypeColors(hpf12_adata, group_colors)
createCellTypeColors(integrated_adata, group_colors)


# # This section changes the obs keys of all datasets to be in the same order to have the same color scheme
# # Define datasets and obs_key
obs_key = "Cell Types" 

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
