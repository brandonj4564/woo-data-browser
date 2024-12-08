from pathlib import Path
import matplotlib.pyplot as plt
from shiny import reactive
from shiny.express import input, render, ui
from shared import adata_dict, obs_metadata_dict, genes_dict, hpf10_genes
import scanpy as sc
import numpy as np
import io

app_dir = Path(__file__).parent
ui.include_css(app_dir / "styles.css")

dataSet = reactive.value("8 hours past fertilization")
ui.page_opts(title="Woo Lab Dataset", fillable=True)

# ui.panel_title("Woo Lab Dataset")

with ui.sidebar():
    ui.input_select(
        "hpf_data",
        "Select data set",
        choices=[
            "5 hours past fertilization",
            "8 hours past fertilization",
            "10 hours past fertilization",
            "12 hours past fertilization",
        ],
        selected="8 hours past fertilization",
    )

    @reactive.effect
    @reactive.event(input.hpf_data)
    def _():
        # This updates the dataSet var when the select menu changes
        dataSet.set(input.hpf_data())
        
    @render.ui
    def obs_meta_ui():
        obs_metadata = obs_metadata_dict[dataSet.get()]
        return ui.input_select(
            "obs_meta",
            "Select gene metadata",
            choices=obs_metadata,
            selected="celltypes",
        )

    @render.ui
    def gene_list_ui():
        gene_list = genes_dict[dataSet.get()]
        return ui.input_selectize(
            "specific_gene",
            "Select gene",
            choices=gene_list,
            multiple=True,
            selected="sox32",
        )
        
    @render.download(label="Download PNG", filename="figure.png")
    def download_png():
        tab = input.tab()
        adata = adata_dict[dataSet.get()]
        markers = list(input.specific_gene())
        
        if tab == "Gene UMAP":
            color_column = input.obs_meta()
            sc.pl.umap(
                adata,
                color=color_column,
                show=False,
                title=f"hpf5: UMAP colored by {color_column}",
            )
            
        elif tab == "Violin Plot":
            if input.specific_gene() != ():
                sc.pl.stacked_violin(
                    adata,
                    markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Dot Plot":
            if input.specific_gene() != ():
                sc.pl.dotplot(
                    adata,
                    var_names=markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Feature Plot":
            if input.specific_gene() != ():
                sc.pl.umap(
                    adata,
                    color=markers,
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
            
        with io.BytesIO() as buf:
            plt.savefig(buf, format="png", bbox_inches="tight")
            yield buf.getvalue()

    @render.download(label="Download JPEG", filename="figure.jpeg")
    def download_jpeg():
        tab = input.tab()
        adata = adata_dict[dataSet.get()]
        markers = list(input.specific_gene())
        
        if tab == "Gene UMAP":
            color_column = input.obs_meta()
            sc.pl.umap(
                adata,
                color=color_column,
                show=False,
                title=f"hpf5: UMAP colored by {color_column}",
            )
            
        elif tab == "Violin Plot":
            if input.specific_gene() != ():
                sc.pl.stacked_violin(
                    adata,
                    markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Dot Plot":
            if input.specific_gene() != ():
                sc.pl.dotplot(
                    adata,
                    var_names=markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Feature Plot":
            if input.specific_gene() != ():
                sc.pl.umap(
                    adata,
                    color=markers,
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
            
        with io.BytesIO() as buf:
            plt.savefig(buf, format="jpeg", bbox_inches="tight")
            yield buf.getvalue()

    @render.download(label="Download SVG", filename="figure.svg")
    def download_svg():
        tab = input.tab()
        adata = adata_dict[dataSet.get()]
        markers = list(input.specific_gene())
        
        if tab == "Gene UMAP":
            color_column = input.obs_meta()
            sc.pl.umap(
                adata,
                color=color_column,
                show=False,
                title=f"hpf5: UMAP colored by {color_column}",
            )
            
        elif tab == "Violin Plot":
            if input.specific_gene() != ():
                sc.pl.stacked_violin(
                    adata,
                    markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Dot Plot":
            if input.specific_gene() != ():
                sc.pl.dotplot(
                    adata,
                    var_names=markers,
                    groupby=input.obs_meta(),
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
        
        elif tab == "Feature Plot":
            if input.specific_gene() != ():
                sc.pl.umap(
                    adata,
                    color=markers,
                    show=False,
                )
            else:
                yield "No gene selected"
                raise Exception("No gene selected")
            
        with io.BytesIO() as buf:
            plt.savefig(buf, format="svg", bbox_inches="tight")
            yield buf.getvalue()


with ui.navset_pill(id="tab"):

    with ui.nav_panel("Gene UMAP"):
        
        with ui.card():
            @render.plot
            def umap_plot():
                color_column = input.obs_meta()
                adata = adata_dict[dataSet.get()]
                sc.pl.umap(
                    adata,
                    color=color_column,
                    show=False,
                    title=f"hpf5: UMAP colored by {color_column}",
                )
                return plt.gcf()

    with ui.nav_panel("Violin Plot"):

        with ui.card():
            @render.plot
            def violin_plot():

                # sc.pl.violin(
                #     adata,
                #     keys="sox32",  # Gene or metadata field
                #     groupby="celltypes",  # Categorical variable for grouping
                #     jitter=True,  # Add jitter to points
                #     scale="width",
                #     show=True,
                # )

                markers = list(input.specific_gene())
                if input.specific_gene() != ():
                    adata = adata_dict[dataSet.get()]
                    sc.pl.stacked_violin(
                        adata,
                        markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    return plt.gcf()

    with ui.nav_panel("Dot Plot"):
        
        with ui.card():
            @render.plot
            def dot_plot():
                markers = list(input.specific_gene())
                if input.specific_gene() != ():
                    adata = adata_dict[dataSet.get()]
                    sc.pl.dotplot(
                        adata,
                        var_names=markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    # plt.title("Dot Plot: Top Marker Genes")
                    # plt.xlabel("Clusters")
                    # plt.ylabel("Genes")

                    return plt.gcf()

    with ui.nav_panel("Feature Plot"):

        with ui.card():
            @render.plot
            def feature_plot():
                markers = list(input.specific_gene())

                if len(markers) > 0:
                    adata = adata_dict[dataSet.get()]
                    sc.pl.umap(
                        adata,
                        color=markers,
                        show=False,
                    )
                    return plt.gcf()