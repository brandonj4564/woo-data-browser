from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from shiny import reactive
from shiny.express import input, render, ui
from shared import adata_dict, obs_metadata_dict, genes_dict, hpf10_genes
import scanpy as sc
import io

dataSet = reactive.value("8 hours past fertilization")

# ui.panel_title("Woo Lab Dataset")
ui.page_opts(title="Woo Lab Dataset", fillable=True)
ui.include_css(Path(__file__).parent / "styles.css")


with ui.sidebar():
    ui.input_select(
        "hpf_data",
        "Select data set",
        choices=[
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


with ui.navset_pill(id="tab"):

    with ui.nav_panel("Gene UMAP"):
        with ui.card():
            # Define colors for the gradient: gray to dark blue/purple
            colors = [
                "#808080",
                "#4B0082",
            ]  # Gray (#808080) to indigo (#4B0082)

            # Create the colormap
            gray_to_purple = LinearSegmentedColormap.from_list("GrayToPurple", colors)

            @render.text
            def no_gene_selected_feature():
                if input.specific_gene() == ():
                    return "Select a gene to generate a plot."

            @render.plot
            def feature_plot():
                markers = list(input.specific_gene())

                if input.specific_gene() != ():
                    adata = adata_dict[dataSet.get()]

                    sc.pl.umap(
                        adata, color=markers, show=False, color_map=gray_to_purple
                    )

                return plt.gcf()

            with ui.layout_columns():

                @render.download(label="Download PNG", filename="feature_plot.png")
                def download_feature_png():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.umap(
                        adata, color=markers, show=False, color_map=gray_to_purple
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="png", bbox_inches="tight")
                        buf.seek(
                            0
                        )  # This forces the color scheme to be applied to the download too
                        yield buf.getvalue()

                @render.download(label="Download PDF", filename="feature_plot.pdf")
                def download_feature_pdf():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.umap(
                        adata, color=markers, show=False, color_map=gray_to_purple
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="pdf", bbox_inches="tight")
                        buf.seek(
                            0
                        )  # This forces the color scheme to be applied to the download too
                        yield buf.getvalue()

                @render.download(label="Download SVG", filename="feature_plot.svg")
                def download_feature_svg():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.umap(
                        adata, color=markers, show=False, color_map=gray_to_purple
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="svg", bbox_inches="tight")
                        buf.seek(
                            0
                        )  # This forces the color scheme to be applied to the download too
                        yield buf.getvalue()

        with ui.card():

            @render.plot()
            def umap_plot():
                color_column = input.obs_meta()
                adata = adata_dict[dataSet.get()]

                sc.pl.umap(
                    adata,
                    color=color_column,
                    show=False,
                    title=f"UMAP colored by {color_column}",
                )
                return plt.gcf()

            with ui.layout_columns():

                @render.download(label="Download PNG", filename="umap.png")
                def download_umap_png():
                    adata = adata_dict[dataSet.get()]

                    color_column = input.obs_meta()
                    adata = adata_dict[dataSet.get()]
                    sc.pl.umap(
                        adata,
                        color=color_column,
                        # save=True,  # Save the plot directly to filepath
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="png", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download PDF", filename="umap.pdf")
                def download_umap_pdf():
                    adata = adata_dict[dataSet.get()]

                    color_column = input.obs_meta()
                    adata = adata_dict[dataSet.get()]
                    sc.pl.umap(
                        adata,
                        color=color_column,
                        # save=True,  # Save the plot directly to filepath
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="pdf", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download SVG", filename="umap.svg")
                def download_umap_svg():
                    adata = adata_dict[dataSet.get()]

                    color_column = input.obs_meta()
                    adata = adata_dict[dataSet.get()]
                    sc.pl.umap(
                        adata,
                        color=color_column,
                        # save=True,  # Save the plot directly to filepath
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="svg", bbox_inches="tight")
                        yield buf.getvalue()

    with ui.nav_panel("Violin Plot"):

        with ui.card():

            @render.text
            def no_gene_selected_violin():
                if input.specific_gene() == ():
                    return "Select a gene to generate a plot."

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

            with ui.layout_columns():

                @render.download(label="Download PNG", filename="violin_plot.png")
                def download_violin_png():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.stacked_violin(
                        adata,
                        markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="png", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download PDF", filename="violin_plot.pdf")
                def download_violin_pdf():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.stacked_violin(
                        adata,
                        markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="pdf", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download SVG", filename="violin_plot.svg")
                def download_violin_svg():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.stacked_violin(
                        adata,
                        markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="svg", bbox_inches="tight")
                        yield buf.getvalue()

    with ui.nav_panel("Dot Plot"):
        with ui.card():

            @render.text
            def no_gene_selected_dot():
                if input.specific_gene() == ():
                    return "Select a gene to generate a plot."

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

                    plt.title("Dot Plot: Top Marker Genes")
                    plt.xlabel("Clusters")
                    plt.ylabel("Genes")

                return plt.gcf()

            with ui.layout_columns():

                @render.download(label="Download PNG", filename="dot_plot.png")
                def download_dot_png():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.dotplot(
                        adata,
                        var_names=markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="png", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download PDF", filename="dot_plot.pdf")
                def download_dot_pdf():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.dotplot(
                        adata,
                        var_names=markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="pdf", bbox_inches="tight")
                        yield buf.getvalue()

                @render.download(label="Download SVG", filename="dot_plot.svg")
                def download_dot_svg():
                    adata = adata_dict[dataSet.get()]
                    markers = list(input.specific_gene())

                    if input.specific_gene() == ():
                        raise Exception("Add at least one gene!")

                    sc.pl.dotplot(
                        adata,
                        var_names=markers,
                        groupby=input.obs_meta(),
                        show=False,
                    )

                    with io.BytesIO() as buf:
                        plt.savefig(buf, format="svg", bbox_inches="tight")
                        yield buf.getvalue()
