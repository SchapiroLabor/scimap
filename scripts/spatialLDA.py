# Import build-in libraries
from pathlib import Path
import argparse

# Import external libraries
import anndata
import numpy as np
import scimap as sm
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns;

sns.set(color_codes=True)


# Gets arguments from the command line
def get_args():
    # Add description
    parser = argparse.ArgumentParser(description="SpatialLDA for MCMICRO")

    # Add input arguments
    inputs_arguments = parser.add_argument_group(title='Input',
                                                 description='Input arguments to run spatialLDA')
    inputs_arguments.add_argument("-i", "--input", dest="input", required=True,
                                  type=str, action="store",
                                  help="Input file for spatialLDA it can be a h5ad file or a csv file")

    # Add module arguments
    module_arguments = parser.add_argument_group(title='SpatialLDA',
                                                 description='SpatialLDA arguments to run spatialLDA')
    module_arguments.add_argument("-p", "--phenotype", dest="phenotype", required=False,
                                  type=str, action="store", default="phenotype",
                                  help="Phenotype column name in the input file")
    module_arguments.add_argument("-x", "--x-coordinate", dest="x", required=False,
                                  type=str, action="store", default="X_centroid",
                                  help="X coordinate column name in the input file")
    module_arguments.add_argument("-y", "--y-coordinate", dest="y", required=False,
                                  type=str, action="store", default="Y_centroid",
                                  help="Y coordinate column name in the input file")
    module_arguments.add_argument("-s", "--sample-id", dest="sampleid", required=False,
                                  type=str, action="store", default="sampleid",
                                  help="Sampleid column name in the input file")
    module_arguments.add_argument("-r", "--radius", dest="radius", required=False,
                                  type=int, action="store", default=30,
                                  help="Radius for the neighborhood")
    module_arguments.add_argument("-n", "--num-motifs", dest="num_motifs", required=False,
                                  type=int, action="store", default=10,
                                  help="Number of motifs to use in spatialLDA")
    module_arguments.add_argument("-k", "--knn", dest="knn", required=False,
                                  type=int, action="store", default=10,
                                  help="Number of neighbors to use in spatialLDA")
    module_arguments.add_argument("-m", "--method", dest="method", required=False,
                                  type=str, action="store", default="radius", choices=["radius", "knn"],
                                  help="Method to use in spatialLDA")

    # Add output arguments
    output_arguments = parser.add_argument_group(title='Output',
                                                 description='Output arguments for spatialLDA')
    output_arguments.add_argument("-o", "--output", dest="output",
                                  required=True, type=str, action="store",
                                  help="Output file for spatialLDA")
    output_arguments.add_argument("-ncp", "--neighborhood-composition-plot", dest="ncp",
                                  required=True, type=str, action="store",
                                  help="Output file for the neighborhood composition plot")
    output_arguments.add_argument("-mlp", "--motif-locations-plot", dest="motif_locations_plot",
                                  required=True, type=str, action="store",
                                  help="Output file for the motif locations plot")

    # Add version
    parser.add_argument("-v", "--version", action="version", version="1.3.3")

    # Parse arguments
    args = parser.parse_args()

    # Standardize paths
    args.input = Path(args.input).resolve()
    args.output = Path(args.output).resolve()
    args.ncp = Path(args.ncp).resolve() if args.ncp else None
    args.motif_locations_plot = Path(args.motif_locations_plot).resolve() if args.motif_locations_plot else None

    return args


## Plotting:
def plot_motifs_CT_composition(adata, ax=None, **kwargs):
    # create empty Axis object
    plt.gca()
    ax = ax

    # create a dataframe with the cell type composition of each cellular neighborhood
    df = pd.DataFrame(adata.obs.groupby('phenotype')['CN'].value_counts()).unstack()
    df = df.fillna(0)
    df = df.astype(int)
    df = df / df.sum(axis=0)

    # create stacked bar plot, where the x axis is the CN, and the y axis is the percentage of cells of each cluster
    # do I pick figure size already?
    df.T.plot(kind='bar', stacked=True, log=False, ax=ax, cmap='tab20')
    ax.get_yaxis().set_visible(False)
    ax.set_ylabel('Cellular composition fraction')
    ax.set_xlabel('Cellular neighborhoods')
    ax.set_title('Cellular composition of cellular neighborhoods')

    #   -------legend---------
    # get box position
    box = ax.get_position()
    # reduce by 20% to make room for legend
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cell clusters')
    #   -------legend---------

    # Annotate the bars with the fractions
    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        if height > 0.015:
            ax.text(x + width / 2,
                    y + height / 2,
                    '{:.0f}%'.format(height * 100),
                    horizontalalignment='center',
                    verticalalignment='center')

    return ax


def plot_sample_motifs(adata, ax=None, **kwargs):
    # create empty Axis object
    plt.gca()
    ax = ax

    df = pd.DataFrame(adata.obs.CN.value_counts())
    df = df.fillna(0)
    df = df.astype(int)
    df = df / df.sum(axis=0)

    df.T.plot(kind='bar', stacked=True, colormap='tab10', log=False, ax=ax)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    ax.set_ylabel('Cellular neighborhoods of tissue')

    #   -------legend---------
    # get box position
    box = ax.get_position()
    # reduce by 20% to make room for legend
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), title='Cellular \nneighborhoods')
    #   -------legend---------

    # annonate the bars with the fractions
    for p in ax.patches:
        width, height = p.get_width(), p.get_height()
        x, y = p.get_xy()
        if height > 0.015:
            ax.text(x + width / 2,
                    y + height / 2,
                    '{:.0f}%'.format(height * 100),
                    horizontalalignment='center',
                    verticalalignment='center')

    return ax

def plot_neighborhood_composition(adata, output):
    # Plot neighborhood composition and cell type composition of the cellular neighborhoods
    # NC= Neighborhood Composition
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True,
                                           gridspec_kw={'width_ratios': [1, 4]}, figsize=(20, 7))

    # Plot the sample motifs
    plot_sample_motifs(adata, ax=ax1)

    # Plot the cell type composition of the motifs
    plot_motifs_CT_composition(adata, ax=ax2)

    # Save the neighborhood composition plot
    fig.savefig(output, dpi=300, bbox_inches='tight')


def plot_motifs_locations(adata, x, y, output):
    # TODO: Prettify
    # TODO: Consider using datashader to handle many points
    # Create a dataframe with the motif locations
    df = pd.DataFrame({"X": adata.obs[x], "Y": adata.obs[y], "CN": adata.obs["CN"]})

    # Create a scatter plot with the motif locations
    fig = px.scatter(df, x="X", y="Y", color="CN", title="Motif locations")

    # Update the marker size
    fig.update_traces(marker=dict(size=15), selector=dict(mode='markers'))

    # Save the plot
    fig.write_html(output)


def main(args):
    # Read input file
    data = pd.read_csv(args.input)

    # Create a random matrix of marker intensities
    data['A'] = np.random.randint(1, 101, size=len(data))
    data['B'] = np.random.randint(1, 101, size=len(data))

    # Get the celltype columns
    celltype_columns = [args.sampleid, args.x, args.y, args.phenotype]

    # Create an AnnData object out of the marker intensities
    # Create an AnnData object out of the marker intensities
    adata = anndata.AnnData(X=data[['A', 'B']].values, obs=data[celltype_columns])

    # Run spatialLDA
    sm.tl.spatial_lda(adata,
                      x_coordinate=args.x,
                      y_coordinate=args.y,
                      phenotype=args.phenotype,
                      method=args.method,
                      radius=args.radius,
                      knn=args.knn,
                      imageid=args.sampleid,
                      num_motifs=args.num_motifs,
                      random_state=42,
                      subset=None,
                      label=f'spatialLDA_{args.num_motifs}_{args.method}_{args.radius}')

    # Add unstructured motifs to structure observations, take the most likely motif
    # CN = Cellular Neighborhoods
    adata.obs['CN'] = adata.uns[f'spatialLDA_{args.num_motifs}_{args.method}_{args.radius}'].idxmax(axis=1)

    # Plot neighborhood composition and cell type composition of the cellular neighborhoods
    plot_neighborhood_composition(adata, args.ncp)

    # Plot motif locations
    plot_motifs_locations(adata, args.x, args.y, args.motif_locations_plot)

    # Save motifs to output csv
    adata.obs.to_csv(args.output, index=False)


if __name__ == "__main__":
    args = get_args()
    main(args)