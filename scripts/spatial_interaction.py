# Import build-in libraries
import argparse

# Import external libraries
import anndata
import numpy as np
import scimap as sm
import pandas as pd
import plotly.express as px
import matplotlib.pyplot as plt
import seaborn as sns
import os

sns.set(color_codes=True, style="white")

# Gets arguments from the command line
def get_args():
    
    # Add description
    parser = argparse.ArgumentParser(description="Spatial_interaction for MCMICRO")

    # Add input arguments
    inputs_arguments = parser.add_argument_group(title='Input', description='Input arguments to run spatial_interaction')
    inputs_arguments.add_argument("-i", "--input", dest="input", required=True, action="store",
                                help="path to dir with files for spatialLDA it can be a h5ad file or a csv file")

    # Add module arguments
    module_arguments = parser.add_argument_group(title='SpatialInteraction',
                                                description='SpatialInteraction arguments to run SpatialInteraction')
    module_arguments.add_argument("-p", "--phenotype", dest="phenotype", required=False,
                                type=str, action="store", default="phenotype",
                                help="Phenotype column name in the input file")
    module_arguments.add_argument("-x", "--x-coordinate", dest="x", required=False,
                                type=str, action="store", default="X_centroid",
                                help="X coordinate column name in the input file")
    module_arguments.add_argument("-y", "--y-coordinate", dest="y", required=False,
                                type=str, action="store", default="Y_centroid",
                                help="Y coordinate column name in the input file")
    module_arguments.add_argument("-z", "--z-coordinate", dest="z", required=False,
                                type=str, action="store", default=None,
                                help="Z coordinate column name in the input file")
    module_arguments.add_argument("-s", "--sample-id", dest="sampleid", required=False,
                                type=str, action="store", default="sampleid",
                                help="Sampleid column name in the input file")
    module_arguments.add_argument("-sub", "--subset", dest="subset", required=False,
                                type=str, action="store", default=None,
                                help="Specific image identifier for targeted analysis.")
    module_arguments.add_argument("-m", "--method", dest="method", required=False,
                                type=str, action="store", default="radius", choices=["radius", "knn"],
                                help="Method to use in spatialLDA")
    module_arguments.add_argument("-r", "--radius", dest="radius", required=False,
                                type=int, action="store", default=30,
                                help="Radius for the neighborhood")
    module_arguments.add_argument("-k", "--knn", dest="knn", required=False,
                                type=int, action="store", default=10,
                                help="Number of neighbors to use in spatial_interaction")
    module_arguments.add_argument("-per", "--permutation", dest="permutation", required=False,
                                type=int, action="store", default=1000,
                                help="Number of permutations for p-value calculation")
    module_arguments.add_argument("-pm", "--pval_method", dest="pval_method", required=False,
                                type=str, action="store", default='zscore', choices=["abs", "zscore"],
                                help="Method for p-value calculation:")
    
    # verbose and label arguments not added
    
    # Add output arguments
    output_arguments = parser.add_argument_group(title='Output',
                                                description='Output arguments for SpatialInteraction')
    output_arguments.add_argument("-o", "--output", dest="output",
                                required=True, type=str, action="store",
                                help="Output file for SpatialInteraction, anndata")
    output_arguments.add_argument("-icp", "--interaction-composition-plot", dest="icp",
                                required=True, type=str, action="store",
                                help="Output file for the neighborhood composition plot")

    # Add version
    parser.add_argument("--version", action="version", version="1.3.3")

    # Parse arguments
    args = parser.parse_args()

    # Standardize paths
    args.input  = os.path.abspath(args.input)
    args.output = os.path.abspath(args.output)
    args.icp    = os.path.abspath(args.icp) if args.icp else None

    return args

def main():
    # Get arguments
    args = get_args()

    list_of_files = os.listdir(args.input)
    full_paths = [os.path.join(args.input, file) for file in list_of_files]

    if list_of_files[0].endswith('.csv'):

        adata = sm.pp.mcmicro_to_scimap(
            feature_table_path=full_paths,
            remove_dna=True,
            remove_string_from_name=None,
            log=False, #not default
            drop_markers=None,
            random_sample=None,
            unique_CellId=True,
            CellId='CellID',
            split='X_centroid',
            custom_imageid=None,
            min_cells=None,
            output_dir=None)

    elif list_of_files[0].endswith('.h5ad'):
        adata = anndata.read_h5ad(list_of_files[0])
    # TODO has to be tested with anndata saved file
    
    #run spatial_interaction
    sm.tl.spatial_interaction(
        adata       = adata, 
        x_coordinate= args.x, 
        y_coordinate= args.y, 
        z_coordinate= args.z, 
        phenotype   = args.phenotype, 
        method      = args.method, 
        radius      = args.radius, 
        knn         = args.knn, 
        permutation = args.permutation, 
        imageid     = args.sampleid, 
        subset      = args.subset, 
        pval_method = args.pval_method,
        label       = 'spatial_interaction'
        )

    #save the output
    adata.write(filename=args.output, compression=None, compression_opts=None, as_dense=())
    
    #save the plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    ax = sm.pl.spatial_interaction(
                        adata=adata, 
                        spatial_interaction='spatial_interaction', 
                        summarize_plot=True, 
                        p_val=0.05, 
                        row_cluster=False, 
                        col_cluster=False, 
                        cmap='coolwarm', 
                        nonsig_color='grey', 
                        subset_phenotype=None, 
                        subset_neighbour_phenotype=None, 
                        binary_view=False, 
                        return_data=False)
    
    plt.savefig(args.icp, format='png', dpi=300)

if __name__ == "__main__":
    main()


"""
Example of usage:
python spatial_interaction.py -i /path/to/input.h5ad -o /path/to/output.h5ad -icp /path/to/interaction_composition_plot.png
"""
