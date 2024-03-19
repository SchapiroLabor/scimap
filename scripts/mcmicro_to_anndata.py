#!/usr/bin/env python
import time
import argparse
import os
from pathlib import Path
import pandas as pd
import scimap as sm


def get_args():

    parser = argparse.ArgumentParser(description='Converts McMicro data to SciMAP format')
    parser.add_argument('--input_dir',  type=str, dest='input_path',  help='Path to the feature table')
    parser.add_argument('--output_dir', type=str, dest='output_path', help='Output directory')
    inputs.add_argument("--version", action="version", version="0.1.0")
    args = parser.parse_args()
    # Standardize paths
    args.input_path  = Path(args.input_path)
    args.output_path = Path(args.output_path)
    return args

def call_mcmicro_to_scimap(args):
    adata = sm.pp.mcmicro_to_scimap(
        feature_table_path=args.input_path,
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
        verbose=False, #not default
        output_dir=None)
    return adata

if __name__ == "__main__":
    args = get_args()
    st = time.time()
    adata = call_mcmicro_to_scimap(args)
    adata.write(args.output_path)
    rt = time.time() - st
    print(f"Script finished in {rt // 60:.0f}m {rt % 60:.0f}s")

