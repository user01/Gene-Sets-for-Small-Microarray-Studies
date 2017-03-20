import argparse
import math
import glob
import os

from itertools import combinations
from functools import reduce
import numpy as np
import pandas as pd


parser = argparse.ArgumentParser(
    description='Build gene sets from scores')
parser.add_argument('--seed', type=int, default=451,
                    help='Random seed for operation')
parser.add_argument('--low', type=int, default=15,
                    help='Smallest allowable size of gene set')
parser.add_argument('--high', type=int, default=200,
                    help='Largest allowable size of gene set')
parser.add_argument('--count', type=int, default=10,
                    help='Number of sets to attempt to generate')

parser.add_argument('--name', type=str, required=True,
                    help='Name of cell type to target')
parser.add_argument('--type', type=str, required=True,
                    help='General or specific cell type')

parser.add_argument('--input', type=str, required=True,
                    help='Path to input data')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output set files')

args = parser.parse_args()
# args = parser.parse_args(([
#     '--low', '16',
#     '--high', '100',
#     '--input', 'results',
#     '--seed', '0',
#     '--type', 'General_Cell_Type',
#     '--name', '"Monocyte"',
#     '--output', 'results'
# ]))

# Important since some names have spaces
cell_name = args.name.replace('"', '').replace(' ', '_')

args.type

root_filename = 'score.*.{}.{}.tsv'.format(args.type, cell_name)
root_glob = os.path.join('results', root_filename)
score_paths = glob.glob(root_glob)


def rbind(df_a, df_b):
    return pd.DataFrame.append(df_a, df_b) \
        .reset_index(drop=True)

score_files = map(pd.read_table, score_paths)
score_data = reduce(rbind, score_files)


paired_scores = score_data.assign(weighted_score=score_data.frequency * score_data.score) \
    .groupby(['low', 'high'], as_index=False) \
    .agg({'weighted_score': 'sum'}) \
    .sort_values('weighted_score', ascending=False) \
    .reset_index(drop=True)
