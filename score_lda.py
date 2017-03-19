import argparse
import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis


parser = argparse.ArgumentParser(description='Perform LDA ')
parser.add_argument('--seed', type=int, default=451,
                    help='Random seed for operation')
parser.add_argument('--pairs', type=int, default=200,
                    help='Number of gene pairs to pick as relevant from a component')
parser.add_argument('--bootstrap', type=float, default=0.66,
                    help='Fraction of genes to bootstrap')
parser.add_argument('--name', type=str, required=True,
                    help='Name of cell type to target')
parser.add_argument('--type', type=str, required=True,
                    help='General or specific cell type')

parser.add_argument('--input', type=str, required=True,
                    help='Path to input data')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output pair/scores data')

args = parser.parse_args()
# args = parser.parse_args(([
#     '--name', '450',
#     '--input', 'results/gene_data_vs_cell_type.tsv',
#     '--seed', '0',
#     '--type', 'General_Cell_Type',
#     '--name', "Monocyte",
#     '--output', 'results/score.00000.General_Cell_Type.Monocyte.tsv'
# ]))


args.input

data_raw = pd.read_table(args.input)
data = data_raw.assign(truth=data_raw.Cell_Type if args.type ==
                       'Cell_Type' else data_raw.General_Cell_Type) \
    .drop(['GSM_ID', 'Cell_Type', 'General_Cell_Type'], axis=1)

print(args.seed)
print(data.shape)
