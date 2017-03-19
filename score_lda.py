import argparse
import math
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


data_00 = pd.read_table(args.input)
data_01 = data_00.assign(truth=data_00.Cell_Type if args.type ==
                         'Cell_Type' else data_00.General_Cell_Type) \
    .drop(['GSM_ID', 'Cell_Type', 'General_Cell_Type'], axis=1)

others = (data_01.truth == args.name).apply(
    lambda x: args.name if x else "other")
data_02 = data.assign(truth=others)

data_normal = data_02.query('truth != "other"').reset_index(drop=True)
data_other = data_02.query('truth == "other"')\
    .sample(n=data_normal.shape[0], random_state=args.seed)\
    .reset_index(drop=True)

data_03 = pd.DataFrame.append(data_normal, data_other)\
    .sample(frac=1, random_state=args.seed + 1)\
    .reset_index(drop=True)


def lda_model(data):
    X_train = ignore_one(
        data.drop(['truth'], axis=1), idx).as_matrix()
    y_train = ignore_one(data[['truth']], idx).as_matrix().T[0]

    clf = LinearDiscriminantAnalysis()
    clf.fit(X_train, y_train)
    return clf


lda_model_main = lda_model(data_03)

lda_model_main.scalings_.shape
lda_model_main.scalings_[:, 0].shape

data_03.columns[:-1]
data_03.columns[:-1].shape

top_genes = math.floor(1 / 2 * (math.sqrt(8 * args.pairs + 1) + 1))

top_gene_names = pd.DataFrame({
    'gene': data_03.columns[:-1],
    'loading': np.abs(lda_model_main.scalings_[:, 0])})\
    .sort_values('loading', ascending=False)\
    .reset_index(drop=True)\
    .iloc[:top_genes]


print(args.seed)
print(data.shape)
