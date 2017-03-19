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
    X_train = data.drop(['truth'], axis=1).as_matrix()
    y_train = data[['truth']].as_matrix().T[0]

    clf = LinearDiscriminantAnalysis()
    clf.fit(X_train, y_train)
    return clf


lda_model_main = lda_model(data_03)
top_genes = math.floor(1 / 2 * (math.sqrt(8 * args.pairs + 1) + 1))

top_gene_names = pd.DataFrame({
    'gene': data_03.columns[:-1],
    'loading': np.abs(lda_model_main.scalings_[:, 0])})\
    .sort_values('loading', ascending=False)\
    .reset_index(drop=True)\
    .iloc[:top_genes]

# np.array(top_gene_names.gene)

def ignore_one(df, idx):
    skip_list = [i for i in range(0, df.shape[0]) if i != idx]
    return df.ix[skip_list]


def keep_one(df, idx):
    skip_list = [i for i in range(0, df.shape[0]) if i == idx]
    return df.ix[skip_list]


def leave_one_out(df, idx):
    data_train = ignore_one(df, idx)
    data_test = keep_one(df, idx)
    clf = lda_model(data_train)
    X_test = data_test.drop(['truth'], axis=1).as_matrix()
    y_test = data_test[['truth']].as_matrix().T[0]
    result = clf.predict(X_test)
    return result[0] == y_test[0]

def run(data):
    return list(map(lambda idx: leave_one_out(data, idx),
                    range(0, data.shape[0])))

results = run(data_03)
score = np.sum(results) / len(results)





print(args.seed)
print(data.shape)
