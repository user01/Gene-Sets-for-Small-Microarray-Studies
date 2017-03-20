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

set_minimal = paired_scores.head(args.low)
set_maximal = paired_scores.head(args.high)

# set_maximal.shape
# set_minimal.shape

set_range = set_maximal.shape[0] - set_minimal.shape[0]

set_step = set_range // args.count
set_step = set_step if set_step > 0 else 1  # ensure reasonable step size
head_sizes = list(np.arange(args.count) * set_step + args.low)


def set_contained_in_sets(pair, gene_sets):
    """Pair is contained in any of the gene sets"""
    for gene_set in gene_sets:
        if set_contained_in_set(pair, gene_set):
            return True
    return False


def set_contained_in_set(pair, gene_set):
    """Either of the pair exist in this gene set"""
    return len(pair & gene_set) > 0


def add_pair_to_sets(pair, gene_sets):
    """If a pair is contained in a gene set, union the sets"""
    gene_sets_new = []
    for gene_set in gene_sets:
        if set_contained_in_set(pair, gene_set):
            gene_sets_new.append(gene_set | pair)
    return gene_sets_new


def collapse_sets(gene_sets):
    """Collapse gene sets that can be unioned"""
    return collapse_sets_(gene_sets, [])


def collapse_sets_(unfinished_sets, known_sets):
    """Private - Collapse gene sets that can be unioned"""
    if len(unfinished_sets) < 1:
        return known_sets
    gene_set, *remaining_sets = unfinished_sets

    if set_contained_in_sets(gene_set, remaining_sets):
        # then this set collides with at least one of the
        # remaining_sets
        known_clean = known_sets
        untested_sets = list(map(
            lambda gs: gs | gene_set if set_contained_in_set(gs, gene_set) else gs))
    else:
    # this set stands alone
        known_clean = known_sets + [gene_set]
        untested_sets = remaining_sets

    return collapse_sets_(untested_sets, known_clean)


def pairs_to_sets(pairs):
    """Turn a low/high pandas dataframe into a list of sets"""
    return pairs_to_sets_(pairs[['low','high']], [])


def pairs_to_sets_(pairs, gene_sets):
    """Private - turn a pandas df into a list of sets"""
    if pairs.shape[0] < 1:
        return gene_sets

    pair = set([pairs.iloc[0][0], pairs.iloc[0][1]])

    if set_contained_in_sets(pair, gene_sets):
        gene_sets = add_pair_to_sets(pair, gene_sets)
    else:
        gene_sets.append(pair)

    gene_sets = collapse_sets(gene_sets)

    pairs_tail = pairs.tail(-1)
    return pairs_to_sets_(pairs_tail, gene_sets)
