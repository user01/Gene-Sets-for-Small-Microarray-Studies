import argparse
import math
import glob
import os
import sys

# Note that python doesn't support tail call optimization, so the recursive
#  calls for set creation should be optimized for a loop
sys.setrecursionlimit(2000000)

from itertools import combinations
import numpy as np
import pandas as pd

from utilities.common import rbind_all, pairs_to_sets


parser = argparse.ArgumentParser(
    description='Build gene sets from scores')
parser.add_argument('--low', type=int, default=15,
                    help='Smallest allowable size of gene set')
parser.add_argument('--high', type=int, default=200,
                    help='Largest allowable size of gene set')
parser.add_argument('--count', type=int, default=10,
                    help='Number of sets to attempt to generate')

parser.add_argument('--title', type=str, required=True,
                    help='Title string')
parser.add_argument('--name', type=str, required=True,
                    help='Name of cell type to target')
parser.add_argument('--type', type=str, required=True,
                    help='General or specific cell type')

parser.add_argument('--input', type=str, required=True,
                    help='Path to input score files')
parser.add_argument('--output', type=str, required=True,
                    help='Path to output set id files')
parser.add_argument('--feedback', type=str, required=True,
                    help='Path to feedback file (list of created sets)')

args = parser.parse_args()
# args = parser.parse_args(([
#     '--low', '100',
#     '--high', '200',
#     '--count', '300',
#     '--input', 'results/100_200_main',
#     '--title', '100_200_main',
#     '--type', 'General_Cell_Type',
#     '--name', '"B cell"',
#     '--output', 'results/100_200_main',
#     '--feedback', 'results/feedback.sets.tsv'
# ]))


# Important since some names have spaces
cell_name = args.name.replace('"', '').replace(' ', '_')

root_filename = '{}.score.*.{}.{}.tsv'.format(
    args.title, args.type, cell_name)
root_glob = os.path.join(args.input, root_filename)
score_paths = glob.glob(root_glob)

def temp_read(path):
    try:
        return pd.read_table(path)
    except:
        print(path)
        pass

# score_files = list(map(pd.read_table, score_paths))
score_files = list(map(temp_read, score_paths))
score_data = rbind_all(score_files)


paired_scores = score_data.assign(weighted_score=score_data.frequency * score_data.score) \
    .groupby(['low', 'high'], as_index=False) \
    .agg({'weighted_score': 'sum'}) \
    .sort_values('weighted_score', ascending=False) \
    .reset_index(drop=True)


def set_values(paired_scores, min_size, max_size, head_size):
    """For top head_size paired scores, generate sets that conform to min/max"""
    values = paired_scores.head(head_size)
    scores = values.weighted_score
    pairs_and_sets = pairs_to_sets(values, max_size)
    print("Sizes: ", head_size, list(map(len, pairs_and_sets)))
    current_sets = list(filter(lambda gene_set: len(gene_set) >= min_size and
                               len(gene_set) <= max_size,
                               pairs_and_sets))

    feedback_frames = []
    for gene_set in current_sets:
        feedback_frames.append(pd.DataFrame({
            'genes_considered': [head_size],
            'genes_in_set': [len(gene_set)],
            'score_max': [np.max(scores)],
            'score_min': [np.min(scores)],
            'score_avg': [np.mean(scores)],
            'score_std': [np.std(scores)]
        }))

    return feedback_frames, current_sets


def unique_sets(existing_sets, dfs, new_sets):
    sets_confirmed = []
    dfs_confirmed = []
    for gene_set, df in zip(new_sets, dfs):
        if gene_set not in existing_sets:
            sets_confirmed.append(gene_set)
            dfs_confirmed.append(df)
    return dfs_confirmed, sets_confirmed

all_sets = []
all_dfs = []
head_sizes = list(np.arange(paired_scores.shape[0] // args.low)
                  * args.low + args.low)
print(head_sizes)
print("number of iterations:", len(head_sizes))
# For each given step, generate sets. If unique, add them to the list
for idx, head_size in enumerate(head_sizes):
    print("    - ", math.floor( 100 * idx / len(head_sizes)), "%")
    df, gene_sets = set_values(
        paired_scores, args.low, args.high, head_size)
    df, gene_sets = unique_sets(all_sets, df, gene_sets)
    if len(gene_sets) < 1:
        continue
    all_dfs = all_dfs + df
    all_sets = all_sets + gene_sets
    if len(all_sets) > args.count:
        break

if len(all_dfs) < 1:
    print("No valid gene sets generated")
    sys.exit(1)

feedback = rbind_all(all_dfs)
feedback = feedback.assign(
    cell_type=args.type, cell_name=cell_name, title=args.title)
feedback.index.name = 'index'
feedback.to_csv(args.feedback, sep='\t',
                encoding='utf-8')


for idx, gene_set in enumerate(all_sets):
    gene_set_lst = list(gene_set)
    frame = pd.DataFrame({'gene': gene_set_lst})
    filename_set = '{}.set.{}.{}.{:05d}.tsv'.format(
        args.title, args.type, cell_name, idx)
    path_set = os.path.join(args.output, filename_set)
    frame.to_csv(path_set, sep='\t',
                 encoding='utf-8', index=False)
