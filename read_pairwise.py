
import argparse
import math
import glob
import os
import sys
import re

import numpy as np
import pandas as pd

from utilities.common import rbind_all


parser = argparse.ArgumentParser(
    description='Read set results and choose winners.')

parser.add_argument('--raw', type=str, required=True,
                    help='Path to raw gene data')
parser.add_argument('--input', type=str, required=True,
                    help='Path to input set files')
parser.add_argument('--title', type=str, required=True,
                    help='Title string')
parser.add_argument('--outputfull', type=str, required=True,
                    help='Path to output full results of each set')
parser.add_argument('--outputleader', type=str, required=True,
                    help='Path to output leader results of each set')
parser.add_argument('--outputsets', type=str, required=True,
                    help='Path to output gmt of sets')


args = parser.parse_args()
# args = parser.parse_args(([
#     '--raw',
#     'results/full/gene_data_vs_cell_type.tsv',
#     '--title',
#     'full',
#     '--input',
#     'results/full',
#     '--outputfull',
#     'results/full/full.sets.full.tsv',
#     '--outputleader',
#     'results/full/full.sets.leaders.tsv',
#     '--outputsets',
#     'results/full/full.sets.gmt'
# ]))


raw_data = pd.read_table(args.raw)


def fix_name_spaces(arr):
    """Replace spaces with underscores for filenames"""
    return list(map(lambda s: s.replace(' ', '_'),
                    list(np.unique(arr))))

cell_names = fix_name_spaces(raw_data.Cell_Type)
general_cell_names = fix_name_spaces(raw_data.General_Cell_Type)


def read_results(cell_group, cell_name):
    """Read a results file - a feedback file with the associated scores"""
    results_file = '{}.set.results.{}.{}.tsv'.format(
        args.title, cell_group, cell_name)
    results_path = os.path.join(args.input, results_file)
    if not os.path.isfile(results_path):
        return None
    return pd.read_table(results_path).rename(
        columns={'index': 'set_id'}) \
        .sort_values('base_scores', ascending=False) \
        .reset_index(drop=True)



def read_types(cell_type, cell_names):
    dfs_full = []
    dfs_leader = []
    for cell_name in cell_names:
        df = read_results(cell_type, cell_name)
        if df is None:
            continue
        dfs_full.append(df)
        dfs_leader.append(df.iloc[0:1])
    return dfs_full, dfs_leader


df_general_full, df_general_leader = read_types(
    'General_Cell_Type', general_cell_names)
df_cell_full, df_cell_leader = read_types('Cell_Type', cell_names)

df_full = rbind_all(df_general_full + df_cell_full)
df_leader = rbind_all(df_general_leader + df_cell_leader)

if df_full is None:
    print("No data supplied to test")
    sys.exit(1)

df_full.to_csv(args.outputfull, sep='\t', encoding='utf-8')
df_leader.to_csv(args.outputleader, sep='\t', encoding='utf-8')

# Write out picked sets to a final file
lines = []
for row in df_leader.iterrows():
    content = row[1]
    set_filename = '{}.set.{}.{}.{:05d}.tsv'.format(
        args.title, content['cell_type'], content['cell_name'], content['set_id'])
    set_path = os.path.join(args.input, set_filename)
    set_data = '\t'.join(list(pd.read_table(set_path).gene))
    leader = '\t'.join(
        [content['cell_type'], content['cell_name'], str(content['genes_in_set'])])
    line = '{}\t\t{}'.format(leader, set_data)
    lines.append(line)

with open(args.outputsets, 'w') as file_handler:
    for line in lines:
        file_handler.write("{}\n".format(line))
