
import argparse
import math
import glob
import os
import sys
import re

import numpy as np
import pandas as pd

from utilities.common import rbind_all, rbind


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
#     '--input', 'results',
#     '--title', 'main',
#     '--raw', 'results/gene_data_vs_cell_type.tsv',
#     '--outputfull', 'results/set.full.tsv',
#     '--outputleader', 'results/set.leader.tsv',
#     '--outputsets', 'results/sets.gmt'
# ]))


raw_data = pd.read_table(args.raw)


def fix_name_spaces(arr):
    return list(map(lambda s: s.replace(' ', '_'),
                    list(np.unique(arr))))

cell_names = fix_name_spaces(raw_data.Cell_Type)
general_cell_names = fix_name_spaces(raw_data.General_Cell_Type)


def interpret_results(df, cell_group, cell_name):
    col_resets = {}
    col_resets[cell_group] = 'truth'
    col_resets['{}_Predicted'.format(cell_group)] = 'predicted'
    df_fixed = df.rename(columns=col_resets)[['truth', 'predicted']]
    df_rated = df_fixed.assign(
        correct=df_fixed.truth == df_fixed.predicted)
    accuracy = np.sum(df_rated.correct) / df_fixed.shape[0]
    df_sensitivity = df_rated.query(
        'truth == "{}"'.format(cell_name.replace('_', ' ')))
    sensitivity = np.sum(df_sensitivity.correct) / \
        df_sensitivity.shape[0]

    # True Negative / (True Negative + False Positive)
    df_specificity = df_rated.query(
        'truth != "{}"'.format(cell_name.replace('_', ' ')))
    df_specificity = df_specificity.assign(truenegative = df_specificity.predicted != cell_name)
    specificity = np.sum(df_specificity.truenegative) / \
        df_specificity.shape[0]

    return pd.DataFrame({'group': [cell_group], 'name': [cell_name],
                         'accuracy': [accuracy], 'sensitivity': [sensitivity],
                         'specificity': [specificity],
                         'set_score': [specificity + 2 * sensitivity]})


info_regex = re.compile('.+set.results.\w+.\w+.(\d+).(.\w+).tsv')


def read_results(path, cell_group, cell_name):
    set_id, classifier = info_regex.match(path).groups()
    results = pd.read_table(path)
    return interpret_results(results, cell_group, cell_name) \
        .assign(set_id=int(set_id), classifier=classifier)


def read_glob(current_glob, cell_group, cell_name):
    result_paths = glob.glob(current_glob)
    if len(result_paths) < 1:
        return None
    results_dfs = map(lambda path: read_results(
        path, cell_group, cell_name), result_paths)
    return rbind_all(list(results_dfs))


def read_feedback(cell_group, cell_name):
    feedback_file = '{}.set.feedback.{}.{}.tsv'.format(
        args.title, cell_group, cell_name)
    feedback_path = os.path.join(args.input, feedback_file)
    if not os.path.isfile(feedback_path):
        return None
    feedback = pd.read_table(feedback_path).rename(
        columns={'index': 'set_id'})

    data_paths = map(lambda idx: '{}.set.results.{}.{}.{:05d}.*.tsv'.format(
        args.title, cell_group, cell_name, idx), feedback.index)
    data_paths_full = map(lambda filename: os.path.join(
        args.input, filename), data_paths)
    results = rbind_all(list(map(lambda g: read_glob(
        g, cell_group, cell_name), data_paths_full)))
    if results is None:
        return None

    results_all = pd.DataFrame.merge(feedback, results, on="set_id")
    return results_all.groupby(
        ['group', 'name', 'set_id', 'genes_in_set', 'genes_considered'], as_index=False) \
        .agg({
            'set_score': np.mean,
            'score_std': np.max,
            'score_avg': np.max
        }) \
        .sort_values('set_score', ascending=False) \
        .reset_index(drop=True)


def read_types(cell_type, cell_names):
    dfs_full = []
    dfs_leader = []
    for cell_name in cell_names:
        df = read_feedback(cell_type, cell_name)
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
    sys.exit(0)

df_full.to_csv(args.outputfull, sep='\t', encoding='utf-8')
df_leader.to_csv(args.outputleader, sep='\t', encoding='utf-8')

# Write out picked sets to a final file
lines = []
for row in df_leader.iterrows():
    content = row[1]
    set_filename = '{}.set.{}.{}.{:05d}.tsv'.format(
        args.title, content['group'], content['name'], content['set_id'])
    set_path = os.path.join(args.input, set_filename)
    set_data = '\t'.join(list(pd.read_table(set_path).gene))
    leader = '\t'.join(
        [content['group'], content['name'], str(content['genes_in_set'])])
    line = '{}\t\t{}'.format(leader, set_data)
    lines.append(line)

with open(args.outputsets, 'w') as file_handler:
    for line in lines:
        file_handler.write("{}\n".format(line))
