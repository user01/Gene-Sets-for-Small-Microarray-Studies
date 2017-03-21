
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

# args = parser.parse_args()
args = parser.parse_args(([
    '--input', 'results',
    '--raw', 'results/gene_data_vs_cell_type.tsv'
]))


raw_data = pd.read_table(args.raw)

cell_types = list(np.unique(raw_data.Cell_Type))
general_cell_types = list(map(lambda s: s.replace(' ', '_'),
                              list(np.unique(raw_data.General_Cell_Type))))


def interpret_results(df, cell_group, cell_name):
    col_resets = {}
    col_resets[cell_group] = 'truth'
    col_resets['{}_Predicted'.format(cell_group)] = 'predicted'
    df_fixed = df.rename(columns=col_resets)[['truth', 'predicted']]
    df_rated = df_fixed.assign(
        correct=df_fixed.truth == df_fixed.predicted)
    accuracy = np.sum(df_rated.correct) / df_fixed.shape[0]
    df_sensitivity = df_rated.query('truth == "{}"'.format(cell_name))
    sensitivity = np.sum(df_sensitivity.correct) / \
        df_sensitivity.shape[0]
    return pd.DataFrame({'group': [cell_group], 'name': [cell_name],
                         'accuracy': [accuracy], 'sensitivity': [sensitivity],
                         'set_score': [accuracy * sensitivity]})


info_regex = re.compile('.+set.results.\w+.\w+.(\d+).(.\w+).tsv')


def read_results(path, cell_group, cell_name):
    set_id, classifier = info_regex.match(path).groups()
    results = pd.read_table(path)
    return interpret_results(results, cell_group, cell_name) \
        .assign(set_id=int(set_id), classifier=classifier)


def read_glob(current_glob, cell_group, cell_name):
    result_paths = glob.glob(current_glob)
    results_dfs = map(lambda path: read_results(
        path, cell_group, cell_name), result_paths)
    return rbind_all(list(results_dfs))


def read_feedback(cell_group, cell_name):
    feedback_file = 'set.feedback.{}.{}.tsv'.format(
        cell_group, cell_name)
    feedback_path = os.path.join(args.input, feedback_file)
    feedback = pd.read_table(feedback_path).rename(
        columns={'index': 'set_id'})

    data_paths = map(lambda idx: 'set.results.{}.{}.{:05d}.*.tsv'.format(
        cell_group, cell_name, idx), feedback.index)
    data_paths_full = map(lambda filename: os.path.join(
        args.input, filename), data_paths)
    results = rbind_all(list(map(lambda g: read_glob(
        g, cell_group, cell_name), data_paths_full)))
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

for cell_name in general_cell_types[4:5]:
    pp = read_feedback('General_Cell_Type', cell_name)
    break



pp.iloc[0:1]
