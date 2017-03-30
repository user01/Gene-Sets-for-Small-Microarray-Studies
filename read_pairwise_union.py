
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
    description='Read leader results and combine based on titles.')

parser.add_argument('--input', type=str, required=True,
                    help='Glob to input leader results')

parser.add_argument('--output', type=str, required=True,
                    help='Path to output combined tsv')


args = parser.parse_args()
# args = parser.parse_args(([
#     '--input',
#     '**/**/*.sets.leaders.tsv',
#     '--output',
#     'results/combined.tsv'
# ]))

paths = glob.glob(args.input.replace('"',''))

def leader_conversion(path):
    df = pd.read_table(path) \
        .drop(['Unnamed: 0', 'set_id'], 1)
    title = df.title[0]
    return df \
        .assign(name = df.cell_type + '|' + df.cell_name) \
        .rename(columns={
            'genes_considered': 'genes_considered_' + title,
            'genes_in_set': 'genes_in_set_' + title,
            'score_train_Macrophage': 'train_score_Macrophage_' + title,
            'score_validation_Macrophage': 'validation_score_Macrophage_' + title,
            'score_train_Microglia': 'train_score_Microglia_' + title,
            'score_validation_Microglia': 'validation_score_Microglia_' + title,
            'score_train_Neutrophil': 'train_score_Neutrophil_' + title,
            'score_validation_Neutrophil': 'validation_score_Neutrophil_' + title,
            'score_train_Monocyte': 'train_score_Monocyte_' + title,
            'score_validation_Monocyte': 'validation_score_Monocyte_' + title,
            'base_score': 'train_score_' + title,
            'validation_score': 'validation_score_' + title
        }) \
        .drop(['title', 'cell_type', 'cell_name'], 1) \
        .set_index('name')


def leader_all(paths):
    dfs = list(map(leader_conversion, paths))
    df_all = pd.concat(dfs, axis=1)
    return df_all[sorted(df_all.columns)]

leader_all(paths).to_csv(args.output, sep='\t', encoding='utf-8')
