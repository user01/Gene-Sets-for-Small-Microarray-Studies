
import pandas as pd


def rbind(df_a, df_b):
    return pd.DataFrame.append(df_a, df_b) \
        .reset_index(drop=True)


def rbind_all(lst):
    return reduce(rbind, lst)
