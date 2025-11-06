import pandas as pd
from qemd.io.ingest_any import load_expression_any

def load_expression(path):
    return load_expression_any(path)

def normalize_tpm(df):
    return df.apply(lambda x: (x - x.mean()) / (x.std(ddof=0) + 1e-8), axis=1)

def filter_pathway(df, genes):
    return df.loc[df.index.intersection(genes)]
