import pandas as pd
import json
import numpy as np
from pathlib import Path

def load_expression_any(path_or_buf):
    if isinstance(path_or_buf, (str, Path)):
        p = str(path_or_buf)
        if p.lower().endswith(('.tsv', '.txt')):
            df = pd.read_csv(p, sep='\t', index_col=0)
        elif p.lower().endswith('.csv'):
            try:
                df = pd.read_csv(p, index_col=0)
            except Exception:
                df = pd.read_csv(p, sep='\t', index_col=0)
        elif p.lower().endswith('.json'):
            with open(p, 'r') as f:
                data = json.load(f)
            df = normalize_json_to_matrix(data)
        else:
            df = pd.read_csv(p, sep=None, engine='python', index_col=0)
    else:
        try:
            df = pd.read_csv(path_or_buf, index_col=0)
        except Exception:
            df = pd.read_csv(path_or_buf, sep='\t', index_col=0)
    df = df.apply(pd.to_numeric, errors='coerce').fillna(0.0)
    return df

def normalize_json_to_matrix(data):
    if isinstance(data, list):
        df = pd.DataFrame(data)
        for gcol in ['gene','hugoGeneSymbol','symbol']:
            if gcol in df.columns: gene_col = gcol; break
        else: gene_col = None
        val_col = None
        for v in ['expression','value','TPM','FPKM','count','counts']:
            if v in df.columns: val_col = v; break
        samp_col = None
        for s in ['sample','sampleId','barcode']:
            if s in df.columns: samp_col = s; break
        if gene_col and val_col and samp_col:
            mat = df.pivot(index=gene_col, columns=samp_col, values=val_col)
            return mat
        return df.set_index(df.columns[0])
    elif isinstance(data, dict):
        first = next(iter(data.values()))
        if isinstance(first, dict):
            inner_keys = list(first.keys())
            if any('-' in k or k.isupper() for k in inner_keys):
                return pd.DataFrame(data)
            else:
                return pd.DataFrame.from_dict(data, orient='index').T
        else:
            raise ValueError("Unsupported JSON dict shape")
    else:
        raise ValueError("Unsupported JSON root type")
