import numpy as np
import pandas as pd
import streamlit as st
from pathlib import Path
from qemd.mapping.batch import simulate_sample
from qemd.io.graph_io import load_graph
from qemd.io.ingest_any import load_expression_any

st.set_page_config(layout='wide', page_title='Quantum Bioenergetics Mapping')
st.title('Quantum Energy Mapping: Mitochondrial Transport')

graph_path = st.text_input('Graph JSON path', 'data/graph.json')
expr_path = st.text_input('Expression path (CSV/TSV/JSON; rows=genes, cols=samples)', 'data/sample_expr.tsv')

df_expr = None
if Path(expr_path).exists():
    try:
        df_expr = load_expression_any(expr_path)
        st.success(f'Parsed matrix: {df_expr.shape[0]} genes × {df_expr.shape[1]} samples')
        st.dataframe(df_expr.iloc[:10, :5])
    except Exception as e:
        st.error(f'Failed to parse: {e}')

sample = ''
if df_expr is not None:
    samples = list(df_expr.columns)
    sample = st.selectbox('Sample', samples)

col1, col2, col3 = st.columns(3)
with col1:
    sigma = st.slider('Static disorder σ', 0.0, 0.05, 0.01, 0.002)
    k_sink = st.slider('k_sink', 0.0, 0.5, 0.1, 0.01)
with col2:
    k_loss = st.slider('k_loss', 0.0, 0.1, 0.01, 0.005)
    eps0 = st.slider('ε0', -0.1, 0.1, 0.0, 0.005)
with col3:
    alpha = st.slider('α (energy slope)', 0.0, 0.1, 0.02, 0.002)
    J0 = st.slider('J0 (coupling scale)', 0.0, 0.1, 0.02, 0.002)

Jmax = 0.05
source_site = st.number_input('Source site index', 0, 20, 0)
sink_site = st.number_input('Sink site index', 0, 20, 5)
sweep_gammas = np.round(np.arange(0.00, 0.05+1e-12, 0.0025), 6)

if st.button('Run simulation') and Path(graph_path).exists() and df_expr is not None and sample:
    expr_vec = df_expr[sample]
    G = load_graph(graph_path)
    etc_genes = sorted({g for genes in G['node_genes'] for g in genes})
    overlap = set(expr_vec.index).intersection(etc_genes)
    if len(overlap) == 0:
        st.error('No overlap between your matrix genes and ETC mapping table genes. Please use an expression matrix, not clinical JSON.')
    else:
        df, ETE_peak, gamma_star = simulate_sample(expr_vec, graph_path, sweep_gammas, sigma, k_sink, k_loss, source_site, sink_site, eps0, alpha, J0, Jmax)
        st.line_chart(df.set_index('gamma')['ETE'])
        st.write(f'ETE_peak: {ETE_peak:.3f} at γ*={gamma_star:.3f}')
        st.line_chart(df.set_index('gamma')['tau_c'])
        st.json(df.to_dict(orient='records'))
