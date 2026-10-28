"""Data loading for the compiled laboratory dataset.

The dataset is not distributed with this repository (see data/README.md).
Set its location with the environment variable WQ_DATA, or place it at data/dataset.xlsx.
"""
import os
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.environ.get('WQ_DATA', os.path.join(HERE, 'data', 'dataset.xlsx'))
OUT = os.environ.get('WQ_OUT', os.path.join(HERE, 'outputs'))

COLS = {
    'Suhu': 'Suhu', 'Zat Padat Terlarut (TDS)': 'TDS', 'Zat Padat Tersuspensi (TSS)': 'TSS',
    'Ammonia (NH₃-N)': 'NH3N', 'BOD': 'BOD', 'COD': 'COD', 'Detergen (MBAS)': 'MBAS', 'DO': 'DO',
    'Fosfat (PO₄)': 'PO4', 'Minyak dan Lemak': 'MinyakLemak', 'Nitrat sebagai N (NO₃-N)': 'NO3N',
    'pH': 'pH', 'Sulfida': 'Sulfida', 'Timbal (Pb) Terlarut': 'Pb', 'Total Fecal Coliform': 'FecalColi',
}


def data_available():
    return os.path.exists(RAW)


def parse(v):
    """Return (value, flag): flag 0 quantified, -1 '<LOQ' (value stored as LOQ/2), +1 '>upper bound'."""
    if isinstance(v, (int, float, np.number)) and not pd.isna(v):
        return float(v), 0
    s = str(v).strip().replace(' ', '').replace(',', '.')
    if s.startswith('<'):
        return float(s[1:]) / 2, -1
    if s.startswith('>'):
        return float(s[1:]), 1
    return float(s), 0


def load(path=None):
    """Load the dataset into (values, flags) DataFrames with normalised parameter and river names."""
    df = pd.read_excel(path or RAW, 'Sheet1')
    out = pd.DataFrame({'lokasi': df['LOCATIONS'], 'tahun': df['YEAR']})
    flags = pd.DataFrame(index=df.index)
    for c, k in COLS.items():
        pv = df[c].map(parse)
        out[k] = pv.map(lambda t: t[0]); flags[k] = pv.map(lambda t: t[1])
    out['y'] = (df['TINGKAT PENCEMARAN'].str.strip() == 'CEMAR SEDANG').astype(int)
    out['sungai'] = df['LOCATIONS'].str.upper().str.replace('SUNAI', 'SUNGAI').str.replace('PALAWANGAN', 'PLAWANGAN') \
        .str.extract(r'SUNGAI\s+([A-Z]+)')[0]
    return out, flags
