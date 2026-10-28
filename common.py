import re, numpy as np, pandas as pd

RAW = '/home/claude/data/ds.xlsx'
COLS = {
    'Suhu': 'Suhu', 'Zat Padat Terlarut (TDS)': 'TDS', 'Zat Padat Tersuspensi (TSS)': 'TSS',
    'Ammonia (NH₃-N)': 'NH3N', 'BOD': 'BOD', 'COD': 'COD', 'Detergen (MBAS)': 'MBAS', 'DO': 'DO',
    'Fosfat (PO₄)': 'PO4', 'Minyak dan Lemak': 'MinyakLemak', 'Nitrat sebagai N (NO₃-N)': 'NO3N',
    'pH': 'pH', 'Sulfida': 'Sulfida', 'Timbal (Pb) Terlarut': 'Pb', 'Total Fecal Coliform': 'FecalColi',
}
P4 = ['Suhu', 'TDS', 'DO', 'pH']                                    # parameter lapangan
P8 = P4 + ['TSS', 'BOD', 'COD', 'NH3N']                             # + lab rutin
P15 = list(COLS.values())                                           # seluruh parameter
FSETS = {'P4': P4, 'P8': P8, 'P15': P15}


def parse(v):
    """return (nilai, flag) flag: 0 normal, -1 < LOD (disubstitusi LOD/2), +1 > batas ukur"""
    if isinstance(v, (int, float, np.number)) and not pd.isna(v):
        return float(v), 0
    s = str(v).strip().replace(' ', '').replace(',', '.')
    if s.startswith('<'):
        return float(s[1:]) / 2, -1
    if s.startswith('>'):
        return float(s[1:]), 1
    return float(s), 0


def load():
    df = pd.read_excel(RAW, 'Sheet1')
    out = pd.DataFrame({'lokasi': df['LOCATIONS'], 'tahun': df['YEAR']})
    flags = pd.DataFrame(index=df.index)
    for c, k in COLS.items():
        pv = df[c].map(parse)
        out[k] = pv.map(lambda t: t[0]); flags[k] = pv.map(lambda t: t[1])
    out['y'] = (df['TINGKAT PENCEMARAN'].str.strip() == 'CEMAR SEDANG').astype(int)
    # nama sungai tanpa nomor titik, untuk uji per lokasi
    out['sungai'] = df['LOCATIONS'].str.upper().str.replace('SUNAI', 'SUNGAI').str.replace('PALAWANGAN', 'PLAWANGAN') \
        .str.extract(r'SUNGAI\s+([A-Z]+)')[0]
    return out, flags


def scaled(df, cols, n, seed=0):
    """Simulasi skala: bootstrap baris data asli + derau Gaussian 5% simpangan baku per parameter."""
    rng = np.random.default_rng(seed)
    if n == len(df):
        return df[cols].to_numpy(float), df['y'].to_numpy()
    idx = rng.integers(0, len(df), n)
    X = df[cols].to_numpy(float)[idx]
    X = X + rng.normal(0, 1, X.shape) * df[cols].std().to_numpy() * 0.05
    return X, df['y'].to_numpy()[idx]


def make_model(name, seed=0):
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier
    from sklearn.ensemble import RandomForestClassifier, StackingClassifier
    from sklearn.neural_network import MLPClassifier
    from xgboost import XGBClassifier
    from catboost import CatBoostClassifier
    m = {
        'LR': lambda: make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000)),
        'SVM': lambda: make_pipeline(StandardScaler(), SVC(kernel='rbf', probability=False, random_state=seed)),
        'DT': lambda: DecisionTreeClassifier(random_state=seed),
        'RF': lambda: RandomForestClassifier(n_estimators=300, n_jobs=1, random_state=seed),
        'XGBoost': lambda: XGBClassifier(n_estimators=300, max_depth=4, learning_rate=0.1, n_jobs=1,
                                         random_state=seed, verbosity=0),
        'CatBoost': lambda: CatBoostClassifier(iterations=300, depth=6, thread_count=1, random_seed=seed,
                                               verbose=False, allow_writing_files=False),
        'MLP': lambda: make_pipeline(StandardScaler(), MLPClassifier((100, 50), max_iter=500, random_state=seed)),
    }
    if name == 'Stacking':
        base = [(k, make_model(k, seed)) for k in ['RF', 'XGBoost', 'CatBoost', 'SVM']]
        return StackingClassifier(base, final_estimator=LogisticRegression(max_iter=2000), cv=5, n_jobs=1)
    return m[name]()


MODELS = ['LR', 'SVM', 'DT', 'RF', 'XGBoost', 'CatBoost', 'MLP', 'Stacking']
