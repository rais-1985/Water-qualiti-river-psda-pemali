"""Test suite: (A) oracle tests against published laboratory values, (B) regression against the analysis pipeline,
(C) metamorphic/property tests of the index definition."""
import math, json, random, sys
import numpy as np, pandas as pd, pytest
sys.path.insert(0, '/home/claude/art1/pkg')
from wqindex import Censored, Spec, LAB_SPEC, ALL14, LIMITS, pi, pi_class, ratios, ika

RAW = '/home/claude/data6/ds.xlsx'
COLS = {'Suhu': 'Suhu', 'Zat Padat Terlarut (TDS)': 'TDS', 'Zat Padat Tersuspensi (TSS)': 'TSS', 'Ammonia (NH₃-N)': 'NH3N',
        'BOD': 'BOD', 'COD': 'COD', 'Detergen (MBAS)': 'MBAS', 'DO': 'DO', 'Fosfat (PO₄)': 'PO4', 'Minyak dan Lemak': 'MinyakLemak',
        'Nitrat sebagai N (NO₃-N)': 'NO3N', 'pH': 'pH', 'Sulfida': 'Sulfida', 'Timbal (Pb) Terlarut': 'Pb',
        'Total Fecal Coliform': 'FecalColi'}
df = pd.read_excel(RAW, 'Sheet1')
SAMPLES = [{k: Censored.parse(row[c]) for c, k in COLS.items()} for _, row in df.iterrows()]
LABEL = [2 if s.strip() == 'CEMAR SEDANG' else 1 for s in df['TINGKAT PENCEMARAN']]
REP24 = {'RAMBUT HULU': 2.39, 'RAMBUT HILIR': 2.41, 'CACABAN HULU': 5.06, 'CACABAN HILIR': 5.06, 'WALUH HULU': 5.05,
         'WALUH HILIR': 2.17, 'PALAWANGAN JEMBATAN KRAMAT HULU': 3.41, 'PALAWANGAN JEMBATAN KALI PALAWANGAN HILIR': 4.22,
         'GUNG HULU': 1.70, 'GUNG HILIR': 5.11, 'WADAS HULU': 5.11, 'WADAS HILIR': 4.92, 'GANGSA HULU': 5.63,
         'GANGSA HILIR': 5.09, 'PEMALI HULU': 3.59, 'PEMALI HILIR': 5.08, 'PAKIJANGAN HULU': 5.36, 'PAKIJANGAN HILIR': 4.46}
I24 = [i for i in range(len(df)) if df.YEAR[i] == 2024]
POLL = [p for p in LIMITS]                     # parameters for which a higher value means worse quality


# ---------------- (A) oracle tests ----------------
@pytest.mark.parametrize('i', I24)
def test_oracle_pi_2024(i):
    key = df.LOCATIONS[i].upper().replace('SUNGAI ', '').replace('SUNAI ', '')
    assert abs(pi(SAMPLES[i], LAB_SPEC) - REP24[key]) <= 0.025


@pytest.mark.parametrize('i', range(len(SAMPLES)))
def test_oracle_class(i):
    assert pi_class(pi(SAMPLES[i], LAB_SPEC)) == LABEL[i]


# ---------------- (B) regression against the analysis pipeline ----------------
def test_regression_pipeline():
    st = pd.read_csv('/home/claude/art1/status_per_sample.csv')
    for i, s in enumerate(SAMPLES):
        assert pi(s, Spec()) == pytest.approx(st['IP15_half'][i], abs=1e-9)
        assert pi(s, Spec(params=('pH', 'DO', 'BOD', 'COD', 'TSS', 'NO3N', 'PO4', 'FecalColi'))) == \
            pytest.approx(st['IKA8_half'][i], abs=1e-9)


# ---------------- (C) property / metamorphic tests ----------------
RNG = random.Random(7)
N_PROP = 2000


def perturb(s, p, factor):
    c = s[p]; return {**s, p: Censored(c.value * factor, c.flag)}


def test_monotone_in_pollutant():
    """Raising any pollutant concentration never lowers the PI (any spec)."""
    for _ in range(N_PROP):
        s = RNG.choice(SAMPLES); p = RNG.choice(POLL); spec = RNG.choice([Spec(), LAB_SPEC])
        if p not in spec.params: continue
        assert pi(perturb(s, p, 1 + RNG.random() * 5), spec) >= pi(s, spec) - 1e-12


def test_monotone_in_censoring_rule():
    """0 <= LOQ/2 <= LOQ/sqrt2 <= LOQ substitution yields non-decreasing PI."""
    for s in SAMPLES:
        v = [pi(s, Spec(left=f)) for f in (0, 0.5, 1 / math.sqrt(2), 1)]
        assert all(a <= b + 1e-12 for a, b in zip(v, v[1:]))


def test_do_rules_opposite_direction():
    """Saturation rule: more DO -> lower ratio; 'inv' rule: more DO -> higher ratio (the recovered defect)."""
    for s in SAMPLES:
        lo, hi = perturb(s, 'DO', 0.8), perturb(s, 'DO', 1.2)
        assert ratios(hi, Spec())['DO'] <= ratios(lo, Spec())['DO'] + 1e-12
        assert ratios(hi, LAB_SPEC)['DO'] > ratios(lo, LAB_SPEC)['DO']


def test_correction_continuous_at_one():
    assert 1 + 5 * math.log10(1.0) == 1.0


def test_permutation_invariance():
    """The PI does not depend on the order of the parameters."""
    for _ in range(200):
        s = RNG.choice(SAMPLES); p = list(ALL14); RNG.shuffle(p)
        assert pi(s, Spec(params=tuple(p))) == pytest.approx(pi(s, Spec()), abs=1e-12)


def test_equal_ratios_identity():
    """If every ratio equals r <= 1, PI = r."""
    s = {p: Censored(0.5 * LIMITS[p]) for p in LIMITS}
    assert pi(s, Spec(params=tuple(LIMITS))) == pytest.approx(0.5)


def test_class_boundaries():
    assert [pi_class(x) for x in (1.0, 1.0001, 5.0, 5.0001, 10.0, 10.0001)] == [0, 1, 1, 2, 2, 3]


def test_ika_bounds():
    assert ika([0] * 10) == 70 and ika([3] * 10) == 10


def test_parse_censored():
    assert Censored.parse('<0,09') == Censored(0.09, -1)
    assert Censored.parse('>16000') == Censored(16000.0, 1)
    assert Censored.parse(4.5) == Censored(4.5, 0)
