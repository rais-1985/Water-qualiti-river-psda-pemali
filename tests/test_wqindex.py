"""Test suite for wqindex.

(A) Oracle tests against values published by the laboratory  -> need the dataset (skipped otherwise);
    class oracles use the 53 verified labels only (common.UNVERIFIED).
(B) Regression test against the analysis pipeline               -> needs the dataset and outputs/status_per_sample.csv.
(C) Metamorphic and unit tests                                   -> run on synthetic samples, and also on the
                                                                    dataset when it is available.
Run:  pytest -q           (set WQ_DATA=/path/to/dataset.xlsx to include A and B)
"""
import math, os, random, sys
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import common
from wqindex import Censored, Spec, LAB_SPEC, ALL14, LIMITS, pi, pi_class, ratios, ika

HAVE_DATA = common.data_available()
needs_data = pytest.mark.skipif(not HAVE_DATA, reason='dataset not available (see data/README.md)')

# PI values published in the calculation tables of the 2024 monitoring report (BPDSDA Pemali Comal / CDULab)
REP24 = {'RAMBUT HULU': 2.39, 'RAMBUT HILIR': 2.41, 'CACABAN HULU': 5.06, 'CACABAN HILIR': 5.06, 'WALUH HULU': 5.05,
         'WALUH HILIR': 2.17, 'PALAWANGAN JEMBATAN KRAMAT HULU': 3.41, 'PALAWANGAN JEMBATAN KALI PALAWANGAN HILIR': 4.22,
         'GUNG HULU': 1.70, 'GUNG HILIR': 5.11, 'WADAS HULU': 5.11, 'WADAS HILIR': 4.92, 'GANGSA HULU': 5.63,
         'GANGSA HILIR': 5.09, 'PEMALI HULU': 3.59, 'PEMALI HILIR': 5.08, 'PAKIJANGAN HULU': 5.36, 'PAKIJANGAN HILIR': 4.46}
TOL = 0.025


# ---------------- samples ----------------
def synthetic(n=60, seed=11):
    """Random samples with plausible ranges, including left- and right-censored values."""
    rng = random.Random(seed); out = []
    for _ in range(n):
        s = {'Suhu': Censored(rng.uniform(24, 32)), 'pH': Censored(rng.uniform(6.2, 8.8)), 'DO': Censored(rng.uniform(2, 8))}
        for p, lim in LIMITS.items():
            u = rng.random()
            if u < 0.15: s[p] = Censored(lim * rng.uniform(0.3, 2.0), -1)
            elif p == 'FecalColi' and u > 0.9: s[p] = Censored(16000.0, 1)
            else: s[p] = Censored(lim * rng.lognormvariate(0, 1))
        out.append(s)
    return out


def real():
    df = pd.read_excel(common.RAW, 'Sheet1')
    samples = [{k: Censored.parse(row[c]) for c, k in common.COLS.items()} for _, row in df.iterrows()]
    labels = [2 if s.strip() == 'CEMAR SEDANG' else 1 for s in df['TINGKAT PENCEMARAN']]
    return df, samples, labels


DF, REAL, LABELS = real() if HAVE_DATA else (None, [], [])
I24 = [i for i in range(len(REAL)) if DF.YEAR[i] == 2024]
# labels verified against the reports (the class of the 2023 downstream Rambut point was not verifiable)
IVER = [i for i in range(len(REAL)) if (int(DF.YEAR[i]), str(DF.LOCATIONS[i]).strip().upper()) not in common.UNVERIFIED]
SAMPLES = synthetic() + REAL
POLL = list(LIMITS)
RNG = random.Random(7)


# ---------------- (A) oracle tests ----------------


@needs_data
@pytest.mark.parametrize('i', I24)
def test_oracle_pi_2024(i):
    key = DF.LOCATIONS[i].upper().replace('SUNGAI ', '').replace('SUNAI ', '')
    assert abs(pi(REAL[i], LAB_SPEC) - REP24[key]) <= TOL, key


@needs_data
@pytest.mark.parametrize('i', IVER)
def test_oracle_class(i):
    assert pi_class(pi(REAL[i], LAB_SPEC)) == LABELS[i]


# ---------------- (B) regression ----------------
@needs_data
def test_regression_pipeline():
    path = os.path.join(common.OUT, 'status_per_sample.csv')
    if not os.path.exists(path):
        pytest.skip('run analysis.py first')
    st = pd.read_csv(path)
    for i, s in enumerate(REAL):
        assert pi(s, Spec()) == pytest.approx(st['IP15_half'][i], abs=1e-9)
        assert pi(s, Spec(params=('pH', 'DO', 'BOD', 'COD', 'TSS', 'NO3N', 'PO4', 'FecalColi'))) == \
            pytest.approx(st['IKA8_half'][i], abs=1e-9)


# ---------------- (C) metamorphic and unit tests ----------------
def perturb(s, p, factor):
    c = s[p]; return {**s, p: Censored(c.value * factor, c.flag)}


def test_mr1_monotone_in_pollutant():
    for _ in range(2000):
        s = RNG.choice(SAMPLES); p = RNG.choice(POLL); spec = RNG.choice([Spec(), LAB_SPEC])
        if p not in spec.params: continue
        assert pi(perturb(s, p, 1 + RNG.random() * 5), spec) >= pi(s, spec) - 1e-12


def test_mr2_monotone_in_censoring_rule():
    for s in SAMPLES:
        v = [pi(s, Spec(left=f)) for f in (0, 0.5, 1 / math.sqrt(2), 1)]
        assert all(a <= b + 1e-12 for a, b in zip(v, v[1:]))


def test_mr3_do_rules_opposite_direction():
    for s in SAMPLES:
        lo, hi = perturb(s, 'DO', 0.8), perturb(s, 'DO', 1.2)
        assert ratios(hi, Spec())['DO'] <= ratios(lo, Spec())['DO'] + 1e-12
        assert ratios(hi, LAB_SPEC)['DO'] > ratios(lo, LAB_SPEC)['DO']


def test_mr4_permutation_invariance():
    for _ in range(200):
        s = RNG.choice(SAMPLES); p = list(ALL14); RNG.shuffle(p)
        assert pi(s, Spec(params=tuple(p))) == pytest.approx(pi(s, Spec()), abs=1e-12)


def test_mr5_equal_ratios_identity():
    s = {p: Censored(0.5 * LIMITS[p]) for p in LIMITS}
    assert pi(s, Spec(params=tuple(LIMITS))) == pytest.approx(0.5)


def test_correction_continuous_at_one():
    s = {p: Censored(1.0 * LIMITS[p]) for p in LIMITS}
    assert pi(s, Spec(params=tuple(LIMITS))) == pytest.approx(1.0)


def test_class_boundaries():
    assert [pi_class(x) for x in (1.0, 1.0001, 5.0, 5.0001, 10.0, 10.0001)] == [0, 1, 1, 2, 2, 3]


def test_ika_bounds():
    assert ika([0] * 10) == 70 and ika([3] * 10) == 10


def test_parse_censored():
    assert Censored.parse('<0,09') == Censored(0.09, -1)
    assert Censored.parse('>16000') == Censored(16000.0, 1)
    assert Censored.parse(4.5) == Censored(4.5, 0)
