# Recovering undocumented specifications of a regulatory water quality index

Code for the article *Recovering Undocumented Specifications of a Regulatory Water Quality Index:
A Censoring-Aware and Tested Computational Pipeline Applied to Rivers in Central Java, Indonesia*
(submitted to the Journal of Universal Computer Science).

The study covers 54 samples (2023–2025) from 16 rivers monitored by BPSDA Pemali Comal, Central Java, Indonesia: Boyo, Cacaban, Comal, Gangsa, Gung, Pakijangan, Pekalongan, Pemali,
Plawangan, Rambut, Sambong, Sengkarang, Sragi, Urang, Wadas and Waluh.

## Contents

| File | Purpose |
|---|---|
| `wqindex.py` | Censoring-aware measurement type (`Censored`), index specification (`Spec`), Pollution Index (Decree 115/2003) and IKA (Regulation 27/2021) — Algorithm 1 |
| `specsearch.py` | Exhaustive specification recovery over 786,384 configurations, identifiability, hold-out validation, ablation — Algorithm 2 |
| `analysis.py` | Status, determinants (leave-one-parameter-out), censoring sensitivity, agreement (Cohen's kappa), figures |
| `common.py` | Data loading and parsing of censored values (`<LOQ`, `>upper bound`) |
| `tests/test_wqindex.py` | Oracle, regression and metamorphic tests |
| `CHANGELOG.md` | Changes between releases |
| `data/README.md` | How to obtain the data and the expected file format |

## Installation

```bash
python -m pip install -r requirements.txt
```

Python 3.10 or later.

## Verified labels

The reported class of one sample (2023, downstream Rambut) could not be verified against the laboratory
summary tables. It is listed in `common.UNVERIFIED` and excluded from every label-based analysis, which
therefore uses 53 verified labels.

## Data

The laboratory data are **not distributed** here; they can be requested from BPSDA Pemali Comal with its
written permission (see `data/README.md`). Place the file at `data/dataset.xlsx` or set `WQ_DATA`.

## Usage

```bash
export WQ_DATA=data/dataset.xlsx     # optional if the file is at this default location
python analysis.py                   # writes outputs/res.json, outputs/status_per_sample.csv, figures
python specsearch.py                 # writes outputs/spec.json (about 4 s)
pytest -q                            # 81 tests with the data
```

Without the data, `pytest -q` runs the metamorphic and unit tests on synthetic samples
(9 passed, 3 skipped), so the index implementation can be checked by anyone.

## Specification used in the index

A specification is c = (S, ρ, λ, μ, γ): parameter set, DO rule (`sat` or `inv`), left-censoring factor,
right-censoring multiplier and ratio correction.

```python
from wqindex import Censored, Spec, LAB_SPEC, pi, pi_class
# a sample maps every parameter in the specification (and 'Suhu') to a Censored value
sample = {'Suhu': Censored(28.5), 'pH': Censored(7.4), 'DO': Censored(5.1), 'TDS': Censored(180),
          'TSS': Censored(24), 'NH3N': Censored.parse('<0,09'), 'BOD': Censored(6.2), 'COD': Censored(21),
          'MBAS': Censored.parse('<0,07'), 'PO4': Censored(0.12), 'MinyakLemak': Censored(1.4),
          'NO3N': Censored.parse('<1'), 'Sulfida': Censored(0.06), 'Pb': Censored(0.05),
          'FecalColi': Censored.parse('>16000')}
pi(sample, Spec())        # documented specification (14 parameters, LOQ/2)
pi(sample, LAB_SPEC)      # specification recovered from the laboratory reports
```

## Citation

See `CITATION.cff`. Please cite the article once published.

## License

Code: GPL-3.0 (see `LICENSE`). The data are not covered by this license.
