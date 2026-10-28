# Changelog

## v1.1.0
- Label-based analyses use the 53 verified reported classes. The class of the 2023 downstream Rambut
  point could not be verified against the laboratory summary tables and is listed in
  `common.UNVERIFIED`; it is excluded from objective O1, the class oracle tests, the hold-out
  validation and the agreement statistics (hold-out: 35 independent labels).
- `specsearch.py`: counts equivalence classes of configurations (distinct PI and class vectors).
- `analysis.py`: river-level cluster bootstrap CIs; numerical check of Lemma 1 (bound on the PI change
  when a non-maximal parameter is removed); exceedance counts split into quantified and
  right-censored values; sensitivity of the one sulfide value reported as 0 mg/L.
- Test suite: 81 tests with the data (18 PI oracles, 53 class oracles, 1 regression, 5 metamorphic,
  4 unit); 9 tests run on synthetic samples without the data.

## v1.0.1
- First archived release.
