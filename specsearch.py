"""Specification recovery: exhaustive search over the configuration space of the Pollution Index.

A configuration c = (S, rho_DO, kappa_L, kappa_R, gamma):
  S       non-empty subset of the 14 candidate parameters
  rho_DO  DO ratio rule: 'sat' (saturation-based, Kepmen 115/2003) or 'inv' (C/L = DO/4)
  kappa_L substitution for left-censored values: 0, LOQ/2, LOQ/sqrt2, LOQ
  kappa_R multiplier for right-censored fecal coliform: 1, 2, 10
  gamma   ratio correction: 'log' (1+5log10 r for r>1) or 'none'
Objectives: O1 = number of the 53 verified reported classes reproduced (see common.UNVERIFIED);
            O2 = max |PI_c - PI_report| over the 18 PI values published in the 2024 report.
"""
import sys, json, itertools, time, numpy as np, pandas as pd
import os, common
os.makedirs(common.OUT, exist_ok=True)
d, fl = common.load()
lab = np.where(d.y == 1, 2, 1)
ver = common.verified(d)                          # 53 verified labels
NV = int(ver.sum())
BM = {'TDS': 1000, 'TSS': 50, 'NH3N': 0.2, 'BOD': 3, 'COD': 25, 'MBAS': 0.2, 'PO4': 0.2, 'MinyakLemak': 1,
      'NO3N': 10, 'Sulfida': 0.002, 'Pb': 0.03, 'FecalColi': 1000}
CAND = list(BM) + ['DO', 'pH']                     # 14 candidate parameters
K = len(CAND)
DO_R = ['sat', 'inv']; CL = {'0': 0.0, 'LOQ/2': 0.5, 'LOQ/sqrt2': 1 / np.sqrt(2), 'LOQ': 1.0}; CR = [1, 2, 10]; GAM = ['log', 'none']
REP = {'RAMBUT HULU': 2.39, 'RAMBUT HILIR': 2.41, 'CACABAN HULU': 5.06, 'CACABAN HILIR': 5.06, 'WALUH HULU': 5.05, 'WALUH HILIR': 2.17,
       'PALAWANGAN JEMBATAN KRAMAT HULU': 3.41, 'PALAWANGAN JEMBATAN KALI PALAWANGAN HILIR': 4.22, 'GUNG HULU': 1.70, 'GUNG HILIR': 5.11,
       'WADAS HULU': 5.11, 'WADAS HILIR': 4.92, 'GANGSA HULU': 5.63, 'GANGSA HILIR': 5.09, 'PEMALI HULU': 3.59, 'PEMALI HILIR': 5.08,
       'PAKIJANGAN HULU': 5.36, 'PAKIJANGAN HILIR': 4.46}
i24 = np.where(d.tahun == 2024)[0]
pi24 = np.array([REP[d.lokasi[i].replace('SUNGAI ', '').replace('SUNAI ', '')] for i in i24])


def do_sat(T): return 14.652 - 0.41022 * T + 0.0079910 * T ** 2 - 0.000077774 * T ** 3


def ratio_matrix(rdo, cl, cr, gam):
    """54 x 14 matrix of (corrected) sub-index ratios for one rule combination."""
    M = np.zeros((len(d), K))
    for j, p in enumerate(CAND):
        if p == 'DO':
            if rdo == 'inv': r = d.DO.to_numpy() / 4
            else:
                cs = do_sat(d.Suhu.to_numpy()); r = np.where(d.DO < cs, (cs - d.DO) / (cs - 4), 0.0)
        elif p == 'pH':
            ph = d.pH.to_numpy(); r = (ph - 7.5) / (np.where(ph <= 7.5, 6, 9) - 7.5)
        else:
            v = d[p].to_numpy().copy(); f = fl[p].to_numpy()
            v = np.where(f == -1, 2 * v * CL[cl], v)            # stored value is LOQ/2
            v = np.where(f == 1, v * cr, v)
            r = v / BM[p]
        M[:, j] = r
    if gam == 'log':
        M = np.where(M > 1, 1 + 5 * np.log10(np.maximum(M, 1e-12)), M)
    return M


# all non-empty subsets as a boolean mask matrix (16383 x 14)
S = np.array([[(m >> j) & 1 for j in range(K)] for m in range(1, 2 ** K)], dtype=bool)
size = S.sum(1)


def cls(pi): return np.where(pi <= 1, 0, np.where(pi <= 5, 1, np.where(pi <= 10, 2, 3)))


t0 = time.time(); rows = []; eq_pi = set(); eq_cls = set()
for rdo, cl, cr, gam in itertools.product(DO_R, CL, CR, GAM):
    M = ratio_matrix(rdo, cl, cr, gam)                       # 54 x 14
    mean = (M @ S.T) / size                                  # 54 x 16383
    mx = np.full(mean.shape, -np.inf)
    for j in range(K):
        mx = np.where(S[:, j][None, :], np.maximum(mx, M[:, [j]]), mx)
    PI = np.sqrt((mx ** 2 + mean ** 2) / 2)
    o1 = (cls(PI)[ver] == lab[ver][:, None]).sum(0)
    # equivalence classes: configurations with identical PI (or class) vectors over all samples
    for row in np.unique(np.round(PI.T, 9), axis=0): eq_pi.add(row.tobytes())
    for row in np.unique(cls(PI).T, axis=0): eq_cls.add(row.tobytes())
    o2 = np.abs(PI[i24] - pi24[:, None]).max(0)
    for k in range(S.shape[0]):
        rows.append((rdo, cl, cr, gam, k, int(o1[k]), float(o2[k])))
el = time.time() - t0
R = pd.DataFrame(rows, columns=['do', 'cl', 'cr', 'gam', 'sid', 'o1', 'o2'])
R['params'] = R.sid.map(lambda k: '+'.join(p for p, b in zip(CAND, S[k]) if b))
R['np'] = R.sid.map(lambda k: int(S[k].sum()))
N = len(R)
TOL = 0.025                                                   # reported PI are rounded to 2 decimals
full = R[R.o1 == NV]; num = R[R.o2 <= TOL]; both = R[(R.o1 == NV) & (R.o2 <= TOL)]
out = dict(N=N, NV=NV, n_eq_pi=len(eq_pi), n_eq_cls=len(eq_cls), seconds=el, n_subsets=int(S.shape[0]), n_rules=int(N / S.shape[0]),
           n_label=int(len(full)), n_num=int(len(num)), n_both=int(len(both)), tol=TOL,
           o1_dist={int(k): int(v) for k, v in R.o1.value_counts().sort_index().items()},
           both=both.sort_values('o2').drop(columns='sid').to_dict('records'),
           label_rules=full.groupby(['do', 'cl', 'cr', 'gam']).size().reset_index().values.tolist(),
           label_np=full.np.value_counts().sort_index().to_dict(),
           # standard (documented) specification, all 14 params
           std=R[(R.do == 'sat') & (R.cl == 'LOQ/2') & (R.cr == 1) & (R.gam == 'log') & (R.np == 14)][['o1', 'o2']].to_dict('records'),
           best_label_min_o2=float(full.o2.min()), label_o2_quant=[float(x) for x in full.o2.quantile([.5, .9]).values])
# which parameters are fixed across the numerically consistent set?
if len(both):
    Sb = S[both.sid.values]; out['incl_freq'] = {p: float(Sb[:, j].mean()) for j, p in enumerate(CAND)}
# coordinate ablation around the recovered config: change one dimension, report O1 and O2
P12 = ['TDS', 'TSS', 'NH3N', 'BOD', 'COD', 'MBAS', 'DO', 'MinyakLemak', 'NO3N', 'pH', 'Pb', 'FecalColi']
sid12 = int(np.where((S == np.array([p in P12 for p in CAND])).all(1))[0][0])
def get(do, cl, cr, gam, sid):
    x = R[(R.do == do) & (R.cl == cl) & (R.cr == cr) & (R.gam == gam) & (R.sid == sid)].iloc[0]; return [int(x.o1), round(float(x.o2), 3)]
def sid_of(ps): return int(np.where((S == np.array([p in ps for p in CAND])).all(1))[0][0])
abl = {'recovered': get('inv', 'LOQ', 1, 'log', sid12),
       'DO rule = saturation': get('sat', 'LOQ', 1, 'log', sid12),
       'left-censored = LOQ/2': get('inv', 'LOQ/2', 1, 'log', sid12),
       'left-censored = 0': get('inv', '0', 1, 'log', sid12),
       'FC right-censored x2': get('inv', 'LOQ', 2, 'log', sid12),
       'no ratio correction': get('inv', 'LOQ', 1, 'none', sid12),
       '+ sulfide': get('inv', 'LOQ', 1, 'log', sid_of(P12 + ['Sulfida'])),
       '+ phosphate': get('inv', 'LOQ', 1, 'log', sid_of(P12 + ['PO4'])),
       '+ sulfide + phosphate (documented set)': get('inv', 'LOQ', 1, 'log', sid_of(P12 + ['Sulfida', 'PO4'])),
       'documented specification': get('sat', 'LOQ/2', 1, 'log', sid_of(P12 + ['Sulfida', 'PO4']))}
out['ablation'] = abl
json.dump(out, open(os.path.join(common.OUT, 'spec.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in out.items() if k not in ('both',)}, indent=0)[:4000])
print('BOTH:'); [print(x) for x in out['both'][:20]]

# ---- hold-out check: identify with the 2024 PI values only, evaluate on 2023 and 2025 labels ----
iho = np.where((d.tahun != 2024).to_numpy() & ver)[0]       # 35 independent verified labels
ho = []
for _, x in num.iterrows():
    M = ratio_matrix(x.do, x.cl, x.cr, x.gam); s = S[x.sid]
    a = M[:, s]; PI = np.sqrt((a.max(1) ** 2 + a.mean(1) ** 2) / 2)
    ho.append(dict(do=x.do, cl=x.cl, cr=int(x.cr), gam=x.gam, params=x.params, o2=round(float(x.o2), 4),
                   lab24=int((cls(PI)[i24] == lab[i24]).sum()), hold=int((cls(PI)[iho] == lab[iho]).sum()), n_hold=len(iho),
                   missed=[f'{int(d.tahun[i])} {d.lokasi[i]} (PI {PI[i]:.2f})' for i in iho if cls(PI)[i] != lab[i]]))
out['holdout'] = ho
# search-space reduction by each evidence source
out['reduction'] = dict(all=N, labels=int(len(full)), pi2024=int(len(num)), both=int(len(both)))
json.dump(out, open(os.path.join(common.OUT, 'spec.json'), 'w'), indent=1)
for h in ho: print(h)
