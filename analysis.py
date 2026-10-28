import sys, json, numpy as np, pandas as pd
sys.path.insert(0, '/home/claude/exp')
import common
common.RAW = '/home/claude/data6/ds.xlsx'
from common import load
from sklearn.metrics import cohen_kappa_score
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

d, flags = load()
raw = pd.read_excel(common.RAW, 'Sheet1')
# nama sungai yang dinormalisasi
d['river'] = d['sungai'].str.title()
BM = {'TDS': 1000, 'TSS': 50, 'NH3N': 0.2, 'BOD': 3, 'COD': 25, 'MBAS': 0.2, 'PO4': 0.2, 'MinyakLemak': 1,
      'NO3N': 10, 'Sulfida': 0.002, 'Pb': 0.03, 'FecalColi': 1000}
IKA_P = ['pH', 'DO', 'BOD', 'COD', 'TSS', 'NO3N', 'PO4', 'FecalColi']
NAME = {'TDS': 'TDS', 'TSS': 'TSS', 'NH3N': 'NH3-N', 'BOD': 'BOD', 'COD': 'COD', 'MBAS': 'Detergent (MBAS)', 'PO4': 'Phosphate',
        'MinyakLemak': 'Oil and grease', 'NO3N': 'NO3-N', 'Sulfida': 'Sulfide', 'Pb': 'Dissolved Pb', 'FecalColi': 'Fecal coliform',
        'DO': 'DO', 'pH': 'pH', 'Suhu': 'Temperature'}
UNIT = {'TDS': 'mg/L', 'TSS': 'mg/L', 'NH3N': 'mg/L', 'BOD': 'mg/L', 'COD': 'mg/L', 'MBAS': 'mg/L', 'PO4': 'mg/L', 'MinyakLemak': 'mg/L',
        'NO3N': 'mg/L', 'Sulfida': 'mg/L', 'Pb': 'mg/L', 'FecalColi': 'MPN/100 mL', 'DO': 'mg/L', 'pH': '-', 'Suhu': '°C'}
STD = dict(BM, DO='≥ 4', pH='6–9', Suhu='dev. 3')
MODES = ['half', 'lod', 'sqrt2', 'zero']
MLAB = {'half': 'LOQ/2', 'lod': 'LOQ', 'sqrt2': 'LOQ/√2', 'zero': '0'}
CLS = ['Good', 'Lightly polluted', 'Moderately polluted', 'Heavily polluted']


def do_sat(T): return 14.652 - 0.41022 * T + 0.0079910 * T ** 2 - 0.000077774 * T ** 3


def ratios(i, mode, params):
    row = d.loc[i]; r = {}
    for p in params:
        if p == 'DO':
            cs = do_sat(row['Suhu']); r['DO'] = (cs - row['DO']) / (cs - 4) if row['DO'] < cs else 0.0
        elif p == 'pH':
            ph = row['pH']; r['pH'] = (ph - 7.5) / ((6 if ph <= 7.5 else 9) - 7.5)
        else:
            v = row[p]
            if flags.loc[i, p] == -1:
                lod = v * 2; v = {'half': lod / 2, 'lod': lod, 'sqrt2': lod / np.sqrt(2), 'zero': 0.0}[mode]
            r[p] = v / BM[p]
    return {k: (1 + 5 * np.log10(v) if v > 1 else v) for k, v in r.items()}


def ip(i, mode, params):
    rr = ratios(i, mode, params); a = np.array(list(rr.values()))
    return float(np.sqrt((a.max() ** 2 + a.mean() ** 2) / 2)), max(rr, key=rr.get), rr


def cls(x): return 0 if x <= 1 else 1 if x <= 5 else 2 if x <= 10 else 3
def ika(c): c = np.asarray(c); return float(sum((c == k).mean() * w for k, w in zip(range(4), [70, 50, 30, 10])))
def ika_cat(x): return 'Very good' if x >= 90 else 'Good' if x >= 70 else 'Moderate' if x >= 50 else 'Poor' if x >= 25 else 'Very poor'


ALL = list(BM) + ['DO', 'pH']
R = {}
# ---------- deskriptif ----------
desc = []
for p in ['Suhu', 'pH', 'DO', 'TDS', 'TSS', 'BOD', 'COD', 'NH3N', 'NO3N', 'PO4', 'Sulfida', 'MBAS', 'MinyakLemak', 'Pb', 'FecalColi']:
    v = d[p]; f = flags[p]
    lods = sorted(set((v[f == -1] * 2).round(3)))
    if p in BM:
        exc = (v[f != -1] > BM[p]).sum()
        # tersensor kiri dengan LOD > baku mutu: status tak tentu
        lod_gt = int(((f == -1) & (v * 2 > BM[p])).sum())
    elif p == 'DO': exc = (v < 4).sum(); lod_gt = 0
    elif p == 'pH': exc = ((v < 6) | (v > 9)).sum(); lod_gt = 0
    else: exc = None; lod_gt = 0
    desc.append(dict(p=NAME[p], unit=UNIT[p], std=STD[p], cens=int((f == -1).sum()), rcens=int((f == 1).sum()),
                     lods=', '.join(f'{x:g}' for x in lods), med=float(v.median()), mn=float(v.min()), mx=float(v.max()),
                     exc=None if exc is None else int(exc), lod_gt=lod_gt, ika=p in IKA_P))
R['desc'] = desc

# ---------- status per skenario ----------
lab = np.where(d.y == 1, 2, 1)
for mode in MODES:
    for nm, params in [('IP15', ALL), ('IKA8', IKA_P)]:
        res = [ip(i, mode, params) for i in d.index]
        d[f'{nm}_{mode}'] = [r[0] for r in res]; d[f'{nm}c_{mode}'] = [cls(r[0]) for r in res]
        d[f'{nm}arg_{mode}'] = [r[1] for r in res]
sc = {}
for nm in ['IP15', 'IKA8']:
    for mode in MODES:
        c = d[f'{nm}c_{mode}']
        sc[f'{nm}_{mode}'] = dict(dist=[int((c == k).sum()) for k in range(4)], ika=ika(c),
                                 changed=int((c != d[f'{nm}c_half']).sum()),
                                 kappa=float(cohen_kappa_score(lab, c)), agree=float((lab == c).mean()),
                                 medIP=float(d[f'{nm}_{mode}'].median()))
R['scen'] = sc
# determinan
R['det_IP15'] = {NAME[k]: int(v) for k, v in d['IP15arg_half'].value_counts().items()}
R['det_IKA8'] = {NAME[k]: int(v) for k, v in d['IKA8arg_half'].value_counts().items()}
R['det_IP15_zero'] = {NAME[k]: int(v) for k, v in d['IP15arg_zero'].value_counts().items()}
# rata-rata rasio terkoreksi per parameter (kontribusi)
rr = pd.DataFrame([ip(i, 'half', ALL)[2] for i in d.index])
R['ratio_mean'] = {NAME[k]: float(v) for k, v in rr.mean().sort_values(ascending=False).items()}
R['ratio_gt1'] = {NAME[k]: int((rr[k] > 1).sum()) for k in rr.columns}
# leave-one-parameter-out
lopo = {}
for p in ALL:
    keep = [q for q in ALL if q != p]
    c = np.array([cls(ip(i, 'half', keep)[0]) for i in d.index])
    lopo[NAME[p]] = dict(changed=int((c != d['IP15c_half']).sum()), dist=[int((c == k).sum()) for k in range(4)])
R['lopo'] = dict(sorted(lopo.items(), key=lambda kv: -kv[1]['changed']))
# coliform tersensor kanan: >16000 dianggap 16000 (batas bawah); kepekaan bila 2x
fc = flags['FecalColi'] == 1
R['fc_right'] = int(fc.sum())
d2 = d.copy()
# per tahun & per sungai (skenario utama)
by = []
for y, g in d.groupby('tahun'):
    by.append(dict(year=int(y), n=len(g), IP15=[int((g['IP15c_half'] == k).sum()) for k in range(4)], ika15=ika(g['IP15c_half']),
                   IKA8=[int((g['IKA8c_half'] == k).sum()) for k in range(4)], ika8=ika(g['IKA8c_half']), ika8_cat=ika_cat(ika(g['IKA8c_half'])),
                   medIP=float(g['IP15_half'].median())))
R['year'] = by
br = []
for r, g in d.groupby('river'):
    br.append(dict(river=r, n=len(g), years=', '.join(str(x) for x in sorted(g.tahun.unique())),
                   IP15=[int((g['IP15c_half'] == k).sum()) for k in range(4)], medIP=float(g['IP15_half'].median()),
                   ika8=ika(g['IKA8c_half']), lab=[int((g.y == 0).sum()), int((g.y == 1).sum())],
                   det=NAME[g['IP15arg_half'].mode()[0]]))
R['river'] = sorted(br, key=lambda x: -x['n'])
R['ika_all'] = dict(ika8=ika(d['IKA8c_half']), cat=ika_cat(ika(d['IKA8c_half'])), ika15=ika(d['IP15c_half']))
R['n'] = len(d); R['n_rivers'] = d.river.nunique(); R['n_points'] = d.lokasi.nunique()
R['label'] = dict(ringan=int((d.y == 0).sum()), sedang=int((d.y == 1).sum()))
R['xtab'] = pd.crosstab(pd.Series(lab, name='lab'), pd.Series(d['IP15c_half'], name='ip')).to_dict()
# bootstrap CI untuk proporsi cemar sedang+ (IP15 utama)
rng = np.random.default_rng(1)
c = d['IP15c_half'].to_numpy()
bs = [((c[rng.integers(0, len(c), len(c))]) >= 2).mean() for _ in range(5000)]
R['p_mod'] = dict(p=float((c >= 2).mean()), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)))
bs = [ika(c8[rng.integers(0, len(c8), len(c8))]) for c8 in [d['IKA8c_half'].to_numpy()] for _ in range(5000)]
R['ika_ci'] = [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))]
json.dump(R, open('/home/claude/art1/res.json', 'w'), indent=1, ensure_ascii=False)
d.to_csv('/home/claude/art1/status_per_sample.csv', index=False)

# ---------- gambar ----------
plt.rcParams.update({'font.family': 'serif', 'font.serif': ['Times New Roman', 'DejaVu Serif'], 'font.size': 9})
SEQ = ['#D6E4F0', '#8FB4D9', '#3D78B5', '#0B3C73']  # ordinal: good -> heavy
INK = '#333333'
# Fig 2: determinan (dua set parameter)
fig, ax = plt.subplots(figsize=(5.6, 2.6))
keys = ['Sulfide', 'Fecal coliform', 'BOD', 'DO', 'Phosphate', 'TSS'][::-1]
a = [R['det_IP15'].get(k, 0) for k in keys]; b = [R['det_IKA8'].get(k, 0) for k in keys]
y = np.arange(len(keys)); h = 0.36
ax.barh(y + h / 2 + 0.02, a, height=h, color='#3D78B5', label='PI, 14 parameters')
ax.barh(y - h / 2 - 0.02, b, height=h, color='#A9C4E2', label='IKA parameter set (8)')
for yy, v in zip(y + h / 2, a):
    if v: ax.text(v + 0.4, yy, str(v), va='center', fontsize=7.5, color=INK)
for yy, v in zip(y - h / 2, b):
    if v: ax.text(v + 0.4, yy, str(v), va='center', fontsize=7.5, color=INK)
ax.set_yticks(y, keys); ax.set_xlabel('Samples in which the parameter has the maximum corrected ratio')
for s_ in ['top', 'right']: ax.spines[s_].set_visible(False)
ax.grid(axis='x', alpha=.25); ax.legend(frameon=False, fontsize=7.5, loc='lower right')
fig.tight_layout(); fig.savefig('/home/claude/art1/fig_det.png', dpi=300); plt.close(fig)
# Fig 3: skenario LOD (stacked, IP15 dan IKA8)
fig, axs = plt.subplots(1, 2, figsize=(6.0, 2.6), sharey=True)
for ax, nm, tt in zip(axs, ['IP15', 'IKA8'], ['(a) PI, 14 parameters', '(b) IKA parameter set']):
    left = np.zeros(4)
    for k in range(4):
        vals = np.array([sc[f'{nm}_{m}']['dist'][k] for m in MODES])
        ax.barh([MLAB[m] for m in MODES], vals, left=left, color=SEQ[k], edgecolor='white', linewidth=1.5, label=CLS[k], height=0.6)
        for j, (l, v) in enumerate(zip(left, vals)):
            if v >= 4: ax.text(l + v / 2, j, str(v), ha='center', va='center', fontsize=7.5, color='white' if k >= 2 else INK)
        left += vals
    ax.set_title(tt, fontsize=9); ax.set_xlim(0, 52); ax.set_xlabel('Number of samples')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
axs[0].set_ylabel('Substitution of values < LOQ')
h, l = axs[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=4, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.08, 1, 1)); fig.savefig('/home/claude/art1/fig_lod.png', dpi=300); plt.close(fig)
print(json.dumps({k: R[k] for k in ['scen', 'det_IP15', 'det_IKA8', 'det_IP15_zero', 'lopo', 'year', 'ika_all', 'n_rivers', 'n_points', 'p_mod', 'ika_ci', 'fc_right', 'xtab', 'ratio_gt1']}, indent=0, ensure_ascii=False))
for x in desc: print(x)
for x in R['river']: print(x)

# ---------- kepekaan coliform tersensor kanan ----------
fcr = {}
base = d['FecalColi'].copy()
for mult in [1, 2, 5, 10]:
    d['FecalColi'] = np.where(fc, base * mult, base)
    c15 = np.array([cls(ip(i, 'half', ALL)[0]) for i in d.index]); c8 = np.array([cls(ip(i, 'half', IKA_P)[0]) for i in d.index])
    fcr[mult] = dict(IP15=[int((c15 == k).sum()) for k in range(4)], IKA8=[int((c8 == k).sum()) for k in range(4)],
                     ika8=ika(c8), ch15=int((c15 != d['IP15c_half']).sum()), ch8=int((c8 != d['IKA8c_half']).sum()))
d['FecalColi'] = base
R['fc_sens'] = fcr
# sulfida tersensor dengan LOD > baku mutu: status bila sulfida dikeluarkan untuk sampel tersebut
R['sulf_cens_idx'] = [int(i) for i in d.index[flags['Sulfida'] == -1]]
R['sulf_cens_status'] = [int(x) for x in d.loc[flags['Sulfida'] == -1, 'IP15c_half']]
json.dump(R, open('/home/claude/art1/res.json', 'w'), indent=1, ensure_ascii=False)
print('FC', fcr); print('sulf cens', R['sulf_cens_status'], d.loc[flags['Sulfida'] == -1, 'IP15arg_half'].tolist())
