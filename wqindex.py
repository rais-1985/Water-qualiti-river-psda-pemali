"""wqindex: Pollution Index (PI) and IKA computation with censored measurements and configurable specifications.

Each measurement is a Censored value (value, flag, limit). A Spec describes one index configuration;
pi() evaluates it. The default Spec is the documented specification (Kepmen LH 115/2003, class II limits
of PP 22/2021); LAB_SPEC is the specification recovered from laboratory reports.
"""
from dataclasses import dataclass, field
import math

LIMITS = {'TDS': 1000, 'TSS': 50, 'NH3N': 0.2, 'BOD': 3, 'COD': 25, 'MBAS': 0.2, 'PO4': 0.2, 'MinyakLemak': 1,
          'NO3N': 10, 'Sulfida': 0.002, 'Pb': 0.03, 'FecalColi': 1000}
ALL14 = tuple(LIMITS) + ('DO', 'pH')
CLASSES = ('good', 'light', 'moderate', 'heavy')


@dataclass(frozen=True)
class Censored:
    """A laboratory result: flag -1 = '<limit' (left-censored), +1 = '>limit' (right-censored), 0 = quantified."""
    value: float
    flag: int = 0

    @staticmethod
    def parse(x):
        if isinstance(x, (int, float)) and not (isinstance(x, float) and math.isnan(x)):
            return Censored(float(x), 0)
        s = str(x).strip().replace(' ', '').replace(',', '.')
        if s.startswith('<'): return Censored(float(s[1:]), -1)
        if s.startswith('>'): return Censored(float(s[1:]), 1)
        return Censored(float(s), 0)

    def resolve(self, left=0.5, right=1.0):
        """Numeric value used in the index: left-censored -> left*limit, right-censored -> right*limit."""
        if self.flag == -1: return left * self.value
        if self.flag == 1: return right * self.value
        return self.value


@dataclass(frozen=True)
class Spec:
    params: tuple = ALL14
    do_rule: str = 'sat'        # 'sat' (saturation-based) or 'inv' (C/L = DO/4)
    left: float = 0.5           # substitution factor for left-censored values (0, 0.5, 1/sqrt2, 1)
    right: float = 1.0          # multiplier for right-censored values
    correct: bool = True        # 1 + 5 log10(r) for r > 1


LAB_SPEC = Spec(params=('TDS', 'TSS', 'NH3N', 'BOD', 'COD', 'MBAS', 'DO', 'MinyakLemak', 'NO3N', 'pH', 'Pb', 'FecalColi'),
                do_rule='inv', left=1.0, right=1.0, correct=True)


def do_sat(t):
    return 14.652 - 0.41022 * t + 0.0079910 * t ** 2 - 0.000077774 * t ** 3


def ratios(sample, spec=Spec()):
    """sample: dict name -> Censored (plus 'Suhu' for the DO rule). Returns dict of (corrected) ratios."""
    r = {}
    for p in spec.params:
        if p == 'DO':
            do = sample['DO'].resolve()
            if spec.do_rule == 'inv':
                r[p] = do / 4
            else:
                cs = do_sat(sample['Suhu'].resolve())
                r[p] = (cs - do) / (cs - 4) if do < cs else 0.0
        elif p == 'pH':
            ph = sample['pH'].resolve(); r[p] = (ph - 7.5) / ((6 if ph <= 7.5 else 9) - 7.5)
        else:
            r[p] = sample[p].resolve(spec.left, spec.right) / LIMITS[p]
    if spec.correct:
        r = {k: (1 + 5 * math.log10(v) if v > 1 else v) for k, v in r.items()}
    return r


def pi(sample, spec=Spec()):
    a = list(ratios(sample, spec).values())
    return math.sqrt((max(a) ** 2 + (sum(a) / len(a)) ** 2) / 2)


def pi_class(x):
    return 0 if x <= 1 else 1 if x <= 5 else 2 if x <= 10 else 3


def ika(classes):
    w = (70, 50, 30, 10); n = len(classes)
    return sum(w[c] for c in classes) / n
