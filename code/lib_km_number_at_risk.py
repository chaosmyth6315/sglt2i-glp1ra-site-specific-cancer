"""Number at risk recovered from a TriNetX Kaplan–Meier export (Greenwood back-calculation).

The export (Outcome_N_Result_b_KM_graph.csv) gives, per day, the survival probability S(t) and its 95% CI
for each arm. The CI is symmetric on the log(−log S) scale (checked on the exports), i.e.
    log(−log L) − log(−log S) = z · sqrt(G(t)) / |ln S(t)|,   G(t) = Σ d_i / (n_i (n_i − d_i))   (Greenwood).
On an event day, a = d/n = 1 − S(t)/S(t−) and ΔG = d / (n (n − d)) = a / (n − d), so
    n = a / (ΔG (1 − a)),   d = a n.
The curve starts at the first day of the outcome window with S = 1 and G = 0 before it, so the first event day
is recovered too. Between event days n changes only through censoring, so the number at risk at a time point t
is taken as n on the first event day ≥ t (a lower bound that differs from the exact value by the censorings in
between). When no event occurs on or after t (the end of the window), n on the last event day before t is used
(an upper bound).

read_km(path) -> (days, S1, L1, U1, S2, L2, U2); risk_table(path, times) -> {arm: {t: n}}
Usage (check): python3 -I lib_km_number_at_risk.py <KM_graph.csv> [times...]
"""
import csv
import math
import sys
from decimal import Decimal, ROUND_HALF_UP

Z = 1.959963984540054


def read_km(path):
    rows = [r for r in csv.reader(open(path, encoding="utf-8-sig", newline="")) if r and r[0].strip().isdigit()]
    cols = list(zip(*[(int(r[0]),) + tuple(float(x) for x in r[1:7]) for r in rows]))
    return cols  # days, S1, L1, U1, S2, L2, U2


def at_risk_series(days, S, L, U):
    """[(day, n, d)] on event days."""
    out, prevS, prevG = [], 1.0, 0.0   # state before the first day of the window
    for t, s, lo, hi in zip(days, S, L, U):
        if s <= 0 or s >= 1 or lo <= 0:
            prevS, prevG = s, None
            continue
        h = math.log(-math.log(lo)) - math.log(-math.log(s))
        G = (h * abs(math.log(s)) / Z) ** 2
        if prevS is not None and prevG is not None and s < prevS:
            a = 1 - s / prevS
            dG = G - prevG
            if dG > 0 and a < 1:
                n = a / (dG * (1 - a))
                out.append((t, n, a * n))
        prevS, prevG = s, G
    return out


def risk_table(path, times):
    days, S1, L1, U1, S2, L2, U2 = read_km(path)
    res = {}
    for arm, (S, L, U) in (("S", (S1, L1, U1)), ("G", (S2, L2, U2))):
        ser = at_risk_series(days, S, L, U)
        res[arm] = {}
        for t in times:
            nxt = next((n for d, n, _ in ser if d >= t), None)
            if nxt is None and ser:  # no event on or after t (end of window): nearest preceding event day
                nxt = [n for d, n, _ in ser if d < t][-1]
            res[arm][t] = None if nxt is None else int(Decimal(str(nxt)).quantize(Decimal("1"), rounding=ROUND_HALF_UP))
    return res


if __name__ == "__main__":
    times = [int(x) for x in sys.argv[2:]] or [365, 730, 1095, 1460, 1825, 2190, 2555, 2920, 3285, 3650]
    rt = risk_table(sys.argv[1], times)
    days, S1, L1, U1, *_ = read_km(sys.argv[1])
    ser = at_risk_series(days, S1, L1, U1)
    print("events recovered (arm 1): sum d =", round(sum(d for _, _, d in ser)), "over", len(ser), "event days")
    for arm in ("S", "G"):
        print(arm, [rt[arm][t] for t in times])
