#!/usr/bin/env python3
"""
F07_ulc_covariance.py — Covariance structure UNDER ULC conditioning (c_M = 1)

Key questions:
  1. Does Cov(c_j, g_i h_{j-i} | c_M=1) = 1/8 still hold?
  2. What is the deviation η_{j,i} = Cov_ULC - 1/8?
  3. Does P(total_j odd | c_M=1) = 1/2 still hold?
  4. What is E[carry_j | c_M=1] compared to (j-1)/4?

Uses the MSB-first carry chain (same convention as Lemma D):
  g_0 = h_0 = 1 (MSBs), g_i, h_i ~ Ber(1/2) for i >= 1
  conv_j = Σ_{i=0}^{j} g_i h_{j-i}
  carry_0 = 0, carry_{j+1} = floor((conv_j + carry_j) / 2)

The ULC condition: carry_M = 1 where M is the highest nonzero carry position.
For semiprimes with d-bit factors: M = D-1 or lower, and carry_M = 1.

We compute for the MSB-first chain of length K+1 (positions 0..K).
"""

import sys
from fractions import Fraction

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def exact_ulc_analysis(K):
    """Exact analysis comparing unconditional vs ULC-conditioned covariances.
    MSB-first convention: g[0]=h[0]=1, g[i],h[i] ~ Ber(1/2) for i=1..K.
    """
    d = K + 1
    n_configs = 1 << (2 * K)

    n_ulc = 0

    # Per-position accumulators (positions j = 0..K)
    unc_Sc = [0] * d
    unc_Sv = [0] * d
    unc_Scv = [0] * d
    unc_Sc2 = [0] * d
    unc_Podd = [0] * d

    ulc_Sc = [0] * d
    ulc_Sv = [0] * d
    ulc_Scv = [0] * d
    ulc_Podd = [0] * d

    # Individual off-diagonal term accumulators: key (j, i) -> [S_cp, S_c, S_p]
    unc_term = {}
    ulc_term = {}

    for cfg in range(n_configs):
        g = [1] + [(cfg >> (2*K - 1 - idx)) & 1 for idx in range(K)]
        h = [1] + [(cfg >> (K - 1 - idx)) & 1 for idx in range(K)]

        carries = [0] * (d + 1)
        convs = [0] * d
        for j in range(d):
            s = 0
            for i in range(j + 1):
                if i < d and j - i < d:
                    s += g[i] * h[j - i]
            convs[j] = s
            carries[j + 1] = (s + carries[j]) >> 1

        # Find M_top: highest position with nonzero carry
        M_top = 0
        for pos in range(d, 0, -1):
            if carries[pos] > 0:
                M_top = pos
                break
        is_ulc = (M_top > 0 and carries[M_top] == 1)

        # Accumulate unconditional stats
        for j in range(d):
            c = carries[j]
            v = convs[j]
            total = c + v
            unc_Sc[j] += c
            unc_Sv[j] += v
            unc_Scv[j] += c * v
            unc_Sc2[j] += c * c
            if total % 2 == 1:
                unc_Podd[j] += 1

            for i in range(1, j):
                ji = j - i
                if i < d and ji < d and i != ji:
                    p = g[i] * h[ji]
                    key = (j, i)
                    if key not in unc_term:
                        unc_term[key] = [0, 0, 0]
                    unc_term[key][0] += c * p
                    unc_term[key][1] += c
                    unc_term[key][2] += p

        # Accumulate ULC-conditioned stats
        if is_ulc:
            n_ulc += 1
            for j in range(d):
                c = carries[j]
                v = convs[j]
                total = c + v
                ulc_Sc[j] += c
                ulc_Sv[j] += v
                ulc_Scv[j] += c * v
                if total % 2 == 1:
                    ulc_Podd[j] += 1

                for i in range(1, j):
                    ji = j - i
                    if i < d and ji < d and i != ji:
                        p = g[i] * h[ji]
                        key = (j, i)
                        if key not in ulc_term:
                            ulc_term[key] = [0, 0, 0]
                        ulc_term[key][0] += c * p
                        ulc_term[key][1] += c
                        ulc_term[key][2] += p

    N = n_configs
    pr(f"\n  ULC valid: {n_ulc}/{N} = {n_ulc/N:.6f}")

    # === PARITY LEMMA ===
    pr(f"\n  === PARITY LEMMA: P(total_j odd) ===")
    pr(f"  {'j':>4}  {'Uncond':>12}  {'=1/2?':>6}  {'ULC':>12}  {'Δ_ULC':>12}")
    for j in range(1, d):
        p_unc = Fraction(unc_Podd[j], N)
        p_ulc = Fraction(ulc_Podd[j], n_ulc) if n_ulc > 0 else Fraction(0)
        is_half = (p_unc == Fraction(1, 2))
        delta = float(p_ulc) - 0.5
        pr(f"  {j:4d}  {float(p_unc):12.8f}  {'✓' if is_half else '✗':>6}  "
           f"{float(p_ulc):12.8f}  {delta:+12.8f}")

    # === E[carry_j] ===
    pr(f"\n  === E[carry_j] ===")
    pr(f"  {'j':>4}  {'Uncond':>12}  {'(j-1)/4':>10}  {'ULC':>12}  {'Δ_ULC':>12}")
    for j in range(1, d):
        ec_unc = Fraction(unc_Sc[j], N)
        ec_ulc = Fraction(ulc_Sc[j], n_ulc) if n_ulc > 0 else Fraction(0)
        target = Fraction(j - 1, 4)
        match_unc = "✓" if ec_unc == target else ""
        delta = float(ec_ulc) - float(target)
        pr(f"  {j:4d}  {float(ec_unc):12.8f} {match_unc:1s}  "
           f"{float(target):10.6f}  {float(ec_ulc):12.8f}  {delta:+12.8f}")

    # === TOTAL Cov(c_j, v_j) ===
    pr(f"\n  === Cov(c_j, v_j) ===")
    pr(f"  {'j':>4}  {'Uncond':>14}  {'(j-1)/8':>10}  {'ULC':>14}  {'Δ_ULC':>12}  "
       f"{'parity':>6}")
    for j in range(2, d):
        ec = Fraction(unc_Sc[j], N)
        ev = Fraction(unc_Sv[j], N)
        ecv = Fraction(unc_Scv[j], N)
        cov_unc = ecv - ec * ev

        ec_u = Fraction(ulc_Sc[j], n_ulc)
        ev_u = Fraction(ulc_Sv[j], n_ulc)
        ecv_u = Fraction(ulc_Scv[j], n_ulc)
        cov_ulc = ecv_u - ec_u * ev_u

        target = Fraction(j - 1, 8)
        delta = float(cov_ulc) - float(target)
        parity = "odd" if j % 2 == 1 else "even"

        unc_match = ""
        if j % 2 == 1 and cov_unc == target:
            unc_match = "✓"
        pr(f"  {j:4d}  {float(cov_unc):14.10f}{unc_match:1s}  {float(target):10.6f}  "
           f"{float(cov_ulc):14.10f}  {delta:+12.8f}  {parity:>6}")

    # === INDIVIDUAL OFF-DIAGONAL TERMS ===
    pr(f"\n  === Individual Cov(c_j, g_i h_{{j-i}}) ===")
    pr(f"  {'(j,i)':>8}  {'Uncond':>14}  {'=1/8?':>6}  {'ULC':>14}  "
       f"{'Δ=ULC-1/8':>14}  {'η/⅛':>10}")

    positions = sorted(set(k[0] for k in ulc_term.keys()))
    for j in positions:
        if j < 3:
            continue
        parity = "odd" if j % 2 == 1 else "even"
        pr(f"\n  j={j} ({parity}):")
        for i in range(1, j):
            ji = j - i
            if ji <= 0 or ji >= d or i >= d or i == ji or i > ji:
                continue
            key = (j, i)
            if key not in ulc_term:
                continue

            # Unconditional
            ud = unc_term.get(key, [0, 0, 0])
            cov_unc = Fraction(ud[0], N) - Fraction(ud[1], N) * Fraction(ud[2], N)

            # ULC
            cd = ulc_term[key]
            cov_ulc = Fraction(cd[0], n_ulc) - Fraction(cd[1], n_ulc) * Fraction(cd[2], n_ulc)

            unc_is_eighth = (cov_unc == Fraction(1, 8))
            delta = float(cov_ulc) - 0.125
            eta_rel = delta / 0.125

            pr(f"  ({j:2d},{i:2d})  {float(cov_unc):14.10f}  "
               f"{'✓' if unc_is_eighth else '✗':>6}  {float(cov_ulc):14.10f}  "
               f"{delta:+14.10f}  {eta_rel:+10.6f}")

    # === PATTERN ANALYSIS ===
    pr(f"\n  === PATTERN: η(j,i) = Cov_ULC(j,i) - 1/8 ===")
    pr(f"  For each j, how does η depend on the 'distance from boundary'?")
    for j in positions:
        if j < 5:
            continue
        parity = "odd" if j % 2 == 1 else "even"
        pr(f"\n  j={j} ({parity}), distance from top = {K - j}:")
        for i in range(1, j):
            ji = j - i
            if ji <= 0 or ji >= d or i >= d or i == ji or i > ji:
                continue
            key = (j, i)
            if key not in ulc_term:
                continue
            cd = ulc_term[key]
            cov_ulc = Fraction(cd[0], n_ulc) - Fraction(cd[1], n_ulc) * Fraction(cd[2], n_ulc)
            delta = float(cov_ulc) - 0.125
            dist_i = min(i, ji)
            pr(f"    i={i}, j-i={ji}, min_dist={dist_i}: η = {delta:+.10f}")


def main():
    pr("=" * 78)
    pr("  ULC-CONDITIONED COVARIANCE ANALYSIS (MSB-first convention)")
    pr("=" * 78)

    for K in [8, 10]:
        pr(f"\n{'═' * 78}")
        pr(f"  K = {K}  (d={K+1}, positions j=0..{K}, configs={1<<(2*K):,})")
        pr(f"{'═' * 78}")
        exact_ulc_analysis(K)

    pr("\n" + "=" * 78)
    pr("  END")
    pr("=" * 78)


if __name__ == '__main__':
    main()
