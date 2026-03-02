#!/usr/bin/env python3
"""
F02_covariance_extended.py — Extended exact covariance computation

Computes Cov(carry_j, conv_j) exactly for j=2..K using integer accumulation.
Tests the conjectured formula Cov = (j-1)/8 for odd j, and searches for
the even-j pattern with data up to j=14.

Also computes sensitivity weights conditioned on c_top = 1.
"""

import sys
import math
from fractions import Fraction

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def exact_covariance(K):
    """Exact Cov(carry_j, conv_j) for j=0..K via integer accumulation."""
    n_bits = 2 * K
    n_configs = 1 << n_bits

    sum_c = [0] * (K + 1)
    sum_v = [0] * (K + 1)
    sum_cv = [0] * (K + 1)
    sum_cc = [0] * (K + 1)
    sum_vv = [0] * (K + 1)

    report_interval = max(1, n_configs // 20)

    for cfg in range(n_configs):
        if cfg % report_interval == 0 and cfg > 0:
            pr(f"    K={K}: {cfg}/{n_configs} ({100*cfg//n_configs}%)")

        g = [0] * (K + 1)
        h = [0] * (K + 1)
        g[0] = 1
        h[0] = 1
        for i in range(K):
            g[i + 1] = (cfg >> (2 * K - 1 - i)) & 1
            h[i + 1] = (cfg >> (K - 1 - i)) & 1

        carry = 0
        for j in range(K + 1):
            conv_j = 0
            for i in range(j + 1):
                if i <= K and (j - i) <= K:
                    conv_j += g[i] * h[j - i]

            sum_c[j] += carry
            sum_v[j] += conv_j
            sum_cv[j] += carry * conv_j
            sum_cc[j] += carry * carry
            sum_vv[j] += conv_j * conv_j

            carry = (conv_j + carry) >> 1

    N = n_configs
    results = []
    for j in range(K + 1):
        Ec = Fraction(sum_c[j], N)
        Ev = Fraction(sum_v[j], N)
        Ecv = Fraction(sum_cv[j], N)
        Ecc = Fraction(sum_cc[j], N)
        Evv = Fraction(sum_vv[j], N)
        cov = Ecv - Ec * Ev
        var_c = Ecc - Ec * Ec
        var_v = Evv - Ev * Ev
        results.append((j, Ec, Ev, cov, var_c, var_v))

    return results


def sensitivity_analysis(K):
    """Compute E[c_{top-k} | c_top=1] for a K+1 bit multiplication."""
    n_bits = 2 * K
    n_configs = 1 << n_bits
    d = K + 1
    D_max = 2 * d

    sum_c_top_k = [0] * (K + 2)
    n_valid = 0

    report_interval = max(1, n_configs // 10)

    for cfg in range(n_configs):
        if cfg % report_interval == 0 and cfg > 0:
            pr(f"    Sensitivity K={K}: {cfg}/{n_configs} ({100*cfg//n_configs}%)")

        g = [0] * d
        h = [0] * d
        g[0] = 1
        h[0] = 1
        for i in range(K):
            g[i + 1] = (cfg >> (2 * K - 1 - i)) & 1
            h[i + 1] = (cfg >> (K - 1 - i)) & 1

        conv = [0] * D_max
        for i in range(d):
            for jj in range(d):
                pos = i + jj
                if pos < D_max:
                    conv[pos] += g[d - 1 - i] * h[d - 1 - jj]

        carries = [0] * (D_max + 1)
        for pos in range(D_max):
            carries[pos + 1] = (conv[pos] + carries[pos]) >> 1

        M_top = D_max
        while M_top > 0 and carries[M_top] == 0:
            M_top -= 1

        if carries[M_top] != 1:
            continue

        n_valid += 1
        for k in range(min(K + 2, M_top + 1)):
            sum_c_top_k[k] += carries[M_top - k]

    return n_valid, n_configs, sum_c_top_k


def main():
    pr("=" * 78)
    pr("  EXTENDED COVARIANCE ANALYSIS")
    pr("=" * 78)

    # Part 1: Covariance for K=8,10,12
    all_results = {}
    for K in [8, 10, 12]:
        pr(f"\n{'═' * 78}")
        pr(f"  K={K} — {1 << (2*K):,} configurations")
        pr(f"{'═' * 78}")

        results = exact_covariance(K)
        all_results[K] = results

        pr(f"\n  {'j':>3}  {'E[c]':>10}  {'E[v]':>10}  {'Cov(c,v)':>18}  "
           f"{'Cov float':>12}  {'(j-1)/8':>10}  {'Δ':>14}")

        for j, Ec, Ev, cov, var_c, var_v in results:
            if j < 2:
                continue
            target = Fraction(j - 1, 8)
            delta = cov - target
            flag = "  ✓ EXACT" if delta == 0 else ""
            pr(f"  {j:3d}  {str(Ec):>10}  {str(Ev):>10}  {str(cov):>18}  "
               f"{float(cov):12.8f}  {float(target):10.6f}  {str(delta):>14}{flag}")

    # Part 2: Pattern analysis with K=12 data
    pr(f"\n{'═' * 78}")
    pr(f"  PATTERN ANALYSIS — EVEN j")
    pr(f"{'═' * 78}")

    results_12 = all_results[12]
    pr(f"\n  Even-j covariances and correction from (j-1)/8:")
    pr(f"  {'j':>3}  {'Cov':>18}  {'(j-1)/8':>12}  {'correction':>18}  "
       f"{'corr float':>14}  {'1/2^?':>10}")

    corrections = []
    for j, Ec, Ev, cov, _, _ in results_12:
        if j < 2 or j % 2 != 0:
            continue
        target = Fraction(j - 1, 8)
        corr = cov - target
        corrections.append((j, cov, corr))
        if corr != 0 and corr.numerator > 0:
            log2_inv = math.log2(float(corr.denominator / corr.numerator))
        else:
            log2_inv = float('nan')
        pr(f"  {j:3d}  {str(cov):>18}  {str(target):>12}  {str(corr):>18}  "
           f"{float(corr):14.10f}  2^{-log2_inv:.2f}")

    pr(f"\n  Testing alternative formulas for even j:")
    for j, cov, corr in corrections:
        k = j // 2
        candidates = [
            ("(2k-1)!!²/2^(4k)", Fraction(math.prod(range(1, 2*k, 2))**2, 1 << (4*k))),
            ("(j-1)/8 + 1/2^(j+2)", Fraction(j - 1, 8) + Fraction(1, 1 << (j + 2))),
            ("(j-1)/8 + 1/4^k", Fraction(j - 1, 8) + Fraction(1, 4**k)),
            ("j²/(8j+8)", Fraction(j*j, 8*j + 8)),
            ("(4k²-1)/(16k)", Fraction(4*k*k - 1, 16*k) if k > 0 else Fraction(0)),
        ]
        for name, val in candidates:
            if cov == val:
                pr(f"    j={j}: {name} ✓✓✓ MATCH")

    # Odd-j verification
    pr(f"\n  Odd-j formula Cov = (j-1)/8:")
    all_odd_match = True
    for j, Ec, Ev, cov, _, _ in results_12:
        if j < 3 or j % 2 != 1:
            continue
        target = Fraction(j - 1, 8)
        if cov != target:
            pr(f"    j={j}: MISMATCH! Cov={cov}, target={target}")
            all_odd_match = False
    if all_odd_match:
        pr(f"    ✓ Cov(carry_j, conv_j) = (j-1)/8 for ALL odd j=3,5,...,11")

    # Induction step check
    pr(f"\n{'═' * 78}")
    pr(f"  INDUCTION STEP: Cov(j+2) - Cov(j)")
    pr(f"{'═' * 78}")

    pr(f"\n  Odd→Odd steps:")
    for j, Ec, Ev, cov, _, _ in results_12:
        if j < 3 or j % 2 != 1 or j + 2 > 12:
            continue
        cov_next = [c for jj, _, _, c, _, _ in results_12 if jj == j + 2][0]
        step = cov_next - cov
        pr(f"    Cov({j+2}) - Cov({j}) = {step} = {float(step):.8f}  "
           f"(target 1/4 = {float(Fraction(1,4)):.8f})")

    pr(f"\n  Even→Even steps:")
    for j, Ec, Ev, cov, _, _ in results_12:
        if j < 2 or j % 2 != 0 or j + 2 > 12:
            continue
        cov_next = [c for jj, _, _, c, _, _ in results_12 if jj == j + 2][0]
        step = cov_next - cov
        target_step = Fraction(1, 4)
        pr(f"    Cov({j+2}) - Cov({j}) = {step} = {float(step):.8f}  "
           f"(Δ from 1/4: {float(step - target_step):+.8f})")

    # Part 3: Sensitivity
    pr(f"\n{'═' * 78}")
    pr(f"  SENSITIVITY CONDITIONED ON c_top = 1")
    pr(f"{'═' * 78}")

    for K_sens in [8, 10]:
        pr(f"\n  K={K_sens}:")
        n_valid, n_total, sums = sensitivity_analysis(K_sens)
        pr(f"  Valid (c_top=1): {n_valid}/{n_total} "
           f"({100*n_valid/n_total:.2f}%)")

        if n_valid > 0:
            pr(f"\n  {'k':>4}  {'E[c_top-k|ULC]':>18}  {'float':>12}  {'α_k':>18}  {'α float':>12}")
            alpha_prev = Fraction(0)
            for k in range(min(K_sens + 1, len(sums))):
                Eck = Fraction(sums[k], n_valid)
                if k == 0:
                    pr(f"  {k:4d}  {str(Eck):>18}  {float(Eck):12.8f}  {'—':>18}  {'—':>12}")
                else:
                    alpha_k = Eck - alpha_prev
                    dev = float(alpha_k) - 0.25
                    pr(f"  {k:4d}  {str(Eck):>18}  {float(Eck):12.8f}  "
                       f"{str(alpha_k):>18}  {float(alpha_k):12.8f}  (Δ: {dev:+.6f})")
                alpha_prev = Eck

            c1_val = Fraction(sums[1], n_valid) - Fraction(1)
            pr(f"\n  c₁(K={K_sens}) = E[c_top-1|ULC] - 1 = {c1_val} = {float(c1_val):.10f}")
            pr(f"  target ln(2)/4 = {math.log(2)/4:.10f}")
            pr(f"  error = {float(c1_val) - math.log(2)/4:+.6e}")

    # Part 4: Carry-carry correlation at boundary (Diaconis-Fulman check)
    pr(f"\n{'═' * 78}")
    pr(f"  CARRY-CARRY AUTOCORRELATION (bulk vs boundary)")
    pr(f"{'═' * 78}")

    results_12 = all_results[12]
    K = 12
    n_bits = 2 * K
    n_configs = 1 << n_bits

    sum_cc_lag1 = [0] * (K + 1)
    sum_c = [0] * (K + 1)
    sum_csq = [0] * (K + 1)

    for cfg in range(n_configs):
        g = [0] * (K + 1)
        h = [0] * (K + 1)
        g[0] = 1
        h[0] = 1
        for i in range(K):
            g[i + 1] = (cfg >> (2 * K - 1 - i)) & 1
            h[i + 1] = (cfg >> (K - 1 - i)) & 1

        carries = [0] * (K + 2)
        carry = 0
        for j in range(K + 1):
            conv_j = 0
            for i in range(j + 1):
                if i <= K and (j - i) <= K:
                    conv_j += g[i] * h[j - i]
            carries[j] = carry
            carry = (conv_j + carry) >> 1
        carries[K + 1] = carry

        for j in range(K + 1):
            sum_c[j] += carries[j]
            sum_csq[j] += carries[j] * carries[j]
            if j < K + 1:
                sum_cc_lag1[j] += carries[j] * carries[j + 1]

    pr(f"\n  Corr(c_j, c_{{j+1}}) — theory: -1/2 for bulk (base 2)")
    for j in range(2, K):
        Ec_j = Fraction(sum_c[j], n_configs)
        Ec_j1 = Fraction(sum_c[j + 1], n_configs)
        Ecc = Fraction(sum_cc_lag1[j], n_configs)
        Vc_j = Fraction(sum_csq[j], n_configs) - Ec_j ** 2
        Vc_j1 = Fraction(sum_csq[j + 1], n_configs) - Ec_j1 ** 2
        cov = Ecc - Ec_j * Ec_j1
        if Vc_j > 0 and Vc_j1 > 0:
            corr = float(cov) / (float(Vc_j) * float(Vc_j1)) ** 0.5
        else:
            corr = 0
        pr(f"    j={j:2d}: Corr = {corr:+.6f}  (theory: -0.500000)")

    pr("\n" + "=" * 78)
    pr("  DONE")
    pr("=" * 78)


if __name__ == '__main__':
    main()
