#!/usr/bin/env python3
"""
F03_fast_covariance.py — Fast covariance computation and closed-form search

From F01 we know:
  - E[carry_j] = (j-1)/4 exactly (unconditional)
  - Cov(carry_j, conv_j) = (j-1)/8 for odd j >= 3
  - For even j = 2n: ratio = 8·Cov(2n)/(2n-1) = (D_n + 1)/D_n
    where D_1=2, D_2=12, D_3=80, D_4=224, D_5=1536

This script:
  PART A: Fast numpy enumeration for K=12,14 to extend the D_n sequence
  PART B: Pattern search on D_n 
  PART C: Conditioned sensitivity (c_top = 1) via fast enumeration
  PART D: Total c₁ from weighted covariance sum
"""

import sys
import math
import numpy as np
from fractions import Fraction

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def fast_covariance(K):
    """Compute Cov(carry_j, conv_j) using numpy with exact integer arithmetic."""
    n_bits = 2 * K
    n_configs = 1 << n_bits

    sum_carry = np.zeros(K + 1, dtype=np.float64)
    sum_conv = np.zeros(K + 1, dtype=np.float64)
    sum_cc = np.zeros(K + 1, dtype=np.float64)

    batch_size = min(n_configs, 1 << 20)
    n_batches = (n_configs + batch_size - 1) // batch_size

    for b_idx in range(n_batches):
        start = b_idx * batch_size
        end = min(start + batch_size, n_configs)
        cfgs = np.arange(start, end, dtype=np.int64)

        g = np.zeros((end - start, K), dtype=np.int8)
        h = np.zeros((end - start, K), dtype=np.int8)
        for i in range(K):
            g[:, i] = (cfgs >> (2 * K - 1 - i)) & 1
            h[:, i] = (cfgs >> (K - 1 - i)) & 1

        g_full = np.ones((end - start, K + 1), dtype=np.int8)
        h_full = np.ones((end - start, K + 1), dtype=np.int8)
        g_full[:, 1:] = g
        h_full[:, 1:] = h

        carries = np.zeros((end - start, K + 1), dtype=np.int32)
        convs = np.zeros((end - start, K + 1), dtype=np.int32)

        for j in range(K + 1):
            s = np.zeros(end - start, dtype=np.int32)
            for i in range(j + 1):
                if i <= K and (j - i) <= K:
                    s += g_full[:, i].astype(np.int32) * h_full[:, j - i].astype(np.int32)
            convs[:, j] = s
            if j < K:
                carries[:, j + 1] = (s + carries[:, j]) // 2

        for j in range(K + 1):
            sum_carry[j] += np.sum(carries[:, j].astype(np.float64))
            sum_conv[j] += np.sum(convs[:, j].astype(np.float64))
            sum_cc[j] += np.sum(carries[:, j].astype(np.float64) * convs[:, j].astype(np.float64))

    N = float(n_configs)
    E_carry = sum_carry / N
    E_conv = sum_conv / N
    E_cc = sum_cc / N
    Cov = E_cc - E_carry * E_conv

    return E_carry, E_conv, Cov


def exact_fraction_covariance(K):
    """Exact Fraction arithmetic for small K."""
    n_bits = 2 * K
    n_configs = 1 << n_bits

    E_carry = [Fraction(0)] * (K + 1)
    E_conv = [Fraction(0)] * (K + 1)
    E_cc = [Fraction(0)] * (K + 1)

    for cfg in range(n_configs):
        g = [(cfg >> (2 * K - 1 - i)) & 1 for i in range(K)]
        h = [(cfg >> (K - 1 - i)) & 1 for i in range(K)]

        g_full = [1] + g
        h_full = [1] + h

        carry = [0] * (K + 1)
        conv_arr = [0] * (K + 1)

        for j in range(K + 1):
            s = 0
            for i in range(j + 1):
                if i <= K and (j - i) <= K:
                    s += g_full[i] * h_full[j - i]
            conv_arr[j] = s
            if j < K:
                carry[j + 1] = (s + carry[j]) // 2

        for j in range(K + 1):
            E_carry[j] += Fraction(carry[j])
            E_conv[j] += Fraction(conv_arr[j])
            E_cc[j] += Fraction(carry[j] * conv_arr[j])

    N = Fraction(n_configs)
    Cov = [Fraction(0)] * (K + 1)
    for j in range(K + 1):
        E_carry[j] /= N
        E_conv[j] /= N
        E_cc[j] /= N
        Cov[j] = E_cc[j] - E_carry[j] * E_conv[j]

    return E_carry, E_conv, Cov


def sensitivity_conditioned(K):
    """Compute E[carry_j | c_top = 1] and related quantities.
    Uses the FULL carry chain for semiprimes (g_0=h_0=1, D=2(K+1) digits)."""

    n_bits = 2 * K
    n_configs = 1 << n_bits
    d = K + 1
    D = 2 * d

    alpha_sum = np.zeros(K + 2, dtype=np.float64)
    cov_sum = np.zeros(K + 2, dtype=np.float64)
    n_valid = 0

    batch_size = min(n_configs, 1 << 18)
    n_batches = (n_configs + batch_size - 1) // batch_size

    for b_idx in range(n_batches):
        start = b_idx * batch_size
        end = min(start + batch_size, n_configs)
        cfgs = np.arange(start, end, dtype=np.int64)
        n_cur = end - start

        g_bits = np.zeros((n_cur, K), dtype=np.int8)
        h_bits = np.zeros((n_cur, K), dtype=np.int8)
        for i in range(K):
            g_bits[:, i] = (cfgs >> (2 * K - 1 - i)) & 1
            h_bits[:, i] = (cfgs >> (K - 1 - i)) & 1

        g_full = np.zeros((n_cur, d), dtype=np.int8)
        h_full = np.zeros((n_cur, d), dtype=np.int8)
        g_full[:, 0] = 1
        h_full[:, 0] = 1
        g_full[:, 1:] = g_bits
        h_full[:, 1:] = h_bits

        g_rev = g_full[:, ::-1]
        h_rev = h_full[:, ::-1]

        conv = np.zeros((n_cur, D), dtype=np.int32)
        for pos in range(D):
            for i in range(d):
                j_idx = pos - i
                if 0 <= j_idx < d:
                    conv[:, pos] += g_rev[:, i].astype(np.int32) * h_rev[:, j_idx].astype(np.int32)

        carries = np.zeros((n_cur, D + 1), dtype=np.int32)
        for pos in range(D):
            carries[:, pos + 1] = (conv[:, pos] + carries[:, pos]) // 2

        M_top = np.zeros(n_cur, dtype=np.int32)
        for pos in range(D, 0, -1):
            mask = (M_top == 0) & (carries[:, pos] > 0)
            M_top[mask] = pos

        valid = (carries[np.arange(n_cur), M_top] == 1)
        valid_idx = np.where(valid)[0]
        n_v = len(valid_idx)
        n_valid += n_v

        if n_v == 0:
            continue

        for k in range(min(K + 2, D)):
            positions = M_top[valid_idx] - k
            vals = np.where(positions >= 0,
                          carries[valid_idx, np.clip(positions, 0, D)], 0)
            alpha_sum[k] += np.sum(vals.astype(np.float64))

    if n_valid == 0:
        return None, None

    E_top_k = alpha_sum / n_valid
    alpha = np.diff(E_top_k)

    return E_top_k, alpha, n_valid


def main():
    pr("=" * 78)
    pr("  FAST COVARIANCE ANALYSIS AND CLOSED-FORM SEARCH")
    pr("=" * 78)

    # ══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART A: EXACT COVARIANCE (Fraction) FOR K=10")
    pr("═" * 78)

    K = 10
    pr(f"\n  Computing exact Fraction covariance for K={K}...")
    E_c, E_v, Cov = exact_fraction_covariance(K)

    pr(f"\n  E[carry_j] vs (j-1)/4:")
    for j in range(K + 1):
        expected = Fraction(j - 1, 4) if j >= 1 else Fraction(0)
        match = "✓" if E_c[j] == expected else f"✗ ({E_c[j]})"
        if j <= 2 or j >= K - 1:
            pr(f"    j={j}: E[c] = {E_c[j]} vs {expected} {match}")

    pr(f"\n  Exact covariances:")
    pr(f"  {'j':>3s}  {'Cov (fraction)':>20s}  {'Cov (float)':>14s}  {'(j-1)/8':>10s}  "
       f"{'δ = Cov-(j-1)/8':>18s}")

    Dn_list = []
    for j in range(2, K + 1):
        target = Fraction(j - 1, 8)
        delta = Cov[j] - target
        pr(f"  {j:3d}  {str(Cov[j]):>20s}  {float(Cov[j]):14.10f}  "
           f"{float(target):10.6f}  {str(delta):>18s}")
        if j % 2 == 0 and delta != 0:
            ratio = Cov[j] / target
            Dn = Fraction(1) / (ratio - 1)
            Dn_list.append((j, Dn))

    pr(f"\n  Even-j: ratio = Cov(j) / ((j-1)/8) and D_n = 1/(ratio - 1):")
    for j, Dn in Dn_list:
        n = j // 2
        pr(f"    n={n} (j={j}): D_n = {Dn} = {float(Dn):.6f}")

    # ══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART A2: NUMPY COVARIANCE FOR K=12,14")
    pr("═" * 78)

    for K in [12, 14]:
        pr(f"\n  K={K} ({1 << (2*K)} configs, numpy)...")
        E_c, E_v, Cov_np = fast_covariance(K)

        pr(f"\n  E[carry_j] = (j-1)/4:")
        all_good = True
        for j in range(1, K + 1):
            expected = (j - 1) / 4.0
            err = abs(E_c[j] - expected)
            if err > 1e-10:
                pr(f"    j={j}: MISMATCH E[c]={E_c[j]:.10f} vs {expected:.10f}")
                all_good = False
        if all_good:
            pr(f"    ✓ All match to <1e-10")

        pr(f"\n  Covariances (even j only, since odd = (j-1)/8 exact):")
        pr(f"  {'j':>3s}  {'Cov':>16s}  {'(j-1)/8':>12s}  {'δ':>14s}  {'ratio':>14s}")

        for j in range(2, K + 1):
            target = (j - 1) / 8.0
            delta = Cov_np[j] - target

            if j % 2 == 1:
                if abs(delta) < 1e-8:
                    continue
                else:
                    pr(f"  {j:3d}  {Cov_np[j]:16.12f}  {target:12.8f}  "
                       f"{delta:14.2e}  ODD MISMATCH!")
            else:
                ratio = Cov_np[j] / target if target > 0 else float('inf')
                pr(f"  {j:3d}  {Cov_np[j]:16.12f}  {target:12.8f}  "
                   f"{delta:14.2e}  {ratio:14.12f}")

    # ══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART B: D_n PATTERN SEARCH")
    pr("═" * 78)

    D_exact = {
        1: Fraction(2),
        2: Fraction(12),
        3: Fraction(80),
        4: Fraction(224),
        5: Fraction(1536),
    }

    pr(f"\n  Known exact D_n values:")
    for n, d in sorted(D_exact.items()):
        pr(f"    D_{n} = {d} = {d.numerator}")

    pr(f"\n  Factorizations:")
    for n, d in sorted(D_exact.items()):
        val = d.numerator
        factors = []
        v = val
        for p in [2, 3, 5, 7, 11, 13]:
            while v % p == 0:
                factors.append(p)
                v //= p
        if v > 1:
            factors.append(v)
        pr(f"    D_{n} = {val} = {'·'.join(map(str, factors))}")

    pr(f"\n  Testing formulas:")

    pr(f"  Test: D_n = 2^(2n-1) · (2n-1)!! / (2n-1) ?")
    for n in range(1, 6):
        dbl_fact = 1
        for k in range(1, 2 * n, 2):
            dbl_fact *= k
        test_val = Fraction(2 ** (2 * n - 1) * dbl_fact, 2 * n - 1)
        pr(f"    n={n}: 2^{2*n-1}·{dbl_fact}/{2*n-1} = {test_val} (D_n = {D_exact[n]}) "
           f"{'✓' if test_val == D_exact[n] else '✗'}")

    pr(f"\n  Test: D_n = 4^n · Catalan(n-1)?")
    for n in range(1, 6):
        cat = math.comb(2 * (n - 1), n - 1) // n
        test_val = 4 ** n * cat
        pr(f"    n={n}: 4^{n}·C_{n-1} = {test_val} (D_n = {D_exact[n]}) "
           f"{'✓' if test_val == D_exact[n].numerator else '✗'}")

    pr(f"\n  Test: D_n = binom(2n,n) · 2^? ?")
    for n in range(1, 6):
        bc = math.comb(2 * n, n)
        ratio = D_exact[n] / Fraction(bc)
        pr(f"    n={n}: binom({2*n},{n}) = {bc}, D_n/binom = {ratio} = {float(ratio):.6f}")

    pr(f"\n  Test: D_n = (2n)! / (n! · 2^n) · something ?")
    for n in range(1, 6):
        base = math.factorial(2 * n) // (math.factorial(n) * 2 ** n)
        ratio = D_exact[n] / Fraction(base)
        pr(f"    n={n}: (2n)!/(n!·2^n) = {base}, D_n/base = {ratio} = {float(ratio):.6f}")

    pr(f"\n  Test: D_n as a product:")
    for n in range(1, 6):
        val = D_exact[n].numerator
        test_prod = 1
        for k in range(1, n + 1):
            test_prod *= (4 * k - 2)
        pr(f"    n={n}: ∏(4k-2) = {test_prod}, D_n = {val}, D_n/prod = {Fraction(val, test_prod)}")

    pr(f"\n  Recursive check: D_{'{n+1}'}/D_n:")
    for n in range(1, 5):
        ratio = D_exact[n + 1] / D_exact[n]
        pr(f"    D_{n+1}/D_n = {ratio} = {float(ratio):.8f}")

    pr(f"\n  Test: D_n · (2n-1) = ?")
    for n in range(1, 6):
        val = D_exact[n] * (2 * n - 1)
        pr(f"    D_{n}·{2*n-1} = {val}")

    pr(f"\n  Test WALLIS: D_n = 2·∏_{{k=1}}^{{n-1}} [4k(k+1)/(2k+1)] ?")
    for n in range(1, 6):
        prod = Fraction(2)
        for k in range(1, n):
            prod *= Fraction(4 * k * (k + 1), 2 * k + 1)
        pr(f"    n={n}: val = {prod} = {float(prod):.6f} (D_n = {D_exact[n]}) "
           f"{'✓' if prod == D_exact[n] else '✗'}")

    pr(f"\n  Numerator search: D_n = 2 · a_n, a_n = ?")
    for n in range(1, 6):
        a = D_exact[n] // 2
        pr(f"    a_{n} = {a}")
    pr(f"    Sequence: 1, 6, 40, 112, 768")
    pr(f"    = 1, 2·3, 8·5, 16·7, 256·3")

    # Now compute with Fraction for δ_n = (2n-1)/(8·D_n)
    pr(f"\n  UNIFIED FORMULA: Cov(2n) = (2n-1)/8 + (2n-1)/(8·D_n)")
    for n in range(1, 6):
        j = 2 * n
        cov_formula = Fraction(2 * n - 1, 8) + Fraction(2 * n - 1, 8) / D_exact[n]
        pr(f"    j={j}: (2n-1)/8·(1 + 1/D_n) = {cov_formula}")

    # ══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART C: SENSITIVITY CONDITIONED ON c_top = 1")
    pr("═" * 78)

    for K in [8, 10]:
        pr(f"\n  K={K}: Computing conditioned carry profile...")
        result = sensitivity_conditioned(K)
        if result is None:
            pr(f"    No valid samples!")
            continue

        E_top_k, alpha, n_valid = result
        n_total = 1 << (2 * K)
        pr(f"  Valid (c_top=1): {n_valid}/{n_total} = {n_valid/n_total:.6f}")

        pr(f"\n  E[c_{{top-k}} | c_top = 1] (the carry profile from MSB):")
        for k in range(min(K + 1, len(E_top_k))):
            if E_top_k[k] > 0 or k <= 2:
                pr(f"    k={k:2d}: {E_top_k[k]:.10f}")

        pr(f"\n  α_k = E[c_{{top-k}}] - E[c_{{top-k+1}}] (increment per position):")
        for k in range(min(K, len(alpha))):
            dev = alpha[k] - 0.25
            pr(f"    α_{k+1:2d} = {alpha[k]:.10f}  (Δ from 1/4: {dev:+.10f})")

        c1_K = E_top_k[1] - E_top_k[0]
        pr(f"\n  c₁(K={K}) = α₁ = {c1_K:.10f}")
        pr(f"  target: ln(2)/4 = {math.log(2)/4:.10f}")

        pr(f"\n  SENSITIVITY WEIGHTS: ∂c₁/∂(Cov at position j)")
        pr(f"  The conditioned covariance structure tells us how much each")
        pr(f"  position contributes to c₁ = E[c_{{M-1}} | ULC] - 1.")
        pr(f"  Key question: does α_k → 1/4 exponentially? At what rate?")

        if K >= 6:
            for k in range(2, min(K, len(alpha))):
                if abs(alpha[k] - 0.25) > 1e-8:
                    eff_rate = abs(alpha[k] - 0.25) / abs(alpha[k - 1] - 0.25) if abs(alpha[k - 1] - 0.25) > 1e-8 else float('inf')
                    pr(f"    |α_{k+1} - 1/4| / |α_{k} - 1/4| = {eff_rate:.6f}")

        S_boundary = sum(alpha[k] - 0.25 for k in range(len(alpha)))
        pr(f"\n  S_boundary = Σ(α_k - 1/4) = {S_boundary:.10f}")
        pr(f"  Expected: c₁ - 1/4 = {math.log(2)/4 - 0.25:.10f}")
        pr(f"  c₁ = 1/4 + S_boundary = {0.25 + S_boundary:.10f}")

    # ══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART D: PERTURBATION DECOMPOSITION OF c₁")
    pr("═" * 78)

    pr(f"""
  The trace anomaly can be written as:
    c₁ = E[c_{{M-1}} | c_top=1] - 1

  Using Theorem 5 (backward recursion):
    c_{{top-1}} = 2·c_top + f_{{m-1}} - conv_{{m-1}}
    E[c_{{top-1}} | c_top=1] = 2 + E[f_{{m-1}}] - E[conv_{{m-1}} | c_top=1]

  Since f_{{m-1}} and conv_{{m-1}} are CONDITIONED on the carry chain reaching
  c_top = 1, the non-Markovian correction enters here.

  The Markov approximation ignores this conditioning:
    c₁^Markov ≈ E[carry_j/j·j/4] ... simplified model gives ~0.1353

  The NON-Markovian correction δ_NM = c₁ - c₁^Markov ≈ 0.038
  comes from Cov(carry_j, conv_j | c_top=1) ≠ Cov(carry_j, conv_j) (unconditional)

  THE KEY QUESTION: What is Cov(carry_j, conv_j | c_top=1)?
  This is DIFFERENT from the unconditional covariance we computed above.
  The conditioning on c_top = 1 modifies the distribution of all carries.
  """)

    pr(f"  Computing CONDITIONED covariance (c_top = 1) via exact enumeration...")

    K = 8
    n_bits = 2 * K
    n_configs = 1 << n_bits
    d = K + 1
    D = 2 * d

    cond_E_carry = np.zeros(D + 1)
    cond_E_conv = np.zeros(D)
    cond_E_cc = np.zeros(D)
    n_valid = 0

    for cfg in range(n_configs):
        g = [1] + [(cfg >> (2 * K - 1 - i)) & 1 for i in range(K)]
        h = [1] + [(cfg >> (K - 1 - i)) & 1 for i in range(K)]

        g_rev = g[::-1]
        h_rev = h[::-1]

        conv = [0] * D
        for pos in range(D):
            for i in range(d):
                j_idx = pos - i
                if 0 <= j_idx < d:
                    conv[pos] += g_rev[i] * h_rev[j_idx]

        carries = [0] * (D + 1)
        for pos in range(D):
            carries[pos + 1] = (conv[pos] + carries[pos]) // 2

        M_top = 0
        for pos in range(D, 0, -1):
            if carries[pos] > 0:
                M_top = pos
                break

        if carries[M_top] != 1:
            continue

        n_valid += 1

        for pos in range(D):
            j_from_top = M_top - 1 - pos
            if j_from_top >= 0:
                cond_E_carry[pos] += carries[pos]
                cond_E_conv[pos] += conv[pos]
                cond_E_cc[pos] += carries[pos] * conv[pos]

    pr(f"  Valid (c_top=1): {n_valid}/{n_configs} = {n_valid/n_configs:.6f}")
    pr(f"\n  Conditioned Cov(carry_pos, conv_pos | c_top=1):")
    pr(f"  Positions from LSB (pos=0 is least significant):")

    cond_E_carry /= n_valid
    cond_E_conv /= n_valid
    cond_E_cc /= n_valid
    cond_Cov = cond_E_cc - cond_E_carry * cond_E_conv

    pr(f"  {'pos':>4s}  {'E[carry]':>12s}  {'E[conv]':>12s}  {'Cov(c,v)':>14s}  "
       f"{'uncond (j-1)/8':>16s}")
    for pos in range(1, D - 1):
        if cond_E_carry[pos] > 0.01:
            uncond = (pos - 1) / 8.0 if pos > 0 else 0.0
            pr(f"  {pos:4d}  {cond_E_carry[pos]:12.8f}  {cond_E_conv[pos]:12.8f}  "
               f"{cond_Cov[pos]:14.10f}  {uncond:16.8f}")

    pr("\n" + "=" * 78)
    pr("  END")
    pr("=" * 78)


if __name__ == '__main__':
    main()
