#!/usr/bin/env python3
"""
F01_covariance_proof.py — Analytical proof of boundary covariance structure

THEOREM (to prove):
  For the binary multiplication carry chain starting from carry_0 = 0:
  1. carry_j = floor((conv_{j-1} + carry_{j-1}) / 2)
  2. E[carry_j] = (j-1)/4 exactly (unconditional bulk)
  3. Cov(carry_j, conv_j) = (j-1)/8 for odd j >= 3 (observed in F01)

APPROACH:
  Use symbolic/exact computation with all 2^{2K} bit configurations.
  For each position j, express carry_j as an EXACT polynomial in bits g_1,...,h_1,...
  Then compute E[carry_j], E[conv_j], E[carry_j · conv_j] exactly.

PART 1: Symbolic carry computation (verify diagonal product structure)
PART 2: Exact covariance for j=2..20
PART 3: Sensitivity analysis conditioned on c_top = 1 (the ULC constraint)
"""

import sys
import math
import numpy as np
from fractions import Fraction

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def exact_covariance_enumeration(K_max):
    """Compute exact E[carry_j], E[conv_j], Cov(carry_j, conv_j) for j=0..K_max
    by enumerating all 2^{2K_max} bit configurations.
    Returns Fraction-valued results for exact arithmetic."""

    n_bits = 2 * K_max
    n_configs = 1 << n_bits

    E_carry = [Fraction(0)] * (K_max + 1)
    E_conv = [Fraction(0)] * (K_max + 1)
    E_carry_conv = [Fraction(0)] * (K_max + 1)
    E_carry_sq = [Fraction(0)] * (K_max + 1)
    E_conv_sq = [Fraction(0)] * (K_max + 1)

    for cfg in range(n_configs):
        g = [0] * K_max
        h = [0] * K_max
        for i in range(K_max):
            g[i] = (cfg >> (2 * K_max - 1 - i)) & 1
            h[i] = (cfg >> (K_max - 1 - i)) & 1

        conv = [0] * (K_max + 1)
        for j in range(K_max + 1):
            s = 0
            for i in range(j + 1):
                gi = 1 if i == 0 else g[i - 1]
                hji = 1 if (j - i) == 0 else h[j - i - 1]
                if i <= K_max and (j - i) <= K_max:
                    s += gi * hji
            conv[j] = s

        carry = [0] * (K_max + 1)
        for j in range(K_max):
            carry[j + 1] = (conv[j] + carry[j]) // 2

        for j in range(K_max + 1):
            E_carry[j] += Fraction(carry[j])
            E_conv[j] += Fraction(conv[j])
            E_carry_conv[j] += Fraction(carry[j] * conv[j])
            E_carry_sq[j] += Fraction(carry[j] * carry[j])
            E_conv_sq[j] += Fraction(conv[j] * conv[j])

    N = Fraction(n_configs)
    for j in range(K_max + 1):
        E_carry[j] /= N
        E_conv[j] /= N
        E_carry_conv[j] /= N
        E_carry_sq[j] /= N
        E_conv_sq[j] /= N

    Cov = [Fraction(0)] * (K_max + 1)
    Var_c = [Fraction(0)] * (K_max + 1)
    Var_v = [Fraction(0)] * (K_max + 1)
    for j in range(K_max + 1):
        Cov[j] = E_carry_conv[j] - E_carry[j] * E_conv[j]
        Var_c[j] = E_carry_sq[j] - E_carry[j] ** 2
        Var_v[j] = E_conv_sq[j] - E_conv[j] ** 2

    return E_carry, E_conv, Cov, Var_c, Var_v


def symbolic_carry_proof():
    """Prove the carry structure analytically for positions j=1,2,3,4."""
    pr("═" * 78)
    pr("  SYMBOLIC PROOF: carry structure at boundary positions")
    pr("═" * 78)

    pr("""
  Setup: g_0 = h_0 = 1 (MSBs). g_i, h_i ~ Bernoulli(1/2) for i >= 1.
  conv_j = Σ_{i+k=j} g_i · h_k   (with g_0=h_0=1)
  carry_0 = 0
  carry_{j+1} = floor((conv_j + carry_j) / 2)

  POSITION j=0:
    conv_0 = g_0 · h_0 = 1 (constant)
    carry_0 = 0
    carry_1 = floor((1 + 0) / 2) = floor(1/2) = 0

  → carry_1 = 0  (always)  ■

  POSITION j=1:
    conv_1 = g_0·h_1 + g_1·h_0 = h_1 + g_1
    carry_1 = 0
    carry_2 = floor((h_1 + g_1 + 0) / 2) = floor((g_1 + h_1) / 2)
    Since g_1, h_1 ∈ {0,1}: g_1 + h_1 ∈ {0,1,2}
      g_1+h_1=0: carry_2 = 0
      g_1+h_1=1: carry_2 = 0
      g_1+h_1=2: carry_2 = 1
    So carry_2 = g_1 · h_1  (product of independent bits)

  → carry_2 = g_1 · h_1    ■

  POSITION j=2:
    conv_2 = g_0·h_2 + g_1·h_1 + g_2·h_0 = h_2 + g_1·h_1 + g_2
    carry_2 = g_1·h_1
    conv_2 + carry_2 = h_2 + 2·g_1·h_1 + g_2
    carry_3 = floor((h_2 + 2·g_1·h_1 + g_2) / 2)

    Since g_1·h_1 is integer ∈ {0,1}:
      carry_3 = g_1·h_1 + floor((h_2 + g_2) / 2)
              = g_1·h_1 + g_2·h_2

  → carry_3 = g_1·h_1 + g_2·h_2  (sum of diagonal products)  ■

  POSITION j=3:
    conv_3 = h_3 + g_1·h_2 + g_2·h_1 + g_3
    carry_3 = g_1·h_1 + g_2·h_2 ∈ {0, 1, 2}

    conv_3 + carry_3 = h_3 + g_1·h_2 + g_2·h_1 + g_3 + g_1·h_1 + g_2·h_2
    
    Case carry_3 = 0 (g_1·h_1=0 AND g_2·h_2=0, prob 9/16):
      carry_4 = floor((h_3 + g_1·h_2 + g_2·h_1 + g_3) / 2)
      conv_3 has 4 independent Bernoulli terms, each with E = 1/4 or 1/2.
      Actually: h_3~Ber(1/2), g_1·h_2~Ber(1/4), g_2·h_1~Ber(1/4), g_3~Ber(1/2)
      Sum = h_3 + g_3 + g_1·h_2 + g_2·h_1, range [0,4]
      
    Case carry_3 = 1 (g_1·h_1 XOR g_2·h_2 = 1, prob 6/16):
      carry_4 = floor((conv_3 + 1) / 2)
      
    Case carry_3 = 2 (g_1·h_1=1 AND g_2·h_2=1, prob 1/16):
      carry_4 = floor((conv_3 + 2) / 2) = 1 + floor(conv_3 / 2)

    The structure becomes more complex at j=4 because carry_3 can be 2.
    The simple "diagonal product" structure breaks for carry_4.
    """)


def main():
    pr("=" * 78)
    pr("  BOUNDARY COVARIANCE PROOF AND SENSITIVITY ANALYSIS")
    pr("=" * 78)

    symbolic_carry_proof()

    # ═══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART 1: EXACT COVARIANCE ENUMERATION")
    pr("═" * 78)

    for K in [8, 10, 12]:
        pr(f"\n  K = {K} ({1 << (2*K)} configurations):")
        if K > 12:
            pr("    Skipping (too large)")
            continue

        E_c, E_v, Cov, Var_c, Var_v = exact_covariance_enumeration(K)

        pr(f"  {'j':>3s}  {'E[carry]':>12s}  {'E[conv]':>12s}  "
           f"{'Cov(c,v)':>16s}  {'Cov float':>12s}  {'(j-1)/8':>10s}  {'match?':>7s}")

        for j in range(K + 1):
            cov_f = float(Cov[j])
            target = Fraction(j - 1, 8) if j >= 1 else Fraction(0)
            target_f = float(target)
            match = "✓" if Cov[j] == target and j % 2 == 1 and j >= 3 else ""
            if j >= 2 and Cov[j] != Fraction(0):
                pr(f"  {j:3d}  {str(E_c[j]):>12s}  {str(E_v[j]):>12s}  "
                   f"{str(Cov[j]):>16s}  {cov_f:12.8f}  {target_f:10.6f}  {match:>7s}")

        pr(f"\n  Verification: E[carry_j] = (j-1)/4 for j>=1 ?  (Theorem 3)")
        all_match = True
        for j in range(K + 1):
            expected = Fraction(j - 1, 4) if j >= 1 else Fraction(0)
            if E_c[j] != expected:
                pr(f"    j={j}: E[carry] = {E_c[j]} ≠ {expected} !")
                all_match = False
        if all_match:
            pr(f"    ✓ E[carry_j] = (j-1)/4 exactly for all j=1..{K}, carry_0=0")

        pr(f"\n  Covariance pattern analysis:")
        pr(f"  Odd j (j=3,5,...): Cov = (j-1)/8 ?")
        for j in range(3, K + 1, 2):
            target = Fraction(j - 1, 8)
            match = "✓" if Cov[j] == target else f"✗ ({Cov[j]})"
            pr(f"    j={j}: Cov = {Cov[j]} = {float(Cov[j]):.8f}, (j-1)/8 = {target} {match}")

        pr(f"\n  Even j (j=2,4,...): pattern search")
        even_covs = []
        for j in range(2, K + 1, 2):
            pr(f"    j={j}: Cov = {Cov[j]} = {float(Cov[j]):.10f}")
            even_covs.append((j, Cov[j]))

        if len(even_covs) >= 2:
            pr(f"\n  Testing even-j formulas:")
            for j, c in even_covs:
                tests = [
                    (f"(j-1)/8", Fraction(j-1, 8)),
                    (f"(j²-1)/32", Fraction(j*j-1, 32)),
                    (f"(2j-1)(2j-3)/128", Fraction((2*j-1)*(2*j-3), 128)),
                    (f"j(j-1)/16", Fraction(j*(j-1), 16)),
                    (f"(3j²-4j+1)/32", Fraction(3*j*j-4*j+1, 32)),
                ]
                for name, val in tests:
                    if c == val:
                        pr(f"    j={j}: Cov = {name} ✓")
                        break

    # ═══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART 2: CLOSED FORM SEARCH FOR EVEN-j COVARIANCES")
    pr("═" * 78)

    pr(f"\n  Using K=12 data for the even-j covariance sequence:")
    E_c12, E_v12, Cov12, _, _ = exact_covariance_enumeration(12)

    even_data = []
    for j in range(2, 13, 2):
        pr(f"    j={j}: Cov = {Cov12[j]}")
        even_data.append((j, Cov12[j]))

    pr(f"\n  Numerators and denominators:")
    for j, c in even_data:
        pr(f"    j={j}: {c.numerator} / {c.denominator}")

    pr(f"\n  Successive ratios:")
    for i in range(1, len(even_data)):
        j1, c1 = even_data[i-1]
        j2, c2 = even_data[i]
        ratio = c2 / c1
        pr(f"    Cov({j2})/Cov({j1}) = {ratio} = {float(ratio):.8f}")

    pr(f"\n  Second differences of Cov(j)/((j-1)/8):")
    for j in range(2, 13, 2):
        ratio = Cov12[j] / Fraction(j-1, 8)
        pr(f"    j={j}: Cov / ((j-1)/8) = {ratio} = {float(ratio):.8f}")

    # ═══════════════════════════════════════════════════════════════
    pr("\n" + "═" * 78)
    pr("  PART 3: SENSITIVITY CONDITIONED ON c_top = 1")
    pr("═" * 78)

    pr(f"\n  For the trace anomaly, we need E[c_{{top-1}} | c_top = 1].")
    pr(f"  This requires computing the full carry chain for REAL semiprimes")
    pr(f"  and conditioning on the top carry being 1 (ULC).")
    pr(f"\n  Using K=10 with MSB conditioning:")

    K = 10
    n_bits = 2 * K
    n_configs = 1 << n_bits
    D_total = 2 * (K + 1)

    alpha_cond = [Fraction(0)] * (K + 2)
    n_valid = 0

    for cfg in range(n_configs):
        g_bits = [0] * K
        h_bits = [0] * K
        for i in range(K):
            g_bits[i] = (cfg >> (2 * K - 1 - i)) & 1
            h_bits[i] = (cfg >> (K - 1 - i)) & 1

        g_full = [1] + g_bits
        h_full = [1] + h_bits
        d = K + 1

        conv = [0] * D_total
        for i in range(d):
            for j_idx in range(d):
                pos = i + j_idx
                if pos < D_total:
                    conv[pos] += g_full[d - 1 - i] * h_full[d - 1 - j_idx]

        carries = [0] * (D_total + 1)
        for pos in range(D_total):
            carries[pos + 1] = (conv[pos] + carries[pos]) // 2

        M_top = D_total
        while M_top > 0 and carries[M_top] == 0:
            M_top -= 1

        if carries[M_top] != 1:
            continue

        n_valid += 1
        for k in range(min(K + 2, M_top + 1)):
            alpha_cond[k] += Fraction(carries[M_top - k])

    pr(f"  Valid samples (c_top=1): {n_valid} / {n_configs}")

    if n_valid > 0:
        N_valid = Fraction(n_valid)
        for k in range(K + 2):
            alpha_cond[k] /= N_valid

        pr(f"\n  E[c_{{top-k}} | c_top = 1]:")
        for k in range(min(K + 1, len(alpha_cond))):
            pr(f"    k={k:2d}: {alpha_cond[k]} = {float(alpha_cond[k]):.10f}")

        pr(f"\n  α_k = E[c_{{top-k}}] - E[c_{{top-k+1}}] conditioned on c_top=1:")
        for k in range(1, min(K + 1, len(alpha_cond))):
            ak = alpha_cond[k] - alpha_cond[k - 1]
            dev = float(ak) - 0.25
            pr(f"    α_{k:2d} = {ak} = {float(ak):.10f}  (Δ from 1/4: {dev:+.8f})")

        c1_cond = alpha_cond[1] - alpha_cond[0]
        pr(f"\n  c₁(K={K}) = α₁ = {c1_cond} = {float(c1_cond):.10f}")
        pr(f"  target: ln(2)/4 = {math.log(2)/4:.10f}")

    pr("\n" + "=" * 78)
    pr("  END")
    pr("=" * 78)


if __name__ == '__main__':
    main()
