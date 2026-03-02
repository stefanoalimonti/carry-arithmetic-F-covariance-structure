#!/usr/bin/env python3
"""
F05_pairing_analysis.py — Decompose Cov(carry_j, conv_j) into individual
bit-product contributions to verify the pairing cancellation mechanism.

For odd j, the covariance decomposes as:
  Cov(carry_j, conv_j) = Σ_{i=1}^{j-1} Cov(carry_j, g_i h_{j-i})

Under the symmetry g_k ↔ h_k:
  Cov(carry_j, g_i h_{j-i}) = Cov(carry_j, g_{j-i} h_i)

For odd j, ALL terms pair (no diagonal), and we verify the sum.
For even j, the diagonal term i=j/2 is unpaired → the ε_j correction.
"""

import sys
from fractions import Fraction

def pr(*a, **kw):
    print(*a, **kw)
    sys.stdout.flush()


def decompose_covariance(K):
    """Compute Cov(carry_j, g_i h_{j-i}) for each i, all j."""
    n = 1 << (2 * K)

    E_carry = [[Fraction(0)] for _ in range(K + 1)]
    E_prod = {}
    E_cross = {}

    for cfg in range(n):
        g = [1] + [(cfg >> (2*K - 1 - i)) & 1 for i in range(K)]
        h = [1] + [(cfg >> (K - 1 - i)) & 1 for i in range(K)]

        carry = 0
        for j in range(K + 1):
            conv_j = sum(g[i] * h[j - i] for i in range(j + 1)
                        if i <= K and j - i <= K)

            for i in range(1, j):
                if i <= K and j - i <= K:
                    p = g[i] * h[j - i]
                    key = (j, i)
                    if key not in E_prod:
                        E_prod[key] = Fraction(0)
                        E_cross[key] = Fraction(0)
                    E_prod[key] += Fraction(p)
                    E_cross[key] += Fraction(carry * p)

            E_carry[j][0] += Fraction(carry)
            carry = (conv_j + carry) >> 1

    N = Fraction(n)
    Ec = [E_carry[j][0] / N for j in range(K + 1)]

    result = {}
    for (j, i), s in E_cross.items():
        ep = E_prod[(j, i)] / N
        ec = Ec[j]
        cov = s / N - ec * ep
        result[(j, i)] = cov

    return result, Ec


def main():
    K = 10
    pr("=" * 78)
    pr(f"  COVARIANCE DECOMPOSITION (K={K})")
    pr("=" * 78)

    decomp, Ec = decompose_covariance(K)

    for j in range(3, K + 1):
        pr(f"\n  j = {j} ({'odd' if j % 2 == 1 else 'even'}):")
        pr(f"  {'i':>4}  {'j-i':>4}  {'Cov(c_j, g_i h_{j-i})':>24}  "
           f"{'paired (j-i,i)':>24}  {'sum of pair':>24}")

        total_cov = Fraction(0)
        for i in range(1, j):
            key = (j, i)
            if key not in decomp:
                continue
            c = decomp[key]
            total_cov += c

            mirror_key = (j, j - i)
            c_mirror = decomp.get(mirror_key, Fraction(0))

            if i <= j - i:
                if i == j - i:
                    pr(f"  {i:4d}  {j-i:4d}  {str(c):>24}  "
                       f"{'(diagonal)':>24}  {str(c):>24}")
                else:
                    pair_sum = c + c_mirror
                    pr(f"  {i:4d}  {j-i:4d}  {str(c):>24}  "
                       f"{str(c_mirror):>24}  {str(pair_sum):>24}")

        target = Fraction(j - 1, 8)
        delta = total_cov - target
        pr(f"  Total Cov = {total_cov} = {float(total_cov):.10f}")
        pr(f"  Target (j-1)/8 = {target}, δ = {delta}")

    # Now look at the INDUCTION STEP decomposition
    pr(f"\n{'═' * 78}")
    pr(f"  INDUCTION STEP DECOMPOSITION (odd j → odd j+2)")
    pr(f"{'═' * 78}")

    for j in [3, 5, 7]:
        j2 = j + 2
        if j2 > K:
            break

        pr(f"\n  Step j={j} → j+2={j2}:")
        pr(f"  ΔCov = Cov(carry_{j2}, conv_{j2}) - Cov(carry_{j}, conv_{j})")

        cov_j = sum(decomp.get((j, i), Fraction(0)) for i in range(1, j))
        cov_j2 = sum(decomp.get((j2, i), Fraction(0)) for i in range(1, j2))
        delta = cov_j2 - cov_j
        pr(f"  Cov({j}) = {cov_j}, Cov({j2}) = {cov_j2}, Δ = {delta}")

        # Decompose the delta into contributions
        pr(f"\n  Decomposition of ΔCov:")
        pr(f"  (A) New terms (i=j and i=j+1 in conv_{j2}):")
        for i in [j, j + 1]:
            key = (j2, i)
            c = decomp.get(key, Fraction(0))
            mirror = decomp.get((j2, j2 - i), Fraction(0))
            pr(f"      Cov(c_{j2}, g_{i} h_{j2-i}) = {c}")
            pr(f"      Mirror: Cov(c_{j2}, g_{j2-i} h_{i}) = {mirror}")

        pr(f"\n  (B) Evolved terms (i=1..{j-1}):")
        pr(f"      Cov(c_{j2}, g_i h_{{{j2}-i}}) vs Cov(c_{j}, g_i h_{{{j}-i}}):")
        evolved_sum = Fraction(0)
        for i in range(1, j):
            c_old = decomp.get((j, i), Fraction(0))
            c_new = decomp.get((j2, i), Fraction(0))
            diff = c_new - c_old
            evolved_sum += diff
            if i <= (j - 1) // 2:
                mirror_old = decomp.get((j, j - i), Fraction(0))
                mirror_new = decomp.get((j2, j2 - i), Fraction(0))
                pair_diff = (c_new + mirror_new) - (c_old + mirror_old)
                pr(f"      i={i}: pair Δ = {pair_diff}")
        pr(f"      Total evolved: {evolved_sum}")

    pr(f"\n{'═' * 78}")
    pr(f"  PAIR SUMS FOR ODD j (checking constancy)")
    pr(f"{'═' * 78}")

    for j in [3, 5, 7, 9]:
        if j > K:
            break
        pr(f"\n  j = {j}:")
        for i in range(1, (j + 1) // 2):
            c1 = decomp.get((j, i), Fraction(0))
            c2 = decomp.get((j, j - i), Fraction(0))
            pair = c1 + c2
            pr(f"    pair ({i},{j-i}): {str(c1):>12} + {str(c2):>12} = {str(pair):>12}")

    pr(f"\n{'═' * 78}")
    pr(f"  SECOND MOMENT ANALYSIS: E[carry_j²] and Var[carry_j]")
    pr(f"{'═' * 78}")

    K2 = 10
    n2 = 1 << (2 * K2)
    S_c2 = [0] * (K2 + 1)
    S_c = [0] * (K2 + 1)

    for cfg in range(n2):
        g = [1] + [(cfg >> (2*K2 - 1 - i)) & 1 for i in range(K2)]
        h = [1] + [(cfg >> (K2 - 1 - i)) & 1 for i in range(K2)]
        carry = 0
        for j in range(K2 + 1):
            S_c[j] += carry
            S_c2[j] += carry * carry
            conv_j = sum(g[i] * h[j - i] for i in range(j + 1)
                        if i <= K2 and j - i <= K2)
            carry = (conv_j + carry) >> 1

    N2 = Fraction(n2)
    pr(f"\n  {'j':>3}  {'E[c²]':>16}  {'E[c]²':>16}  {'Var[c]':>16}  {'Var exact':>20}")
    for j in range(K2 + 1):
        ec2 = Fraction(S_c2[j], n2)
        ec = Fraction(S_c[j], n2)
        var = ec2 - ec * ec
        pr(f"  {j:3d}  {str(ec2):>16}  {str(ec*ec):>16}  {str(var):>16}  {float(var):20.15f}")

    pr(f"\n  Testing Var[carry_j] = a·j + b pattern:")
    for j in range(2, K2 + 1):
        ec2 = Fraction(S_c2[j], n2)
        ec = Fraction(S_c[j], n2)
        var = ec2 - ec * ec
        slope = var / Fraction(j)
        pr(f"    j={j}: Var/j = {float(slope):.10f}")

    pr("\n" + "=" * 78)
    pr("  END")
    pr("=" * 78)


if __name__ == '__main__':
    main()
