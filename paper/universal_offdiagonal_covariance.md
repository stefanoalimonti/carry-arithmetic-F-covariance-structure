# Exact Covariance Structure of Binary Carry Chains

**Author:** Stefano Alimonti
**Affiliation:** Independent Researcher
**Date:** March 2026

---

**Abstract.**
We establish the exact second-order structure of the carry chain arising from binary multiplication with fixed most-significant bits. Six results are proved: (i) the column expectation E[conv_j] = (j+3)/4 (Lemma A); (ii) the Parity Lemma, that total_j mod 2 is uniform for all j ≥ 1 (Lemma B); (iii) the carry expectation E[carry_j] = (j−1)/4 (Theorem 3); (iv) the universal off-diagonal covariance Cov(carry_j, g_i h_{j−i}) = 1/8 for all off-diagonal pairs (Lemma D); (v) Cov(carry_j, conv_j) = (j−1)/8 for all odd j ≥ 3 (Main Theorem); and (vi) the exact inductive step ΔCov = 1/4 for odd-to-odd transitions (Corollary). For even j, a diagonal correction ε_j > 0 arises from the failure of the Parity Lemma under conditioning; we tabulate these corrections through j = 12.

**Overview.** Section 1 introduces the binary multiplication model. Section 2 establishes three fundamental identities: the column expectation (Lemma A), the Parity Lemma (Lemma B), and the carry expectation (Theorem 3). Section 3 proves the universal off-diagonal covariance (Lemma D). Section 4 derives the main covariance formula and its corollary. Section 5 treats even-j corrections. Sections 6–7 present computational verification and discussion.

---

## 1. Setup and Notation

We study the carry chain arising from binary multiplication of two n-bit integers G = g₀ g₁ ⋯ g_{n−1} and H = h₀ h₁ ⋯ h_{n−1} in base 2.

**Bit model.** Fix the most-significant bits g₀ = h₀ = 1. For i ≥ 1, each g_i and h_i is drawn independently from Bernoulli(1/2).

**Convolution sums.** For each column j ≥ 0:

$$\text{conv}_j = \sum_{i=0}^{j} g_i \, h_{j-i}.$$

**Carry recurrence.** Set carry₀ = 0, and define

$$\text{carry}_{j+1} = \left\lfloor \frac{\text{conv}_j + \text{carry}_j}{2} \right\rfloor.$$

Write total_j = conv_j + carry_j for the full column sum. The output bit at position j is total_j mod 2, and the carry propagated forward is ⌊total_j / 2⌋.

**Terminology.** A bit-product term g_i h_{j−i} appearing in conv_j is *off-diagonal* if i ≠ j − i (equivalently, j is odd, or j is even and i ≠ j/2). It is *diagonal* if i = j − i = j/2.

---

## 2. Fundamental Identities

### 2.1 Lemma A (Column Expectation)

**Lemma A.** *For j ≥ 1,*

$$E[\text{conv}_j] = \frac{j+3}{4}.$$

*Proof.* The convolution sum conv_j contains j + 1 terms. Two are border terms involving the fixed MSBs: g₀ h_j = h_j with E[h_j] = 1/2, and g_j h₀ = g_j with E[g_j] = 1/2. The remaining j − 1 terms are interior products g_i h_{j−i} with both factors random, so E[g_i h_{j−i}] = 1/4. Therefore:

$$E[\text{conv}_j] = 2 \cdot \frac{1}{2} + (j-1) \cdot \frac{1}{4} = 1 + \frac{j-1}{4} = \frac{j+3}{4}. \qquad \square$$

### 2.2 Lemma B (Parity Lemma)

**Lemma B** (Parity Lemma). *For all j ≥ 1,*

$$P(\text{total}_j \text{ is odd}) = \frac{1}{2}.$$

*Proof.* Write total_j = g_j · h₀ + W_j, where

$$W_j = g_0 \cdot h_j + \sum_{i=1}^{j-1} g_i h_{j-i} + \text{carry}_j.$$

Since h₀ = 1, we have total_j = g_j + W_j. The digit g_j first enters the computation at column j: it does not appear in conv₀, …, conv_{j−1} and hence is independent of carry_j. It is also independent of every term in W_j. Since g_j ∼ Ber(1/2):

$$P(\text{total}_j \text{ odd}) = \tfrac{1}{2}\,P(W_j \text{ even}) + \tfrac{1}{2}\,P(W_j \text{ odd}) = \frac{1}{2}. \qquad \square$$

**Remark** (General base). For base b ≥ 2, the analogous statement is that total_j mod b is uniform on {0, …, b−1}. The proof uses the same decomposition total_j = g_j · h₀ + W_j: as g_j ranges over {0, …, b−1}, the map a ↦ a · h₀ mod b is a bijection if and only if gcd(h₀, b) = 1. For base 2 this is automatic since h₀ = 1. For general bases, the coprimality condition gcd(g₀ h₀, b) = 1 suffices, since either g_j · h₀ or g₀ · h_j provides the required uniform residue. When both gcd(g₀, b) > 1 and gcd(h₀, b) > 1, a character-sum argument using the pair (g_j, h_j) simultaneously shows that uniform distribution still holds provided gcd(g₀, b) · gcd(h₀, b) < b, which is satisfied for all valid MSB pairs in bases b ≤ 12 and for all prime-power bases.

### 2.3 Theorem 3 (Carry Expectation)

**Theorem 3.** *For all j ≥ 1,*

$$E[\text{carry}_j] = \frac{j-1}{4}.$$

*Proof.* By induction. **Base:** carry₁ = ⌊g₀ h₀ / 2⌋ = ⌊1/2⌋ = 0 = (1 − 1)/4. ✓

**Inductive step.** Suppose E[carry_j] = (j − 1)/4. Since carry_{j+1} = (total_j − (total_j mod 2)) / 2:

$$E[\text{carry}_{j+1}] = \frac{E[\text{total}_j] - E[\text{total}_j \bmod 2]}{2}.$$

By Lemma A and the inductive hypothesis, E[total_j] = E[conv_j] + E[carry_j] = (j + 3)/4 + (j − 1)/4 = (j + 1)/2. By the Parity Lemma, E[total_j mod 2] = P(total_j odd) = 1/2. Hence:

$$E[\text{carry}_{j+1}] = \frac{(j+1)/2 - 1/2}{2} = \frac{j}{4}. \qquad \square$$

---

## 3. Universal Off-Diagonal Covariance (Lemma D)

**Lemma D.** *For all 1 ≤ i ≤ j − 1 with i ≠ j − i:*

$$\text{Cov}(\text{carry}_j,\; g_i h_{j-i}) = \frac{1}{8}.$$

*Proof.* Since g_i h_{j−i} ∈ {0, 1} with P(g_i h_{j−i} = 1) = 1/4:

$$\text{Cov}(\text{carry}_j,\; g_i h_{j-i}) = \frac{1}{4}\,\delta E[c_j],$$

where δE[c_j] = E[carry_j | g_i = 1, h_{j−i} = 1] − E[carry_j] is the perturbation in the carry expectation caused by conditioning. We establish δE[c_j] = 1/2 in four steps.

**Step 1: The Parity Lemma survives conditioning.** Let a = min(i, j − i) and b = max(i, j − i). Since i ≠ j − i, we have a < b, and the conditioning fixes at most one bit at any given index. Specifically: at column k, the fresh bit g_k is conditioned only if k = i, but then h_k is free (since k = i ≠ j − i means h_k is unconditioned), and h_k can serve as the uniform parity-flipping variable. Similarly, if k = j − i, then g_k is free. At all other columns, g_k is free. Therefore the Parity Lemma holds at every column k ≥ 1 under the conditioning.

**Step 2: Perturbation recursion.** Let δE[c_k] = E[carry_k | g_i = 1, h_{j−i} = 1] − E[carry_k]. Because the Parity Lemma holds under conditioning, the mod-2 correction is unchanged, giving the exact recursion:

$$\delta E[c_{k+1}] = \frac{\delta E[\text{conv}_k] + \delta E[c_k]}{2},$$

with δE[c₀] = 0.

**Step 3: Perturbation profile δE[conv_k].** Conditioning g_i = 1 shifts E[g_i] from 1/2 to 1; conditioning h_{j−i} = 1 shifts E[h_{j−i}] from 1/2 to 1. The term g_i h_{k−i} in conv_k gains +1/2 (border, at k = i) or +1/4 (interior, at k > i); similarly for h_{j−i}. The full profile is:

| Range | δE[conv_k] |
|-------|------------|
| k < a | 0 |
| k = a | 1/2 |
| a < k < b | 1/4 |
| k = b | 3/4 |
| k > b | 1/2 |

At k = a, the first conditioned bit enters as a border term (+1/2). For a < k < b, it contributes as an interior term (+1/4). At k = b, the second conditioned bit enters as a border term (+1/2) while the first continues as interior (+1/4), giving 3/4. For k > b, both contribute interior terms (+1/4 each).

**Step 4: Exact evaluation via generating function.** Define T_k = 2^k · δE[c_k]. The recursion yields T_{k+1} = T_k + 2^k · δE[conv_k], so:

$$T_j = \sum_{k=0}^{j-1} 2^k \,\delta E[\text{conv}_k].$$

Substituting the profile from Step 3:

$$T_j = 2^a \cdot \frac{1}{2} + \frac{1}{4}\sum_{k=a+1}^{b-1} 2^k + \frac{3}{4} \cdot 2^b + \frac{1}{2}\sum_{k=b+1}^{j-1} 2^k.$$

Evaluating the geometric sums:

$$T_j = 2^{a-1} + \frac{2^b - 2^{a+1}}{4} + \frac{3 \cdot 2^b}{4} + \frac{2^j - 2^{b+1}}{2}.$$

Collecting terms:

$$T_j = 2^{a-1} + \frac{4 \cdot 2^b - 2^{a+1}}{4} + \frac{2^j - 2^{b+1}}{2} = 2^{a-1} + 2^b - 2^{a-1} + 2^{j-1} - 2^b = 2^{j-1}.$$

Therefore δE[c_j] = T_j / 2^j = 1/2, and:

$$\text{Cov}(\text{carry}_j,\; g_i h_{j-i}) = \frac{1}{4} \cdot \frac{1}{2} = \frac{1}{8}. \qquad \square$$

---

## 4. The Main Covariance Theorem

**Main Theorem.** *For all odd j ≥ 3:*

$$\text{Cov}(\text{carry}_j, \text{conv}_j) = \frac{j-1}{8}.$$

*Proof.* Expand by linearity of covariance:

$$\text{Cov}(\text{carry}_j, \text{conv}_j) = \sum_{i=0}^{j} \text{Cov}(\text{carry}_j,\; g_i h_{j-i}).$$

We classify each term by type:

- **Border terms (i = 0 and i = j):** Cov(carry_j, g₀ h_j) = Cov(carry_j, h_j) = 0, because h_j first enters the computation at column j and is independent of carry_j (which depends only on columns 0, …, j − 1). Similarly, Cov(carry_j, g_j h₀) = Cov(carry_j, g_j) = 0.

- **Off-diagonal terms (1 ≤ i ≤ j − 1):** Since j is odd, i ≠ j − i for every such i (there is no diagonal term). By Lemma D, each covariance equals 1/8.

There are exactly j − 1 off-diagonal terms, giving:

$$\text{Cov}(\text{carry}_j, \text{conv}_j) = (j-1) \cdot \frac{1}{8} = \frac{j-1}{8}. \qquad \square$$

**Corollary** (Exact inductive step). *For all odd j ≥ 3:*

$$\text{Cov}(\text{carry}_{j+2}, \text{conv}_{j+2}) - \text{Cov}(\text{carry}_j, \text{conv}_j) = \frac{1}{4}.$$

*Proof.* Direct computation:

$$\frac{(j+2) - 1}{8} - \frac{j - 1}{8} = \frac{j+1 - j + 1}{8} = \frac{2}{8} = \frac{1}{4}. \qquad \square$$

This confirms that the ΔCov = 1/4 induction step observed computationally is exact for odd-to-odd transitions. The constancy of the step size reflects the fact that each increase of j by 2 adds exactly two new off-diagonal terms to conv_j, each contributing Cov = 1/8 by Lemma D.

---

## 5. Even-j Corrections

For even j, the convolution conv_j contains a diagonal term g_{j/2} h_{j/2} = g_{j/2}², for which Lemma D does not apply. The decomposition becomes:

$$\text{Cov}(\text{carry}_j, \text{conv}_j) = \underbrace{(j-2) \cdot \frac{1}{8}}_{\text{off-diagonal}} + \underbrace{\text{Cov}(\text{carry}_j,\; g_{j/2} h_{j/2})}_{\text{diagonal}}.$$

For the off-diagonal terms, Lemma D applies as before. The diagonal term requires separate treatment. Conditioning on g_{j/2} = 1 and h_{j/2} = 1 fixes *both* free bits at index j/2 simultaneously. At column j/2, neither g_{j/2} nor h_{j/2} is available to flip the parity of total_{j/2}, so the Parity Lemma fails at that single column. The carry perturbation recursion (Step 2 of the Lemma D proof) acquires a nonzero error at column j/2, which propagates forward and decays geometrically but does not vanish by position j. This introduces a correction:

$$\varepsilon_j = \text{Cov}(\text{carry}_j,\; g_{j/2} h_{j/2}) - \frac{1}{8}.$$

The correction ε_j is strictly positive and decays exponentially with the gap j/2 between the conditioning column and the target.

**Table of verified diagonal corrections.**

| j | ε_j | Cov(carry_j, conv_j) | (j−1)/8 | Difference |
|---|-----|----------------------|---------|------------|
| 2 | 1/16 | 3/16 | 1/8 | 1/16 |
| 4 | 1/32 | 13/32 | 3/8 | 1/32 |
| 6 | 1/128 | 81/128 | 5/8 | 1/128 |
| 8 | 1/256 | 225/256 | 7/8 | 1/256 |
| 10 | 3/4096 | 4611/4096 | 9/8 | 3/4096 |
| 12 | 11/32768 | 45067/32768 | 11/8 | 11/32768 |

The corrections ε₂ = 1/16 and ε₄ = 1/32 have simple forms, but for j ≥ 10 the numerators (3, 11, …) do not follow a simple closed-form pattern. The overall scaling ε_j ≈ C · 2^{−j} is consistent with the geometric convergence rate of the perturbed carry chain after the Parity Lemma breaks at a single column.

For all even j:

$$\text{Cov}(\text{carry}_j, \text{conv}_j) = \frac{j-1}{8} + \varepsilon_j,$$

so (j − 1)/8 serves as a lower bound with exponentially small correction. The numerator pattern for j ≥ 10 is non-trivial (ε₁₀ = 3/4096, ε₁₂ = 11/32768), reflecting the increasingly complex interaction between the diagonal perturbation and the carry chain dynamics.

---

## 6. Computational Verification

All results have been verified by exact rational arithmetic, enumerating all 2^{2K} bit configurations for K up to 12 (corresponding to 2^{24} = 16,777,216 configurations).

**Enumeration method.** For a given K, we enumerate all 2^{2K} assignments of the K free bits in each of G and H (with g₀ = h₀ = 1 fixed), compute the carry chain exactly using integer arithmetic, and accumulate the required moments and cross-moments as exact rational numbers.

**Column and carry expectations.** For all j ≤ 12, E[conv_j] = (j + 3)/4 and E[carry_j] = (j − 1)/4 hold exactly.

**Off-diagonal universality.** At K = 8 (65,536 configurations per bit-pair check), all off-diagonal covariances Cov(carry_j, g_i h_{j−i}) for 1 ≤ i ≤ j − 1, i ≠ j − i, j ≤ 8, evaluate to exactly 1/8.

**Diagonal anomaly.** All diagonal covariances Cov(carry_j, g_{j/2} h_{j/2}) for even j ≤ 12 strictly exceed 1/8, with corrections matching the table in §5.

**Odd-j formula.** For all odd j with 3 ≤ j ≤ 11, Cov(carry_j, conv_j) = (j − 1)/8 exactly.

**ΔCov step.** The step Cov(carry_{j+2}, conv_{j+2}) − Cov(carry_j, conv_j) = 1/4 holds exactly for all consecutive odd pairs (3, 5), (5, 7), (7, 9), (9, 11).

**Parity Lemma under conditioning.** The survival of the Parity Lemma under conditioning on off-diagonal bit pairs (Step 1 of the Lemma D proof) has been verified for all off-diagonal pairs with j ≤ 10: the conditional distribution of total_k mod 2 remains exactly Ber(1/2) at every column, confirming the theoretical argument.

**Even-j diagonal corrections.** The diagonal covariance values Cov(carry_j, g_{j/2} h_{j/2}) for even j ≤ 12 all exceed 1/8, with the excess matching the ε_j values in §5 to full precision.

---

## 7. Discussion

**Independence from the trace anomaly conjecture.** The results in this paper — Lemmas A, B, D, Theorem 3, the Main Theorem, and the Corollary — are proved by elementary means (induction, conditioning, geometric sums). They do not depend on the c₁ = π/18 conjecture, any spectral assumptions, or any unproved distributional hypotheses about the carry chain. Every step relies only on the independence structure of the bit model and the deterministic carry recurrence.

**Connection to the non-Markovian correction.** The covariance Cov(carry_j, conv_j) quantifies the non-Markovian correction δ_NM in the trace anomaly analysis. Specifically, the variance of carry_j decomposes as:

$$\text{Var}(\text{carry}_{j+1}) = \frac{1}{4}\bigl(\text{Var}(\text{conv}_j) + \text{Var}(\text{carry}_j) + 2\,\text{Cov}(\text{carry}_j, \text{conv}_j) - \text{Var}(\text{total}_j \bmod 2)\bigr).$$

The exact evaluation of the covariance term removes one unknown from this recursion, constraining the carry variance growth and thereby the trace anomaly.

**Odd vs. even dichotomy.** The clean formula Cov = (j−1)/8 holds for all odd j because the off-diagonal universality of Lemma D applies to every interior term and no diagonal term exists. For even j, the single diagonal term g_{j/2} h_{j/2} breaks the Parity Lemma under conditioning at column j/2, introducing a correction that decays geometrically with the gap j/2 between the conditioning point and the target column j. This odd/even dichotomy is a structural feature of binary multiplication — it arises from the parity of the convolution length, not from any approximation in the proof method.

**Open question.** Whether the even-j corrections ε_j admit a closed-form expression remains open. The values through j = 12 are consistent with ε_j = f(j) · 2^{−j} for some slowly-varying function f, but the numerator sequence (1, 1, 1, 1, 3, 11, …) does not match any standard combinatorial sequence.

**The perturbation method.** The proof of Lemma D introduces a general technique: condition on a bit pair, track the perturbation to column expectations via an exact linear recursion (enabled by the Parity Lemma), and evaluate a weighted geometric sum. This method extends naturally to higher-order covariances and to other bases, suggesting that analogues of the universal off-diagonal covariance may hold in broader settings.

**Companion papers.** The spectral theory of the carry operator is developed in [A] (*Spectral Theory of Carries in Positional Multiplication*). The trace anomaly and its connection to the c₁ conjecture are analyzed in [E] (*The Trace Anomaly of Binary Multiplication*). The present paper provides the exact covariance inputs required by both analyses. In particular, the Main Theorem supplies the Cov(carry_j, conv_j) term that governs the non-Markovian correction in the variance recursion of [E], and the Parity Lemma provides the key parity input for the spectral analysis of [A].

---

## 8. References

1. P. Diaconis and J. Fulman, "Carries, shuffling, and symmetric functions," *Advances in Applied Mathematics*, vol. 43, no. 2, pp. 176–196, 2009.

2. J. M. Holte, "Carries, combinatorics, and an amazing matrix," *The American Mathematical Monthly*, vol. 104, no. 2, pp. 138–149, 1997.

3. D. E. Knuth, *The Art of Computer Programming, Volume 2: Seminumerical Algorithms*, 3rd ed. Addison-Wesley, 1997.

4. [A] Companion paper: *Spectral Theory of Carries in Positional Multiplication*, this series.

5. [E] Companion paper: *The Trace Anomaly of Binary Multiplication*, this series.

6. [P2] Companion paper: *The Sector Ratio in Binary Multiplication: From Markov Failure to Transcendence*, this series.
