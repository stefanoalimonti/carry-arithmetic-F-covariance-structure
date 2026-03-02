# carry-arithmetic-F-covariance-structure

**Exact Covariance Structure of Binary Carry Chains**

*Author: Stefano Alimonti* · [ORCID 0009-0009-1183-1698](https://orcid.org/0009-0009-1183-1698)

## Main Result

Complete second-order structure of the binary multiplication carry chain in six results:

- **Lemma A**: Column expectation $E[\text{conv}_j] = (j+3)/4$
- **Lemma B** (Parity Lemma): $\text{total}_j \bmod 2$ is uniform for base 2
- **Theorem 3**: Carry expectation $E[\text{carry}_j] = (j-1)/4$
- **Lemma D**: Universal off-diagonal covariance $\text{Cov}(\text{carry}_j, g_i h_{j-i}) = 1/8$
- **Main Theorem**: $\text{Cov}(\text{carry}_j, \text{conv}_j) = (j-1)/8$ for all odd $j \geq 3$
- **Corollary**: Inductive step $\Delta\text{Cov} = 1/4$ exact for odd-to-odd transitions

## Repository Structure

```
paper/universal_offdiagonal_covariance.md   The paper
experiments/
  F01_covariance_proof.py                   Core covariance proof verification
  F02_covariance_extended.py                Extended covariance computation
  F03_fast_covariance.py                    Fast covariance algorithm
  F04_induction_proof.py                    Induction step verification
  F05_pairing_analysis.py                   Pairing structure analysis
  F06_complete_proof.py                     Complete proof verification
  F07_ulc_covariance.py                     ULC covariance computation
```

## Reproduction

```bash
pip install numpy
python experiments/F01_covariance_proof.py
# ... through F07
```

## Dependencies

- Python >= 3.8, NumPy

## Companion Papers

| Label | Title | Repository |
|-------|-------|------------|
| [A] | Spectral Theory of Carries | [`carry-arithmetic-A-spectral-theory`](https://github.com/stefanoalimonti/carry-arithmetic-A-spectral-theory) |
| [E] | The Trace Anomaly of Binary Multiplication | [`carry-arithmetic-E-trace-anomaly`](https://github.com/stefanoalimonti/carry-arithmetic-E-trace-anomaly) |
| [P2] | The Sector Ratio in Binary Multiplication | [`carry-arithmetic-P2-sector-ratio`](https://github.com/stefanoalimonti/carry-arithmetic-P2-sector-ratio) |

### Citation

```bibtex
@article{alimonti2026covariance,
  author  = {Alimonti, Stefano},
  title   = {Exact Covariance Structure of Binary Carry Chains},
  year    = {2026},
  note    = {Preprint},
  url     = {https://github.com/stefanoalimonti/carry-arithmetic-F-covariance-structure}
}
```

## License

Paper: CC BY 4.0. Code: MIT License.
