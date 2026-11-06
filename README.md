# Initial laws with all moments that stay multimodal under stable evolution

Yitzchak Shmalo

Paper: [paper/main.pdf](paper/main.pdf) (source [paper/main.tex](paper/main.tex)).

Hasebe and Ueda (ALEA 15 (2018) 353-374, arXiv:1710.08240, Conjecture 5.5) conjectured that if $\mu$ has a finite
moment of order $2+\alpha$, then the law at time $t$ of a classical, or free, strictly $\alpha$-stable process
started from $\mu$ is unimodal for all large $t$. The paper shows that the conjecture fails in both settings for
the positive stable processes of index $1/2$:

- for $\mu=\sum_{n\ge1}Z^{-1}2^{-n^2}\delta_{-2^n}$, the law $\mu*L_t$ (with $L_t$ the Levy distribution) is not
  unimodal for any $t>0$;
- for $\mu_0=(1-W)\delta_0+\sum_{n\ge1}2^{-n}A_n^{-n}\delta_{-A_n}$ with $A_n=2^{2^n}$, the free convolution
  $\mu_0\boxplus T_t$ with the positive free $1/2$-stable law is not unimodal for any $t>0$.

Both initial laws have finite moments of every order. The free case uses the global inversion theorem of Huang and
Wang (arXiv:2203.00421).

## Code

Python 3.14, standard library only. `code/check.py` writes `code/check.out` next to itself. It checks the
constants in the two proofs in exact rational arithmetic, evaluates the derivative of the density of $\mu*L_t$ at
the test points for $t=0.5,1,3,10$, and, for the free case, computes the subordination function on the gap
$(-256,-16)$ and at the two flanking atoms for $t=1$ and $t=1.7$ (Section 4 of the paper).

## Funding

The research presented in this paper was supported by the European Research Council (ERC) under the European
Union's Horizon Europe research and innovation programme (grant agreement No. 101041711), by the Simons Foundation
as part of the Collaboration on the Mathematical and Scientific Foundations of Deep Learning, by Heights Labs, by
Convex Nexus Capital, by the Israel Science Foundation (grant number 2258/19), and by the Israel Science Foundation
(ISF Grant 4101/25). Computational resources were provided in part by Google Cloud through its education credits
program (credit 530550103).
