# frust

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14205558.svg)](https://doi.org/10.5281/zenodo.14205558)

A C++ implementation of the stochastic series expansion quantum Monte Carlo method [[1]](#1) with abstract loop updates, which was used in the works

- J. L. Jiménez, S. P. G. Crone, E. Fogh, M. E. Zayed, R. Lortz, E. Pomjakushina, K. Conder, A. M. Läuchli, L. Weber, S. Wessel, A. Honecker, B. Normand, C. Rüegg, P. Corboz, H. M. Rønnow, and F. Mila,
  *A quantum magnetic analogue to the critical point of water,*
  Nature **592**, 370-375 (2021)
- L. Weber, A. Honecker, B. Normand, P. Corboz, F. Mila, and S. Wessel,
  *Quantum Monte Carlo simulations in the trimer basis: first-order transitions and thermal critical points in frustrated trilayer magnets*,
  SciPost Phys. **12**, 054 (2022)
- L. Weber, N. Caci, and S. Wessel,
  *Cluster quantum Monte Carlo study of two-dimensional weakly coupled frustrated trimer antiferromagnets*,
  Phys. Rev. B **106**, 035141 (2022)
- L. Weber, A. Y. D. Fache, F. Mila, and S. Wessel,
  *Thermal critical points from competing singlet formations in fully frustrated bilayer antiferromagnets*,
  Phys. Rev. B **106**, 235128 (2022)
- L. Weber, E. Viñas Boström, M. Claassen, A. Rubio, and D. M. Kennes,
  *Cavity-renormalized quantum criticality in a honeycomb bilayer antiferromagnet*,
  Commun Phys **6**, 247 (2023)

There is a Julia version of this code, [StochasticSeriesExpansion.jl](https://github.com/lukas-weber/StochasticSeriesExpansion.jl), which is better documented, more user-friendly and also slightly faster than this C++ version. It is therefore recommended to use the Julia code.

This repository is provided without maintenance or support, solely for the interest of the community. If you still decide to use it for some research, kindly cite the DOI.

## References
<a id="1">[1]</a> A. W. Sandvik, Phys. Rev. B **59**, R14157(R) (1999)
