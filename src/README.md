# src/

Native (C/C++/Fortran) source for Rtractor, wrapped via Rcpp.

- `dfa.cpp`      — direct C++ port of PhysioNet's `dfa.c` (GPL-2-or-later;
  see `inst/COPYRIGHTS`). Validated to match the compiled reference binary
  exactly on synthetic test data.
- `higuchi.cpp`  — clean-room reimplementation of Higuchi (1988) fractal
  dimension. Validated against a transliteration of the MATLAB reference
  on synthetic test data (max abs. difference ~1e-11).

Planned, not yet ported:
- `fractal.cpp` (MFDMA, Chhabra-Jensen multifractal spectrum)
- `lyapunov.cpp`, `multiscale.cpp`, `rqa.cpp`, `embed.cpp`

Package license is `GPL (>= 2)` (see root `LICENSE` and `inst/COPYRIGHTS`)
because some wrapped reference implementations are themselves GPL-licensed.
Reference sources being wrapped or ported should be dropped in
`dev/reference/` (untracked scratch — see root `.gitignore`) rather than
committed verbatim, unless it's the actual wrapped source (like `dfa.c`),
in which case attribution goes in `inst/COPYRIGHTS` and the original
copyright notice is preserved in a header comment in the corresponding
`.cpp` file.
