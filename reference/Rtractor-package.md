# Rtractor: Complexity and Nonlinear Time Series Analysis for Physiological Signals

A staging- and signal-agnostic toolkit of complex systems measures for
physiological time series: entropy (sample, permutation, approximate),
fractal and multifractal analysis (Higuchi dimension, DFA, MFDFA),
Lyapunov exponents, multiscale entropy, and recurrence quantification
analysis. Core algorithms wrap established C/C++/Fortran reference
implementations rather than reimplementing them in pure R, to preserve
numerical parity with the original methods literature. Some wrapped
reference implementations are GPL-licensed (see inst/COPYRIGHTS), which
is why the package as a whole is licensed GPL (\>= 2) rather than MIT.
Designed to sit underneath other Circadia Lab and CoDe-Neuro Lab
packages (mrpheus, hypnoR, zeitR, dynR) as a shared complexity-metrics
layer, while remaining fully usable standalone on any numeric time
series.

## See also

Useful links:

- <https://rtractor.circadia-lab.uk>

- <https://github.com/circadia-bio/Rtractor>

- Report bugs at <https://github.com/circadia-bio/Rtractor/issues>

## Author

**Maintainer**: Lucas França <lucas.franca@northumbria.ac.uk>
([ORCID](https://orcid.org/0000-0003-0853-1319))

Authors:

- Lucas França <lucas.franca@northumbria.ac.uk>
  ([ORCID](https://orcid.org/0000-0003-0853-1319))

- Mario Leocadio-Miguel <mario.miguel@northumbria.ac.uk>
  ([ORCID](https://orcid.org/0000-0002-7248-3529))
