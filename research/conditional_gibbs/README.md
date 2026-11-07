# Conditional Gibbs intervals — GeoForm mathematical research v10

This is an experimental, self-contained **mathematical method**, not an ab initio mineral predictor. It introduces a no-spline alternative to the separately preserved v9 QHA interpolator, without changing the nine-crystal phonon catalog.

## Theorem: sharp conditional intervals

Assume the true molar/formula Helmholtz energy g(V) is defined throughout the *closed sampled volume interval* and obeys an **independently justified global** Lipschitz bound |g(V)-g(W)| <= L|V-W|. Also assume independently justified intervals g(V_i) in [y_i-eps_i,y_i+eps_i]. Define

- lower(V) = max_i [y_i-eps_i-L|V-V_i|]
- upper(V) = min_i [y_i+eps_i+L|V-V_i|].

For any supported pressure P, write p = 0.006241509074460763 P, with pressure GPa, volume A^3 per formula, and energy eV per formula. Then

min_V [lower(V)+pV] <= min_V [g(V)+pV] <= min_V [upper(V)+pV].

The upper and lower envelopes are themselves L-Lipschitz and, when the pairwise interval consistency checks hold, compatible with the supplied sample intervals. Therefore these are mathematically sharp worst-case bounds among functions obeying the assumptions. Because the envelopes are piecewise linear, scanning sample knots and pairwise affine intersections yields the extrema without assuming cubic-spline accuracy, convexity, or a unique equilibrium volume.

**No sample-only procedure proves a global physical L**. The source labels and stated L/error bounds are user assertions; neither the tests nor this code authenticates them. A DFT or other external physical certificate for the pressure/energy derivative bound is necessary for a non-illustrative use. The result covers only a finite supplied volume domain, no unprovided competing phases, and no synthesis kinetics.

## Connection to our independent research

GeoForm v9's 17-phase, six-formula retrospective analysis identified 2 erroneous same-formula orderings among 22 pairs: coesite/tridymite, and clinoenstatite/enstatite. Exact cancellation proves a formula-constant correction can never reverse them. A retrospective LP shows the corresponding minimum L-infinity phase-specific corrections to a tie are 0.145 and 0.265 kJ/mol. These numbers use the known reference values and are **not** prospective forecasts.

This nonparametric conditional-bounds implementation adds a separate decision rule: a winner is reported only when its *upper* minimum-energy bound lies below every comparator's *lower* bound. Otherwise it abstains. Source references, formula and temperature must match before comparison.

## Reproducibility and release boundary

Run from the repository root:

    python -m unittest discover -s tests -p test_conditional_gibbs.py -v

The tests include 200 analytically minimized quadratic cases, exact linear limits, contradictory input, missing provenance and overlapping-interval abstention.

The standalone GeoForm Physics v9 package remains distinct. Merging a clearly scoped mathematical methods module does not establish first-principles energies, true material phase stability, experimental spectra, or deployment of the owner-private GeoForm Site. The nine-case calculated phonon catalog and its withheld-unstable cases remain unchanged.

Background: https://phonopy.github.io/phonopy/qha.html
Materials Project phase methodology: https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/phase-diagrams-pds
