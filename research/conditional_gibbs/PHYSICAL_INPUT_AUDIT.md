# From conditional envelopes to mineral phase comparisons

## What was verified

PR #6's envelope theorem is sharp for a single Helmholtz branch on its closed sampled interval. The interval observations are feasible exactly when every pair obeys
`|y_i-y_j| <= L|V_i-V_j| + eps_i + eps_j`. Necessity follows from the triangle inequality. For sufficiency, this condition makes every lower cone lie below every upper cone by the triangle inequality. Each envelope is itself Lipschitz and lies in the given interval at each knot. Thus each endpoint of the minimum-energy enclosure is attained by an admissible function. Independent phase admissibility also makes the difference interval `[lo_A-hi_B, hi_A-lo_B]` sharp. Correlated errors can tighten it if a joint admissible set is independently justified.

The old implementation could accept mutually inconsistent samples because of a fixed energy tolerance. It could also skip a real cone crossing with two nearly parallel slopes or lose a volume interval to floating-point cancellation. The proposed implementation uses rational arithmetic for input consistency, the pressure tilt and minimization. It preserves the JSON numeric interface with outward-rounded energy endpoints and displays the sampled domain. Floating-point inputs represent their actual binary values. Use integers or `Fraction` for exact decimal inputs. Volume diagnostics are rounded for display. The lower sweep uses prefix and suffix extrema because all cones on each side of a cell have the same slope. The upper minimum is the smallest individual upper-cone minimum. Both require a linear number of rational operations.

Pressure work is positive `+PV` for compression-positive external pressure. Divide cell energy and volume by the same formula-unit count before comparing. Stress is intensive and needs no formula-unit divisor. With the exact SI charge and Avogadro constant, 1 GPa Å³ is approximately 0.006241509074460763 eV, 1 eV per formula unit is 96.4853321233100184 kJ per mole of those units, and 1 GPa cm³/mol is exactly 1 kJ/mol. [NIST SI definitions](https://www.nist.gov/pml/special-publication-330/sp-330-section-2)

## Original silica scans and conditional results

The authors' [r²SCAN archive](https://doi.org/10.5281/zenodo.14925354) supplies static cell energies and stresses. I verified its publisher MD5 and pinned SHA-256. The extractor selects the 26 plain `volsc` rows per branch and excludes rattled structures, separate relaxations and active-learning snapshots. The Si/O counts confirm eight formula units in coesite and four in each tridymite branch. The pressure diagnostic is minus one third the trace of the ASE stress tensor. The static database energy field contains neither zero-point nor thermal vibrational free energy. [Original study and supplement](https://materialsmodeling.org/assets/publications/SveSadWu26.pdf)

| Branch | Volume range (Å³/SiO₂) | Largest gap (Å³/SiO₂) | Largest sampled absolute stress pressure (GPa) | Zero-error energy secant floor (GPa) | Trial cap (GPa) |
|---|---:|---:|---:|---:|---:|
| Coesite | 29.080–37.735 | 0.346 | 23.545 | 22.754 | 24 |
| α-tridymite | 39.732–51.557 | 0.473 | 5.388 | 5.312 | 6 |
| β-tridymite | 42.603–55.282 | 0.507 | 38.124 | 36.656 | 39 |
| MC-tridymite | 39.280–50.970 | 0.468 | 3.999 | 3.887 | 4 |

The trial caps round up the larger sampled diagnostic per branch. They are explicit hypotheses. A finite set of pressures cannot establish an upper bound between knots. With zero knot error, α- and MC-tridymite have 3 and 5 decreases in successive secant slopes, respectively. Each decrease rules out exact convex interpolation of that series. Relaxation branch changes or numerical scatter are possible explanations, and the database does not establish which applies. The trial caps should not be promoted by treating the entire scan as convex.

Under the trial caps and a hypothetical absolute error of **1 meV/SiO₂ per knot**, at **0.25 GPa** the restricted static minimum of coesite minus that of α-tridymite lies between **−42.003 and −3.503 meV/SiO₂**. The β and MC intervals cross zero. The four-branch comparison abstains at 0 and 0.25 GPa, even with the knot errors set to zero. This calculation is conditional on the chosen caps, the identified structure branches, the error budget and their sampled domains.

An expanded sweep evaluates 21 exact pressure values from 0 to 1 GPa in 0.05 GPa increments and four trial knot errors. At 0.50 GPa and 1 meV/SiO₂ per knot, coesite separates from all three scanned tridymite branches. The upper endpoints of coesite minus α, β and MC tridymite are approximately −22.399, −3.750 and −7.991 meV/SiO₂. On this grid, the first four-branch separation is 0.45 GPa with zero error and 0.50 GPa at each of 0.5, 1 and 2 meV/SiO₂. These are decisions at tested pressures on the fixed sampled volume domains. They do not locate a physical phase boundary.

The broader replay adds the two α-quartz and two β-quartz series, α-cristobalite, keatite and stishovite. The pinned database contains 282 plain scan rows across eleven branches. At 0.50 GPa and 1 meV/SiO₂ per knot, both α-quartz scans have lower restricted static minima than coesite under their respective trial caps. The lower endpoints of coesite minus those scans are about +1.399 and +1.507 meV/SiO₂. The eleven-way comparison abstains at all 84 pressure/error scenarios. The earlier coesite–tridymite result therefore remains a statement about that four-branch subset.

The nominal Murnaghan reconstructions from [Demuth et al., Tables 6 and 8](https://doi.org/10.1088/0953-8984/11/19/306) provide a separate test with an analytically justified derivative bound **for each fitted EOS**. Internal pressure decreases monotonically with volume in these models, so the maximum absolute endpoint pressure bounds the entire roughly ±2% volume window. The original authors report restricted fit ranges. The replay treats the printed coefficients as exact model parameters and evaluates the model at 41 rational knots. At zero pressure, coesite minus orthorhombic `C2221` tridymite is enclosed by **[−22.202, −21.909] meV/SiO₂ in LDA** and **[116.494, 116.777] meV/SiO₂ in GGA**. The opposite signs show that tightly bounding the chosen EOS does not settle the material ordering. The model branch, old functionals and static temperature also differ from the retrospective comparison.

## How dense and precise do the inputs need to be?

Let `m_i = min_j(y_j + p V_j)` be the smallest nominal tilted knot value for phase `i`. With a uniform sample error `eps_i`, its upper minimum is at most `m_i + eps_i`. Its lower minimum is at least `m_i - eps_i - q_i h_i`, where `h_i` is its largest knot spacing, `p` is external pressure converted to eV/Å³, and `q_i = (L_i²-p²)/(2 L_i)` for `|p| < L_i`, zero otherwise. In a cell of width `h`, the two adjacent lower cones give deficit ramps `(L-p)x` and `(L+p)(h-x)`. Their smaller value peaks at their crossing and equals `q h`. For `|p| >= L`, the tilted curve is monotone and a domain endpoint minimizes it. Thus a sufficient strict test for A below B is
`m_B - m_A > eps_A + eps_B + q_B h_B`. Only B's grid penalty enters this one-sided test. The full envelope can certify some cases where this shortcut is inconclusive.

At zero pressure the original MC-tridymite scans have a sampled minimum 16.217 meV/SiO₂ below coesite. If that gap persists on refinement, a 24 GPa coesite cap and 1 meV/SiO₂ error per phase would require coesite spacing below **0.190 Å³/SiO₂** by the sufficient test. The current maximum is 0.346 Å³/SiO₂. A better-supported local slope bound or tighter lower envelope can reduce the sampling needed.

[Zhong, Li & John, Table 1](https://www.nature.com/articles/s41524-026-02056-3.pdf) gives equilibrium **enthalpies** at 300 K and 1 bar. Both historical inversions reproduce by exact decimal subtraction. These numbers size a diagnostic and cannot be used as a finite-temperature Helmholtz input.

| A − B | Raw gap (kJ/mol as reported) | Reference gap (kJ/mol as reported) | Minimum equal phase-shift budget to reach raw tie | Per-phase error ceiling to resolve reference gap, with no grid penalty |
|---|---:|---:|---:|---:|
| Coesite − tridymite | −0.29 | +0.06 | 0.145 | <0.030 |
| Clinoenstatite − enstatite | +0.53 | −0.89 | 0.265 | <0.445 |

For silica, the 0.06 kJ/mol reference gap is about 0.622 meV/SiO₂. Allocating half to grid coverage at an illustrative 1 GPa internal-pressure cap and zero external pressure gives spacing below **0.100 Å³/SiO₂**, and leaves less than 0.015 kJ/mol per phase for errors. The pyroxene table lists `MgSiO3` but reports clino/ortho volumes about twice the ~31.3 cm³/mol measured per MgSiO₃ by [Angel & Hugh-Jones](https://doi.org/10.1029/94JB01750). Keep the historical numbers in their reported basis until the enthalpy normalization is confirmed from the original workflow.

The original [Zhong training archive](https://osf.io/784bv/) has 948 clinoenstatite reference configurations and 864 enstatite reference configurations, with 8 and 16 MgSiO₃ units per cell respectively. For example, reference clinoenstatite frames 740 and 68 differ by only 0.0000833 Å³/MgSiO₃ but by 2.753 eV/MgSiO₃. Enstatite reference frames 48 and 673 give 0.0000696 Å³/MgSiO₃ and 3.471 eV/MgSiO₃. Such configurations sample structural degrees of freedom at nearly the same volume. Feeding the archive directly to an envelope for a constrained-relaxed energy would demand errors or slopes far above the resolution budget. The frame normalization does not resolve the published Table 1 molar basis. The companion OSF example contains two trained models and a quartz LAMMPS input, with no Table 1 pyroxene normalization script.

## Inputs still required for physical certification

1. **A defined energy target and pointwise uncertainty.** Use a consistent cell normalization, energy reference and structural identity at every volume. Numerical cutoffs, k meshes, PAW choices, smearing and constrained relaxation need a phase-paired convergence study. For a physical Gibbs claim, bound the difference between the electronic model and the physical energy and account for zero-point, thermal, and where relevant anharmonic terms. An SCF stopping criterion, fit residual or average model error gives no uniform physical enclosure. Even a perfect static replay leaves these terms open.
2. **A global derivative bound on the branch.** One route is a verified bound `M` on `|F''|` across a same-branch interval, with independently bounded pointwise derivative errors `delta_i`. Then `L = max_i(|d_i| + delta_i) + M h/2` works by comparing to the nearest knot. At finite temperature the derivative is that of the full Helmholtz energy and includes vibrational pressure. A fitted EOS can bound its own derivative. A bounded EOS energy residual alone can transfer minimum values by its uniform bound, but cannot transfer the derivative.
3. **Exterior minima and phase completeness.** A positive lower curvature bound on an entire physical branch and derivative intervals that bracket `-p` at the sampled endpoints would confine its minimizer. Otherwise establish a lower bound on all omitted volumes or extend the sampled domain. Include other plausible polymorphs for a stability claim. A local minimum or a sampled endpoint stress sign is insufficient by itself.
4. **Matched pyroxene inputs.** [Yu et al.](https://doi.org/10.1029/2009JB006329) covers low- and high-pressure clinoenstatite, and [the orthoenstatite study](https://doi.org/10.1002/2017GL075647) supplies a separate ortho branch. Separate curves need a common electronic energy reference, controlled structure relaxations, and matched finite-temperature contributions. Yu et al.'s reported cutoff and k-mesh convergence scale of 0.00015 Ry/atom is about 10.2 meV/MgSiO₃ before model and thermal uncertainty. The diffraction EOS independently constrains compression and also reports a change in the orthoenstatite compression mechanism near 4 GPa.

## Replay

Run from the repository root.

```sh
python -m unittest discover -s tests -v
python -m research.conditional_gibbs.sensitivity
python -m research.conditional_gibbs.scan_audit /path/to/reference-structures.db \
  --zhong-training /path/to/RawTrainingData.zip --pressure-sweep --all-silica-branches
```

The last command reads the [original SQLite database](https://zenodo.org/records/14925354/files/reference-structures.db) in read-only mode and requires a matching SHA-256. The optional argument also audits the [original training archive](https://osf.io/download/uqhdc/). Both computations read published numbers and do not execute DFT.
