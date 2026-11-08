"""Conditional Gibbs-energy bounds without interpolation assumptions.

The exact piecewise-linear interval is conditional on a true global
Lipschitz constant and valid energy error bars over the FULL volume domain.
No source is inferred or verified from a provenance string.
Energies: eV/formula; volumes: Angstrom^3/formula; P: GPa.
"""
from dataclasses import dataclass
from fractions import Fraction
import math

# 1 GPa Å³ = 10⁻²¹ J. The elementary charge and Avogadro constant are exact.
GPA_A3_TO_EV_EXACT = Fraction(10**7, 1602176634)
GPA_A3_TO_EV = float(GPA_A3_TO_EV_EXACT)
EV_FORMULA_TO_KJ_MOL = Fraction("96.4853321233100184")
Numeric = int | float | Fraction

def finite(x):
    return type(x) in (int, Fraction) or isinstance(x, float) and math.isfinite(x)

@dataclass(frozen=True)
class Sample:
    volume: Numeric
    energy: Numeric
    error: Numeric

@dataclass(frozen=True)
class Phase:
    name: str
    composition: tuple
    samples: tuple
    lipschitz_ev_per_a3: Numeric
    reference: str
    method: str
    source: str
    error_source: str
    lipschitz_source: str
    temperature_k: float
    evidence: str = "illustrative"

    def validate(self):
        for field in ("name", "reference", "method", "source", "error_source", "lipschitz_source"):
            if not isinstance(getattr(self, field), str) or not getattr(self, field).strip():
                raise ValueError(f"{field} is required")
        if self.evidence not in ("illustrative", "calculated", "reconstructed", "measured"):
            raise ValueError("invalid evidence class")
        if not finite(self.temperature_k) or self.temperature_k < 0:
            raise ValueError("temperature must be nonnegative")
        L = self.lipschitz_ev_per_a3
        if not finite(L) or L < 0:
            raise ValueError("invalid global Lipschitz bound")
        if not self.composition or len(set(k for k, _ in self.composition)) != len(self.composition):
            raise ValueError("unique formula composition required")
        if any(not isinstance(k, str) or not k.strip() or type(v) is not int or v <= 0
               for k, v in self.composition):
            raise ValueError("invalid stoichiometry")
        if len(self.samples) < 2:
            raise ValueError("two or more volume samples needed")
        if any(not isinstance(s, Sample) or not all(finite(x) for x in (s.volume, s.energy, s.error))
               or s.volume <= 0 or s.error < 0 for s in self.samples):
            raise ValueError("invalid samples")
        if any(a.volume >= b.volume for a, b in zip(self.samples, self.samples[1:])):
            raise ValueError("volumes must strictly increase")
        # Each shifted endpoint is checked against all prior samples via
        # running extrema. This is equivalent to pairwise consistency.
        L = Fraction(L)
        prior_lower, prior_upper = -math.inf, math.inf
        for v, e, error in _rational_samples(self):
            if prior_lower > e + error + L * v or prior_upper < e - error - L * v:
                raise ValueError("inconsistent Lipschitz and sample intervals")
            prior_lower = max(prior_lower, e - error + L * v)
            prior_upper = min(prior_upper, e + error - L * v)


def _rational_samples(phase):
    return [(Fraction(s.volume), Fraction(s.energy), Fraction(s.error))
            for s in phase.samples]


def envelope_minimum(phase, p, lower):
    """Return the exact (energy, volume) minimum for a validated phase."""
    samples = _rational_samples(phase)
    L, p = Fraction(phase.lipschitz_ev_per_a3), Fraction(p)
    first, last = samples[0][0], samples[-1][0]
    if not lower:
        # A minimum over volumes and upper cones can be taken in either order.
        candidates = []
        for v, e, error in samples:
            at = first if p >= L else last if p < -L else v
            candidates.append((e + error + L * abs(at - v) + p * at, at))
        return min(candidates)

    # Between knots the lower envelope is the maximum of two lines. All
    # cones from the left share slope -L and those from the right share +L.
    suffix = [None] * len(samples)
    right_intercept = -math.inf
    for i in range(len(samples) - 1, -1, -1):
        v, e, error = samples[i]
        right_intercept = max(right_intercept, e - error - L * v)
        suffix[i] = right_intercept

    left_intercept = -math.inf
    best = (math.inf, first)
    for i, (left, e, error) in enumerate(samples[:-1]):
        right = samples[i + 1][0]
        left_intercept = max(left_intercept, e - error + L * left)
        right_intercept = suffix[i + 1]
        candidates = [left, right]
        if L:
            crossing = (left_intercept - right_intercept) / (2 * L)
            if left < crossing < right:
                candidates.append(crossing)
        for v in candidates:
            value = max(left_intercept - L * v, right_intercept + L * v) + p * v
            best = min(best, (value, v))
    return best


def _outward_float(value, *, upper):
    rounded = float(value)
    if (upper and rounded < value) or (not upper and rounded > value):
        rounded = math.nextafter(rounded, math.inf if upper else -math.inf)
    if not math.isfinite(rounded):
        raise OverflowError("energy interval exceeds finite float output range")
    return rounded


def minimum_bounds(phase, pressure_gpa, *, allow_illustrative=False):
    """Sharp min(F+PV) interval under independently justified assumptions."""
    phase.validate()
    if phase.evidence == "illustrative" and not allow_illustrative:
        raise ValueError("illustrative input requires explicit override")
    if not finite(pressure_gpa) or pressure_gpa < 0:
        raise ValueError("pressure must be finite and nonnegative")
    p = Fraction(pressure_gpa) * GPA_A3_TO_EV_EXACT
    low, lower_argmin = envelope_minimum(phase, p, True)
    high, upper_argmin = envelope_minimum(phase, p, False)
    return {"phase": phase.name, "pressure_gpa": float(pressure_gpa),
            "temperature_k": float(phase.temperature_k),
            "interval_ev_formula": [_outward_float(low, upper=False),
                                    _outward_float(high, upper=True)],
            "volume_domain_a3_formula": [float(phase.samples[0].volume),
                                          float(phase.samples[-1].volume)],
            "lower_envelope_argmin_a3_formula": float(lower_argmin),
            "upper_envelope_argmin_a3_formula": float(upper_argmin),
            "method": "exact Lipschitz-envelope minimization, no spline/EOS extrapolation",
            "conditional": True, "evidence": phase.evidence,
            "source": phase.source, "error_source": phase.error_source,
            "lipschitz_source": phase.lipschitz_source}

def compare(phases, pressure_gpa, *, allow_illustrative=False):
    """Claim only a conditional winner among supplied comparable polymorphs."""
    if len(phases) < 2:
        raise ValueError("two or more phases required")
    for p in phases:
        p.validate()
    if len({p.name for p in phases}) != len(phases):
        raise ValueError("duplicate phase")
    a = phases[0]
    for b in phases[1:]:
        if (dict(b.composition) != dict(a.composition)
            or b.reference != a.reference or b.method != a.method
            or b.temperature_k != a.temperature_k):
            raise ValueError("incompatible composition, reference, method or temperature")
    results = [minimum_bounds(p, pressure_gpa, allow_illustrative=allow_illustrative) for p in phases]
    winners = [r for r in results if all(r["interval_ev_formula"][1] < s["interval_ev_formula"][0]
               for s in results if r is not s)]
    winner = winners[0]["phase"] if len(winners) == 1 else None
    return {"status": "conditional_bounded_winner" if winner else "abstain_overlap",
            "winner": winner, "results": results,
            "scope": "supplied polymorphs and sampled volume domains only; exterior minima, source error and global Lipschitz assumptions unverified"}
