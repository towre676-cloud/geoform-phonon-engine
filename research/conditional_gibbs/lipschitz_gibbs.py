"""Conditional Gibbs-energy bounds without interpolation assumptions.

The exact piecewise-linear interval is conditional on a true global
Lipschitz constant and valid energy error bars over the FULL volume domain.
No source is inferred or verified from a provenance string.
Energies: eV/formula; volumes: Angstrom^3/formula; P: GPa.
"""
from dataclasses import dataclass
import math

GPA_A3_TO_EV = 0.006241509074460763

def finite(x):
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)

@dataclass(frozen=True)
class Sample:
    volume: float
    energy: float
    error: float

@dataclass(frozen=True)
class Phase:
    name: str
    composition: tuple
    samples: tuple
    lipschitz_ev_per_a3: float
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
        # Necessary and sufficient for existence of an L-Lipschitz graph
        # through the error intervals; does NOT prove the true graph is one.
        for i, a in enumerate(self.samples):
            for b in self.samples[i + 1:]:
                if abs(a.energy - b.energy) > L * (b.volume - a.volume) + a.error + b.error + 1e-12:
                    raise ValueError("inconsistent Lipschitz and sample intervals")

def envelope_minimum(phase, p, lower):
    """Exact minimum of max lower cones or min upper cones plus p*V."""
    L = phase.lipschitz_ev_per_a3
    best = (math.inf, None)
    for a, b in zip(phase.samples, phase.samples[1:]):
        left, right = a.volume, b.volume
        midpoint = (left + right) / 2
        lines = []
        for s in phase.samples:
            side = 1 if midpoint >= s.volume else -1
            if lower:
                slope, intercept = -L * side, s.energy - s.error + L * side * s.volume
            else:
                slope, intercept = L * side, s.energy + s.error - L * side * s.volume
            lines.append((slope, intercept))
        candidates = [left, right]
        for i, (m, c) in enumerate(lines):
            for n, d in lines[i + 1:]:
                if abs(m - n) < 1e-15:
                    continue
                root = (d - c) / (m - n)
                if left < root < right:
                    candidates.append(root)
        for v in candidates:
            values = [m * v + c for m, c in lines]
            energy = (max(values) if lower else min(values)) + p * v
            if energy < best[0]:
                best = (energy, v)
    return best

def minimum_bounds(phase, pressure_gpa, *, allow_illustrative=False):
    """Sharp min(F+PV) interval under independently justified assumptions."""
    phase.validate()
    if phase.evidence == "illustrative" and not allow_illustrative:
        raise ValueError("illustrative input requires explicit override")
    if not finite(pressure_gpa) or pressure_gpa < 0:
        raise ValueError("pressure must be finite and nonnegative")
    p = pressure_gpa * GPA_A3_TO_EV
    low, lower_argmin = envelope_minimum(phase, p, True)
    high, upper_argmin = envelope_minimum(phase, p, False)
    if low > high + 1e-9:
        raise ArithmeticError("impossible bounds")
    return {"phase": phase.name, "pressure_gpa": pressure_gpa,
            "temperature_k": phase.temperature_k,
            "interval_ev_formula": [low, high],
            "lower_envelope_argmin_a3_formula": lower_argmin,
            "upper_envelope_argmin_a3_formula": upper_argmin,
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
            "scope": "supplied polymorphs only; source error and global Lipschitz assumptions unverified"}
