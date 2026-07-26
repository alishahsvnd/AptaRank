"""The one and only ViennaRNA entry point (spec §4.5, §4.7, §7.4).

ViennaRNA's API is order-dependent: `fc.pf()` must be called before
`ensemble_defect()` or `bpp()`, or those return garbage rather than an error.
Every fold in this project therefore goes through `fold()` below; the raw
`RNA.fold_compound` API is not called from anywhere else.

Two further non-obvious requirements, both verified against ViennaRNA 2.7.2:

* `md.uniq_ML = 1` must be set *at fold_compound construction time*, otherwise
  `pbacktrack()` silently returns an empty tuple instead of raising.
* `exp_params_rescale(mfe)` before `pf()` keeps the partition function
  numerically stable for longer sequences.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path

warnings.filterwarnings("ignore", category=DeprecationWarning)

import RNA  # noqa: E402  (import after warning filter is deliberate)

from ..errors import FoldingError


@dataclass(frozen=True)
class FoldResult:
    """Everything a single fold tells us about one sequence."""

    sequence: str
    dot_bracket: str
    mfe: float
    mfe_norm: float
    ensemble_free_energy: float
    ensemble_defect: float
    positional_entropy_mean: float
    positional_entropy: tuple[float, ...]  # per position, 5'->3', length n

    @property
    def length(self) -> int:
        return len(self.sequence)

    def feature_dict(self) -> dict[str, float | str]:
        return {
            "dot_bracket": self.dot_bracket,
            "mfe": self.mfe,
            "mfe_norm": self.mfe_norm,
            "ensemble_free_energy": self.ensemble_free_energy,
            "ensemble_defect": self.ensemble_defect,
            "positional_entropy_mean": self.positional_entropy_mean,
        }


def fold(sequence: str) -> FoldResult:
    """MFE structure plus the three ensemble statistics of spec §4.5."""
    return _finish(*_prepared_compound(sequence), sequence=sequence)


def fold_and_sample(
    sequence: str, n_samples: int, seed: int | None = None
) -> tuple[FoldResult, tuple[str, ...]]:
    """Fold, then stochastically sample `n_samples` structures (spec §4.7).

    Sampling is seeded per call so results do not depend on the order in which
    a process pool happens to schedule candidates.
    """
    prepared = _prepared_compound(sequence)
    fc = prepared[0]
    result = _finish(*prepared, sequence=sequence)

    if n_samples <= 0:
        return result, ()

    if seed is not None:
        RNA.init_rand(int(seed))
    samples = tuple(fc.pbacktrack(int(n_samples)))
    for db in samples:
        _validate_dot_bracket(db, len(sequence))

    if not samples:
        # The uniq_ML failure mode: pbacktrack returns an empty tuple with no
        # error when md.uniq_ML was not set at fold_compound construction.
        raise FoldingError(
            f"pbacktrack returned no structures for a sequence of length "
            f"{len(sequence)}. This is the md.uniq_ML failure mode — sampling "
            f"is silently empty without it."
        )
    # Occasionally ViennaRNA returns one or two fewer structures than asked for
    # (a redundancy rejection inside the sampler). That is harmless: the loop
    # statistics are computed over however many structures came back, and the
    # count is recorded in the artifact as n_ensemble_samples.
    return result, samples


def plot_svg(sequence: str, dot_bracket: str, path: str | Path) -> Path:
    """Render a secondary-structure diagram (spec §7.4)."""
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    ok = RNA.svg_rna_plot(sequence, dot_bracket, str(out))
    if not ok or not out.exists():
        raise FoldingError(f"ViennaRNA failed to render a structure diagram to {out}")
    return out


#: Model settings that change every number this module produces. Recorded in
#: the corpus cache key and the run artifact, so a cache can never survive a
#: change to the folding model.
MODEL_SETTING_NAMES = (
    "temperature", "dangles", "noLP", "noGU", "no_closingGU",
    "special_hp", "logML", "uniq_ML", "max_bp_span", "min_loop_size",
)


def model_settings() -> dict[str, float | int]:
    """The effective ViennaRNA model details used for every fold here."""
    md = RNA.md()
    md.uniq_ML = 1
    out: dict[str, float | int] = {}
    for name in MODEL_SETTING_NAMES:
        value = getattr(md, name, None)
        if value is not None:
            out[name] = value if isinstance(value, (int, float)) else str(value)
    return out


def _validate_dot_bracket(dot_bracket: str, expected_length: int) -> None:
    """Reject a structure that is truncated, unbalanced, or not dot-bracket.

    Silently accepting a malformed sample would bias the §4.7 loop-size
    distribution, which is Tier 2's only input from Tier 1.
    """
    if len(dot_bracket) != expected_length:
        raise FoldingError(
            f"sampled structure has length {len(dot_bracket)}, expected {expected_length}"
        )
    depth = 0
    for char in dot_bracket:
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
            if depth < 0:
                raise FoldingError(f"unbalanced sampled structure: {dot_bracket!r}")
        elif char != ".":
            raise FoldingError(f"unexpected character {char!r} in {dot_bracket!r}")
    if depth != 0:
        raise FoldingError(f"unbalanced sampled structure: {dot_bracket!r}")


def _prepared_compound(sequence: str):
    if not sequence:
        raise FoldingError("cannot fold an empty sequence")
    md = RNA.md()
    md.uniq_ML = 1  # REQUIRED for pbacktrack (spec §4.5)
    fc = RNA.fold_compound(sequence, md)
    mfe_struct, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _pf_struct, efe = fc.pf()  # must precede ensemble_defect() and bpp()
    return fc, mfe_struct, float(mfe), float(efe)


def _finish(fc, mfe_struct: str, mfe: float, efe: float, *, sequence: str) -> FoldResult:
    n = len(sequence)
    entropy = fc.positional_entropy()  # length n+1, index 0 unused
    per_position = tuple(float(x) for x in entropy[1 : n + 1])
    return FoldResult(
        sequence=sequence,
        dot_bracket=mfe_struct,
        mfe=mfe,
        mfe_norm=mfe / n,
        ensemble_free_energy=efe,
        ensemble_defect=float(fc.ensemble_defect(mfe_struct)),
        positional_entropy_mean=sum(per_position) / n if n else 0.0,
        positional_entropy=per_position,
    )
