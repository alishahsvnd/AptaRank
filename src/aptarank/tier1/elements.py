"""Structural element features from dot-bracket, via forgi (spec §4.6).

Dot-bracket is a flat encoding. forgi turns it into named elements — stems
(`s`), hairpin loops (`h`), interior loops (`i`), multiloops (`m`) and 5'/3'
tails (`f`/`t`). The loops matter most: they are the unpaired regions that
physically stick out and can contact a target, which is what the domain
experts asked us to surface.

Install gotcha: plain `import forgi` pulls in an optional 3D submodule that
emits a NumPy 2.x compatibility traceback. We import the graph module directly
and never touch the 3D code.
"""

from __future__ import annotations

import contextlib
import io
import warnings
from dataclasses import dataclass, field

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="forgi")

# Importing any forgi submodule executes forgi/__init__.py, which pulls in the
# optional 3D module. That module is compiled against NumPy 1.x, so on NumPy 2
# it prints a multi-line traceback to stderr before forgi swallows the error
# and carries on. The import succeeds and nothing we use is affected — but the
# noise lands in job logs, which the dashboard shows to users when a run fails,
# where an unexplained traceback reads as a crash. Silence it here rather than
# pinning NumPy back for a module this project never touches.
with contextlib.redirect_stderr(io.StringIO()):
    from forgi.graph.bulge_graph import BulgeGraph  # noqa: E402

LOOP_CODES = ("h", "i", "m")


@dataclass(frozen=True)
class ElementFeatures:
    """Named structural elements of one predicted secondary structure."""

    n_hairpins: int
    n_interior: int
    n_multiloop: int
    n_stems: int
    stem_fraction: float
    longest_stem_bp: int
    max_loop_nt: int
    total_unpaired: int
    loop_sizes: tuple[int, ...] = field(default=())
    element_string: str = ""

    def feature_dict(self) -> dict[str, float | int]:
        return {
            "n_hairpins": self.n_hairpins,
            "n_interior": self.n_interior,
            "n_multiloop": self.n_multiloop,
            "n_stems": self.n_stems,
            "stem_fraction": self.stem_fraction,
            "longest_stem_bp": self.longest_stem_bp,
            "max_loop_nt": self.max_loop_nt,
            "total_unpaired": self.total_unpaired,
        }


def parse_elements(dot_bracket: str, sequence: str | None = None) -> ElementFeatures:
    """Extract element features from a dot-bracket string.

    An entirely unpaired structure has no elements forgi can name as loops; it
    is reported with zero stems and zero loops rather than raising, because a
    generated candidate that fails to fold at all is a legitimate (very poor)
    candidate and must still be scored.
    """
    n = len(dot_bracket)
    total_unpaired = dot_bracket.count(".")

    if "(" not in dot_bracket:
        return ElementFeatures(
            n_hairpins=0,
            n_interior=0,
            n_multiloop=0,
            n_stems=0,
            stem_fraction=0.0,
            longest_stem_bp=0,
            max_loop_nt=total_unpaired,
            total_unpaired=total_unpaired,
            loop_sizes=(total_unpaired,) if total_unpaired else (),
            element_string="f" * n,
        )

    bg = BulgeGraph.from_dotbracket(dot_bracket, seq=sequence)
    by_code: dict[str, list[str]] = {}
    for name in bg.defines:
        by_code.setdefault(name[0], []).append(name)

    stems = by_code.get("s", [])
    stem_bp = [int(bg.stem_length(e)) for e in stems]
    stem_nt = sum(2 * bp for bp in stem_bp)

    loop_sizes = [
        int(bg.element_length(e))
        for code in LOOP_CODES
        for e in by_code.get(code, [])
    ]

    return ElementFeatures(
        n_hairpins=len(by_code.get("h", [])),
        n_interior=len(by_code.get("i", [])),
        n_multiloop=len(by_code.get("m", [])),
        n_stems=len(stems),
        stem_fraction=stem_nt / n if n else 0.0,
        longest_stem_bp=max(stem_bp, default=0),
        max_loop_nt=max(loop_sizes, default=0),
        total_unpaired=total_unpaired,
        loop_sizes=tuple(loop_sizes),
        element_string=bg.to_element_string(),
    )


def max_loop_nt(dot_bracket: str) -> int:
    """Cheap path used when sampling hundreds of structures per candidate."""
    return parse_elements(dot_bracket).max_loop_nt
