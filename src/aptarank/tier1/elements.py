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
import math
import warnings
from collections import deque
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
#: Elements whose two strands bridge the same gap between the same pair of stems.
CONNECTING_CODES = ("i", "m")
#: A Gaussian segment of n monomers with bond length b has Rg = b*sqrt(n/6).
GAUSSIAN_RG_DIVISOR = 6.0
#: A rigid rod of length L has Rg = L/sqrt(12).
ROD_RG_DIVISOR = 12.0


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
    #: Overall size of the folded molecule, in Å (see `radius_of_gyration`).
    radius_of_gyration_A: float = 0.0

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
            "radius_of_gyration_A": self.radius_of_gyration_A,
        }


def parse_elements(
    dot_bracket: str,
    sequence: str | None = None,
    a_per_bp_helix: float = 2.8,
    a_per_nt_ss: float = 6.0,
) -> ElementFeatures:
    """Extract element features from a dot-bracket string.

    An entirely unpaired structure has no elements forgi can name as loops; it
    is reported with zero stems and zero loops rather than raising, because a
    generated candidate that fails to fold at all is a legitimate (very poor)
    candidate and must still be scored.

    The two lengths are nucleic-acid geometry, not scoring knobs, and they reach
    here from `tier2.geometry` in the config so that one file remains the single
    source of truth for them. The defaults match it and exist so that this
    function stays usable on its own.
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
            radius_of_gyration_A=free_coil_radius_of_gyration(n, a_per_nt_ss),
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
        radius_of_gyration_A=radius_of_gyration(bg, a_per_bp_helix, a_per_nt_ss),
    )


def radius_of_gyration(
    bulge_graph, a_per_bp_helix: float, a_per_nt_ss: float
) -> float:
    """Overall size of the folded molecule, in Å, from its element graph.

    The textbook mass-weighted definition, evaluated on the secondary structure
    rather than on coordinates we do not have:

        Rg² = Σ (mₑ/M)·rgₑ²  +  (1/2M²)·ΣΣ mᵢmⱼ·d(i,j)²

    Each element carries its nucleotides as mass, distances between elements are
    measured along the structure graph, and each element's own internal spread
    is added by the parallel-axis theorem. Stems are rigid A-form rods; unpaired
    regions are Gaussian segments, the same √L convention Tier 2 already uses for
    loop reach.

    What this buys over "size ∝ length": two 100-nt candidates now differ. One
    long helix is an extended rod and one four-way junction is a compact ball,
    and a surface footprint proxy that called them identical was measuring the
    sequence rather than the fold.

    Known coarseness, stated rather than hidden: the secondary structure fixes
    no angles, so a multiloop's arms are assumed to splay rather than stack
    coaxially. Real junctions stack, which makes real molecules somewhat more
    extended than this estimate.
    """
    info = {
        name: _element_geometry(bulge_graph, name, a_per_bp_helix, a_per_nt_ss)
        for name in bulge_graph.defines
    }
    if not info:
        return 0.0
    edges = {name: set(bulge_graph.edges.get(name, ())) for name in info}
    total_mass = sum(mass for mass, _extent, _rg in info.values())
    if total_mass <= 0:
        return 0.0

    spread = 0.0
    for start in info:
        for other, distance in _graph_distances(start, info, edges).items():
            spread += info[start][0] * info[other][0] * distance * distance
    internal = sum(mass * rg * rg for mass, _extent, rg in info.values())

    return math.sqrt(spread / (2.0 * total_mass**2) + internal / total_mass)


def _element_geometry(
    bulge_graph, name: str, a_per_bp_helix: float, a_per_nt_ss: float
) -> tuple[int, float, float]:
    """(mass in nt, end-to-end extent in Å, internal Rg in Å) for one element."""
    if name.startswith("s"):
        base_pairs = int(bulge_graph.stem_length(name))
        length = a_per_bp_helix * base_pairs
        return 2 * base_pairs, length, length / math.sqrt(ROD_RG_DIVISOR)

    strands = _strand_lengths(bulge_graph, name)
    mass = sum(strands)
    # A connecting element's two strands bridge the same gap between the same
    # two stems, so the shorter one limits how far apart they can sit. A
    # terminal element (hairpin, tail) has only one strand to extend along.
    reach = (
        min(strands) if name[0] in CONNECTING_CODES and len(strands) > 1
        else max(strands, default=0)
    )
    return (
        mass,
        a_per_nt_ss * math.sqrt(max(reach, 0)),
        a_per_nt_ss * math.sqrt(max(mass, 0) / GAUSSIAN_RG_DIVISOR),
    )


def _strand_lengths(bulge_graph, name: str) -> list[int]:
    """Nucleotides per strand, from forgi's [start, end, start, end] ranges."""
    define = bulge_graph.defines[name]
    return [define[i + 1] - define[i] + 1 for i in range(0, len(define) - 1, 2)] or [0]


def _graph_distances(
    start: str, info: dict[str, tuple[int, float, float]], edges: dict[str, set[str]]
) -> dict[str, float]:
    """Centre-to-centre distance from `start` to every reachable element."""
    out = {start: 0.0}
    queue = deque([start])
    while queue:
        node = queue.popleft()
        for nxt in edges[node] - out.keys():
            out[nxt] = out[node] + info[node][1] / 2.0 + info[nxt][1] / 2.0
            queue.append(nxt)
    return out


def free_coil_radius_of_gyration(n_nucleotides: int, a_per_nt_ss: float) -> float:
    """Rg of an unstructured strand: a Gaussian coil, with nothing to fold on."""
    return a_per_nt_ss * math.sqrt(max(n_nucleotides, 0) / GAUSSIAN_RG_DIVISOR)


def max_loop_nt(dot_bracket: str) -> int:
    """Cheap path used when sampling hundreds of structures per candidate."""
    return parse_elements(dot_bracket).max_loop_nt
