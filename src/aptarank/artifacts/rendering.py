"""Secondary-structure diagrams (spec §7.4).

Rendered for the top N candidates only and embedded as SVG text so the run
artifact stays self-contained.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Iterable, Mapping

from ..errors import FoldingError
from ..tier1 import folding


def render_diagrams(
    rows: Iterable[Mapping[str, str]],
    n: int,
    out_dir: str | Path | None = None,
) -> dict[str, str]:
    """Render up to `n` diagrams; returns {candidate_id: svg_text}.

    A candidate whose diagram fails to render is skipped rather than failing
    the run — a missing picture must never cost a whole batch of scores.
    """
    diagrams: dict[str, str] = {}
    if n <= 0:
        return diagrams

    target_dir = Path(out_dir) if out_dir else Path(tempfile.mkdtemp(prefix="aptarank_svg_"))
    target_dir.mkdir(parents=True, exist_ok=True)

    for row in list(rows)[:n]:
        cid = row["candidate_id"]
        path = target_dir / f"{cid}.svg"
        try:
            folding.plot_svg(row["sequence"], row["dot_bracket"], path)
            diagrams[cid] = path.read_text(encoding="utf-8")
        except (FoldingError, OSError):
            continue
    return diagrams
