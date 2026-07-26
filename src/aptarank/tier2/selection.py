"""Which cavity is the functional one? (spec §5.4)

fpocket ranks cavities by its own score, which is not guaranteed to identify
the functionally relevant one, so this is deliberately not fully automated:
literature-confirmed active-site residues select the pocket, and the automatic
fallback is recorded in the bundle so the UI can caveat it.

Zero overlap between a supplied residue list and every pocket is treated as a
build failure by default. It is far more likely to mean a numbering, chain or
parsing mismatch than a genuine finding, and silently falling back to the
top-scoring pocket would hide that.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Iterable, Sequence

from ..errors import TargetError
from .fpocket import Pocket, Residue

TIE_BREAK_ORDER = (
    "overlap_count_desc",
    "fpocket_score_desc",
    "druggability_score_desc",
    "pocket_index_asc",
)


@dataclass(frozen=True)
class ResidueSpec:
    """A requested active-site residue. Name is validation, not identity."""

    chain_id: str
    residue_number: int
    insertion_code: str = ""
    residue_name: str | None = None

    def key(self) -> tuple[str, int, str]:
        return (self.chain_id.strip(), self.residue_number, self.insertion_code.strip())

    def to_dict(self) -> dict[str, Any]:
        return {
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "residue_name": self.residue_name,
        }


def parse_residue_specs(raw: Iterable[Any], default_chain: str) -> list[ResidueSpec]:
    """Accept plain residue numbers (spec §9) or explicit chain/icode dicts."""
    specs: list[ResidueSpec] = []
    for item in raw or []:
        if isinstance(item, int):
            specs.append(ResidueSpec(default_chain, item))
        elif isinstance(item, str) and item.strip().lstrip("-").isdigit():
            specs.append(ResidueSpec(default_chain, int(item.strip())))
        elif isinstance(item, dict):
            specs.append(
                ResidueSpec(
                    chain_id=str(item.get("chain_id", default_chain)),
                    residue_number=int(item["residue_number"]),
                    insertion_code=str(item.get("insertion_code", "")),
                    residue_name=item.get("residue_name"),
                )
            )
        else:
            raise TargetError(
                f"cannot interpret active-site residue {item!r}; use an integer "
                f"residue number or a mapping with residue_number"
            )
    return specs


def select_pocket(
    pockets: Sequence[Pocket],
    requested: Sequence[ResidueSpec] = (),
    structure_residues: Sequence[Residue] = (),
    allow_zero_overlap_fallback: bool = False,
) -> dict[str, Any]:
    """Choose the pocket and record the full evidence for that choice."""
    if not pockets:
        raise TargetError("cannot select a pocket: fpocket detected none")

    present = {r.key() for r in structure_residues}
    missing = [spec for spec in requested if structure_residues and spec.key() not in present]
    if missing:
        # A residue that is not in the cleaned structure at all cannot overlap
        # anything, and almost always means the wrong chain or numbering.
        raise TargetError(
            f"{len(missing)} configured active-site residue(s) are absent from the "
            f"prepared structure: "
            f"{[f'{s.chain_id}{s.residue_number}{s.insertion_code}' for s in missing[:6]]}. "
            f"Check the chain selector and the residue numbering scheme."
        )

    wanted = {spec.key() for spec in requested}
    evidence = []
    for pocket in pockets:
        lining = {r.key(): r for r in pocket.lining_residues}
        overlapping = [lining[key] for key in wanted & set(lining)]
        evidence.append(
            {
                "pocket_index": pocket.index,
                "overlap_count": len(overlapping),
                "overlapping_residues": [r.to_dict() for r in overlapping],
                "fpocket_score": pocket.score,
                "druggability_score": pocket.druggability,
                "selected": False,
            }
        )

    warnings: list[str] = []
    total_overlap = sum(item["overlap_count"] for item in evidence)

    if requested and total_overlap == 0:
        if not allow_zero_overlap_fallback:
            raise TargetError(
                f"none of the {len(requested)} configured active-site residues line "
                f"any of the {len(pockets)} detected pockets. This usually means a "
                f"numbering, chain or preparation mismatch rather than a real "
                f"result. Set tier2.target.allow_zero_overlap_fallback to accept "
                f"automatic selection instead."
            )
        method = "active_site_zero_overlap_fallback"
        warnings.append(
            "active-site residues were supplied but overlapped no pocket; "
            "selection fell back to the highest fpocket score and must NOT be "
            "described as active-site selected"
        )
    elif requested:
        method = "active_site_overlap"
    else:
        method = "automatic_fpocket_score"
        warnings.append(
            "no active-site residues configured; the pocket was chosen by fpocket "
            "score alone, which is not guaranteed to be the functional cavity"
        )

    ordered = sorted(
        evidence,
        key=lambda e: (
            -e["overlap_count"],
            -e["fpocket_score"],
            -(e["druggability_score"] if e["druggability_score"] is not None else float("-inf")),
            e["pocket_index"],
        ),
    )
    chosen = ordered[0]
    chosen["selected"] = True

    return {
        "status": "selected",
        "method": method,
        "selected_pocket_index": chosen["pocket_index"],
        "active_site": {
            "requested": bool(requested),
            "allow_zero_overlap_fallback": allow_zero_overlap_fallback,
            "requested_residues": [spec.to_dict() for spec in requested],
            "n_requested": len(requested),
            "total_overlap": total_overlap,
        },
        "pocket_evidence": evidence,
        "tie_break_order": list(TIE_BREAK_ORDER),
        "warnings": warnings,
    }
