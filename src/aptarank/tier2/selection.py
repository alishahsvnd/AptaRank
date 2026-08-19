"""Which cavity is the functional one? (spec §5.4)

fpocket ranks cavities by its own score, which is not guaranteed to identify
the functionally relevant one, so this is deliberately not fully automated:
literature-confirmed binding-site residues select the pocket, and the automatic
fallback is recorded in the bundle so the UI can caveat it.

"Binding site", not "active site": an active site is an enzyme's catalytic
machinery, and most aptamer targets are not enzymes. The old name quietly
implied the tool only worked on enzymes.

In pocket mode, zero overlap between a supplied residue list and every cavity is
a build failure. It is far more likely to mean a numbering, chain or parsing
mismatch than a genuine finding, and silently falling back to the top-scoring
cavity would hide that. In surface mode the patch is the evidence and a cavity
is only a cross-reference, so the same situation is expected and merely noted.
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
    """A requested binding-site residue. Name is validation, not identity."""

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
    require_overlap: bool = True,
) -> dict[str, Any]:
    """Choose the pocket and record the full evidence for that choice."""
    if not pockets:
        if require_overlap:
            raise TargetError("cannot select a pocket: fpocket detected none")
        return {
            "status": "not_applicable",
            "method": "no_cavity_detected",
            "selected_pocket_index": None,
            "target_site": {
                "requested": bool(requested),
                "requested_residues": [spec.to_dict() for spec in requested],
                "n_requested": len(requested),
                "total_overlap": 0,
            },
            "pocket_evidence": [],
            "tie_break_order": list(TIE_BREAK_ORDER),
            "warnings": [
                "no cavity was detected on this target; in surface mode the "
                "measured patch is the evidence and no cavity is needed"
            ],
        }

    present = {r.key() for r in structure_residues}
    missing = [spec for spec in requested if structure_residues and spec.key() not in present]
    if missing:
        # A residue that is not in the cleaned structure at all cannot overlap
        # anything, and almost always means the wrong chain or numbering.
        raise TargetError(
            f"{len(missing)} configured binding-site residue(s) are absent from the "
            f"prepared structure: "
            f"{[f'{s.chain_id}{s.residue_number}{s.insertion_code}' for s in missing[:6]]}. "
            f"Check the chain selector and the residue numbering scheme."
        )

    wanted = {spec.key() for spec in requested}
    evidence = []
    for pocket in pockets:
        lining = {r.key(): r for r in pocket.lining_residues}
        # sorted(), because iterating a set intersection gives an order that
        # varies between processes (string hashing is randomised per run). The
        # residues would be the same but the recorded evidence would differ, and
        # anything hashing this bundle would call two identical builds different.
        overlapping = [lining[key] for key in sorted(wanted & set(lining))]
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
        if require_overlap and not allow_zero_overlap_fallback:
            raise TargetError(
                f"none of the {len(requested)} configured binding-site residues line "
                f"any of the {len(pockets)} detected cavities. This usually means a "
                f"numbering, chain or preparation mismatch rather than a real "
                f"result. Set tier2.target.allow_zero_overlap_fallback to accept "
                f"automatic selection instead."
            )
        method = "target_site_zero_overlap_fallback"
        warnings.append(
            "binding-site residues were supplied but overlapped no cavity; "
            "selection fell back to the highest fpocket score and must NOT be "
            "described as binding-site selected"
            + ("" if require_overlap else
               " (expected in surface mode, where the patch is the evidence)")
        )
    elif requested:
        method = "target_site_overlap"
    else:
        method = "automatic_fpocket_score"
        warnings.append(
            "no binding-site residues configured; the cavity was chosen by fpocket "
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
        "target_site": {
            "requested": bool(requested),
            "allow_zero_overlap_fallback": allow_zero_overlap_fallback,
            "require_overlap": require_overlap,
            "requested_residues": [spec.to_dict() for spec in requested],
            "n_requested": len(requested),
            "total_overlap": total_overlap,
        },
        "pocket_evidence": evidence,
        "tie_break_order": list(TIE_BREAK_ORDER),
        "warnings": warnings,
    }
