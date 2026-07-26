"""Fetch and prepare a target structure (spec §5.3).

The hetero-atom policy is the part that matters scientifically. The spec says
"keep non-water HETATM by default" because the demo targets have a metal ion
sitting in the functional cavity, and removing it changes the cavity geometry.
But "keep everything" is too broad: crystallisation buffers, cryoprotectants
and co-crystallised inhibitors also arrive as HETATM and can *block* the very
cavity we are trying to measure.

So the policy is explicit per target — `retain_hetero_resnames` /
`remove_hetero_resnames` — and whatever was actually kept or dropped is
recorded residue by residue in the bundle.
"""

from __future__ import annotations

import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Sequence

from ..errors import TargetError
from ..provenance import sha256_file
from .fpocket import Residue

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
WATER_RESNAMES = ("HOH", "DOD", "WAT")
MIN_ATOMS = 100

AMINO_ACIDS = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
    "MSE SEC PYL".split()
)


@dataclass
class PreparedTarget:
    """A single cleaned chain, ready for fpocket, plus what was done to it."""

    pdb_id: str
    chain_id: str
    model_index: int
    path: Path
    source: dict[str, Any]
    applied: dict[str, Any]
    residues: list[Residue] = field(default_factory=list)
    name: str | None = None


def fetch_structure(pdb_id: str, cache_dir: str | Path) -> dict[str, Any]:
    """Download `{pdb_id}.pdb` from RCSB, cached by id."""
    ident = pdb_id.strip().upper()
    if len(ident) != 4 or not ident.isalnum():
        raise TargetError(f"expected a 4-character PDB ID, got {pdb_id!r}")

    directory = Path(cache_dir)
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{ident}.pdb"
    url = RCSB_URL.format(pdb_id=ident)

    if not path.exists():
        try:
            with urllib.request.urlopen(url, timeout=60) as response:
                payload = response.read()
        except (urllib.error.URLError, OSError) as exc:
            raise TargetError(f"could not download {url}: {exc}") from exc
        if len(payload) < 1000:
            raise TargetError(f"{url} returned {len(payload)} bytes; not a structure")
        path.write_bytes(payload)

    return {
        "url": url,
        "format": "pdb",
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def prepare(
    pdb_id: str,
    cache_dir: str | Path,
    work_dir: str | Path,
    chain_id: str | None = None,
    model_index: int = 0,
    retain_hetero_resnames: Sequence[str] = (),
    remove_hetero_resnames: Sequence[str] = WATER_RESNAMES,
    hetero_default: str = "retain",
) -> PreparedTarget:
    """Select one chain, apply the hetero policy, write a cleaned PDB."""
    from Bio.PDB import PDBIO, PDBParser, Select  # imported here: heavy

    source = fetch_structure(pdb_id, cache_dir)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id.upper(), source["path"])

    try:
        model = structure[model_index]
    except KeyError as exc:
        raise TargetError(f"{pdb_id}: no model {model_index}") from exc

    chain = _select_chain(model, chain_id, pdb_id)
    retain = {r.upper() for r in retain_hetero_resnames}
    remove = {r.upper() for r in remove_hetero_resnames}

    kept: list[Residue] = []
    summary: dict[str, dict[str, Any]] = {}
    n_input_atoms = 0
    keep_ids: set[tuple] = set()

    for residue in chain:
        n_input_atoms += len(residue)
        hetflag, resseq, icode = residue.id
        name = residue.get_resname().strip().upper()
        is_water = name in WATER_RESNAMES
        is_protein = name in AMINO_ACIDS and hetflag.strip() in ("", "H_MSE")

        if is_water:
            action, reason = "removed", "water"
        elif is_protein:
            action, reason = "retained", "protein"
        elif name in remove:
            action, reason = "removed", "explicit_remove"
        elif name in retain:
            action, reason = "retained", "explicit_retain"
        else:
            action = "retained" if hetero_default == "retain" else "removed"
            reason = "default_policy"

        entry = summary.setdefault(
            name, {"residue_name": name, "action": action, "reason": reason,
                   "residue_count": 0, "atom_count": 0}
        )
        entry["residue_count"] += 1
        entry["atom_count"] += len(residue)

        if action == "retained":
            keep_ids.add((chain.id, residue.id))
            kept.append(
                Residue(
                    chain_id=chain.id,
                    residue_number=int(resseq),
                    insertion_code=icode.strip(),
                    residue_name=name,
                    record_type="ATOM" if is_protein else "HETATM",
                )
            )

    class _Keep(Select):
        def accept_chain(self, c):  # noqa: D102
            return c.id == chain.id

        def accept_residue(self, r):  # noqa: D102
            return (r.get_parent().id, r.id) in keep_ids

        def accept_atom(self, atom):  # noqa: D102
            # Keep only the first altloc: fpocket has no notion of alternates
            # and would otherwise see overlapping atoms in the same place.
            return atom.get_altloc() in (" ", "A")

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    out_path = work / f"{pdb_id.upper()}_{chain.id}_clean.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), select=_Keep())

    n_output_atoms = sum(
        1 for line in out_path.read_text(encoding="utf-8").splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )
    n_protein = sum(1 for r in kept if r.record_type == "ATOM")

    if n_output_atoms < MIN_ATOMS or n_protein == 0:
        raise TargetError(
            f"{pdb_id} chain {chain.id} produced only {n_output_atoms} atoms / "
            f"{n_protein} protein residues after cleaning. Is this a "
            f"nucleic-acid-only structure, or the wrong chain?"
        )

    applied = {
        "input_atom_count": n_input_atoms,
        "output_atom_count": n_output_atoms,
        "output_protein_residue_count": n_protein,
        "removed_water_residue_count": sum(
            e["residue_count"] for e in summary.values() if e["reason"] == "water"
        ),
        "hetero_summary": [e for e in summary.values() if e["reason"] != "protein"],
        "retained_hetero_residues": [
            r.to_dict() for r in kept if r.record_type == "HETATM"
        ],
        "altloc_policy": "first_altloc_only",
        "warnings": [],
    }
    retained_het = applied["retained_hetero_residues"]
    if retained_het:
        applied["warnings"].append(
            f"{len(retained_het)} non-water hetero residue(s) retained "
            f"({sorted({r['residue_name'] for r in retained_het})}); these change "
            f"the detected cavity geometry"
        )

    return PreparedTarget(
        pdb_id=pdb_id.upper(),
        chain_id=chain.id,
        model_index=model_index,
        path=out_path,
        source=source,
        applied=applied,
        residues=kept,
        name=_structure_name(source["path"]),
    )


def _select_chain(model, chain_id: str | None, pdb_id: str):
    if chain_id:
        for chain in model:
            if chain.id == chain_id:
                return chain
        raise TargetError(
            f"{pdb_id}: chain {chain_id!r} not found; available: "
            f"{[c.id for c in model]}"
        )
    # A PDB file often contains several copies of the same protein; the first
    # chain with protein residues is the conventional choice.
    for chain in model:
        if any(r.get_resname().strip().upper() in AMINO_ACIDS for r in chain):
            return chain
    raise TargetError(f"{pdb_id}: no protein chain found in model")


def _structure_name(path: str | Path) -> str | None:
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("TITLE"):
            return line[10:].strip() or None
        if line.startswith("ATOM"):
            break
    return None
