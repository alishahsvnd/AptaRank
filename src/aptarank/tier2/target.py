"""Fetch and prepare a target structure (spec §5.3, refinements §3).

The user gives an identifier, a chain and a binding mode; everything here runs
server-side. Nobody uploads a prepared structure, and nobody needs Linux tools
on their own machine.

Two parts of this are scientific decisions rather than plumbing:

**Hetero-atom policy.** Non-water HETATM are kept by default because the enzyme
demo targets have a metal ion sitting in the functional cavity, and removing it
changes the cavity geometry. But "keep everything" is too broad: crystallisation
buffers, cryoprotectants and co-crystallised inhibitors also arrive as HETATM and
can *block* the very cavity we are trying to measure. So the policy is explicit
per target and whatever was actually kept or dropped is recorded residue by
residue in the bundle.

**Partner-chain ordering (§3.4).** For an interface target the binding partner
sits *on* the surface we want to measure. The partner is therefore used first —
to confirm the configured binding-site residues really are the interface — and
only then stripped, so the surface is exposed for measurement. The partner's
location informs *where* to look; the partner itself must not be there when we
look.
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..errors import TargetError
from ..provenance import sha256_file
from .fpocket import Residue

RCSB_URL = "https://files.rcsb.org/download/{ident}.pdb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{ident}"
#: Only a fallback. The version suffix moves (v4 was current when the spec was
#: written and 404s now), so the API is asked for the real URL first.
ALPHAFOLD_FILE = "https://alphafold.ebi.ac.uk/files/AF-{ident}-F1-model_v{version}.pdb"
ALPHAFOLD_FALLBACK_VERSIONS = (6, 5, 4)

WATER_RESNAMES = ("HOH", "DOD", "WAT")
MIN_ATOMS = 100
#: Heavy-atom distance at which two residues are called "in contact". The
#: 4.0-5.0 A band is the usual convention for an interface residue.
CONTACT_CUTOFF_A = 4.5

AMINO_ACIDS = frozenset(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL "
    "MSE SEC PYL".split()
)


@dataclass
class PreparedTarget:
    """A single cleaned chain, ready for measurement, plus what was done to it."""

    identifier: str
    chain_id: str
    model_index: int
    path: Path
    source: dict[str, Any]
    applied: dict[str, Any]
    residues: list[Residue] = field(default_factory=list)
    name: str | None = None
    #: "experimental" (PDB) or "predicted" (AlphaFold). Carried all the way to
    #: the run artifact: a measurement taken on a model is not a measurement.
    structure_kind: str = "experimental"
    source_kind: str = "pdb"
    site_residues: list[Residue] = field(default_factory=list)
    partner_evidence: dict[str, Any] = field(default_factory=dict)

    @property
    def pdb_id(self) -> str:
        """Backwards-compatible alias — bundles v1 called this pdb_id."""
        return self.identifier


def fetch_structure(
    identifier: str, cache_dir: str | Path, source: str = "pdb"
) -> dict[str, Any]:
    """Download a structure from RCSB or AlphaFold DB, cached by source and id."""
    if source == "pdb":
        return _fetch_pdb(identifier, cache_dir)
    if source == "alphafold":
        return _fetch_alphafold(identifier, cache_dir)
    raise TargetError(f"unknown target source {source!r}; expected 'pdb' or 'alphafold'")


def _fetch_pdb(identifier: str, cache_dir: str | Path) -> dict[str, Any]:
    ident = identifier.strip().upper()
    if len(ident) != 4 or not ident.isalnum():
        raise TargetError(f"expected a 4-character PDB ID, got {identifier!r}")
    url = RCSB_URL.format(ident=ident)
    path = Path(cache_dir) / f"pdb_{ident}.pdb"
    return _download(url, path, ident, source_kind="pdb", structure_kind="experimental")


def _fetch_alphafold(identifier: str, cache_dir: str | Path) -> dict[str, Any]:
    ident = identifier.strip().upper()
    if not ident.isalnum() or not 6 <= len(ident) <= 10:
        raise TargetError(
            f"expected a UniProt accession for an AlphaFold target, got {identifier!r}"
        )
    path = Path(cache_dir) / f"alphafold_{ident}.pdb"

    url, version = _alphafold_url(ident)
    result = _download(url, path, ident, source_kind="alphafold",
                       structure_kind="predicted")
    result["model_version"] = version
    return result


def _alphafold_url(ident: str) -> tuple[str, int | None]:
    """Ask AlphaFold DB for the current model URL rather than pinning a version."""
    try:
        with urllib.request.urlopen(ALPHAFOLD_API.format(ident=ident), timeout=60) as fh:
            payload = json.loads(fh.read().decode("utf-8"))
        entry = payload[0] if isinstance(payload, list) and payload else {}
        url = entry.get("pdbUrl")
        if url:
            return url, entry.get("latestVersion")
    except (urllib.error.URLError, OSError, ValueError, KeyError, IndexError):
        pass  # fall through to the version sweep below
    for version in ALPHAFOLD_FALLBACK_VERSIONS:
        candidate = ALPHAFOLD_FILE.format(ident=ident, version=version)
        try:
            request = urllib.request.Request(candidate, method="HEAD")
            with urllib.request.urlopen(request, timeout=30):
                return candidate, version
        except (urllib.error.URLError, OSError):
            continue
    raise TargetError(
        f"no AlphaFold model found for {ident}. Check the UniProt accession, or "
        f"use an experimental structure (source: pdb)."
    )


def _download(
    url: str, path: Path, ident: str, source_kind: str, structure_kind: str
) -> dict[str, Any]:
    path.parent.mkdir(parents=True, exist_ok=True)
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
        "identifier": ident,
        "source": source_kind,
        "structure_kind": structure_kind,
        "url": url,
        "format": "pdb",
        "path": str(path),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    }


def prepare(
    identifier: str,
    cache_dir: str | Path,
    work_dir: str | Path,
    source: str = "pdb",
    chain_id: str | None = None,
    partner_chains: Sequence[str] = (),
    target_site_residues: Sequence[int] = (),
    model_index: int = 0,
    strip_hetatm: bool = False,
    retain_hetero_resnames: Sequence[str] = (),
    remove_hetero_resnames: Sequence[str] = WATER_RESNAMES,
    hetero_default: str = "retain",
) -> PreparedTarget:
    """Select one chain, confirm the site, strip partners, write a cleaned PDB."""
    from Bio.PDB import PDBIO, PDBParser, Select  # imported here: heavy

    fetched = fetch_structure(identifier, cache_dir, source=source)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(identifier.upper(), fetched["path"])

    try:
        model = structure[model_index]
    except KeyError as exc:
        raise TargetError(f"{identifier}: no model {model_index}") from exc

    if source == "alphafold" and chain_id not in (None, "A"):
        raise TargetError(
            f"AlphaFold models contain a single chain 'A'; chain {chain_id!r} was "
            f"requested. Use chain A, or an experimental structure."
        )

    chain = _select_chain(model, chain_id, identifier)
    present_chains = [c.id for c in model]

    # -- §3.4: use the partner while it is still here, then remove it --------
    site_residues = _resolve_site_residues(chain, target_site_residues, identifier)
    partners = [c for c in partner_chains if c and c != chain.id]
    missing_partners = [c for c in partners if c not in present_chains]
    if missing_partners:
        raise TargetError(
            f"{identifier}: partner chain(s) {missing_partners} are not in the "
            f"structure; available chains: {present_chains}"
        )
    partner_evidence = _partner_contacts(model, chain, partners, site_residues)

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
        elif strip_hetatm:
            action, reason = "removed", "strip_hetatm"
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
            # Every other chain — partners included — is dropped here. This is
            # a label filter, not structural surgery.
            return c.id == chain.id

        def accept_residue(self, r):  # noqa: D102
            return (r.get_parent().id, r.id) in keep_ids

        def accept_atom(self, atom):  # noqa: D102
            # Keep only the first altloc: fpocket has no notion of alternates
            # and would otherwise see overlapping atoms in the same place.
            return atom.get_altloc() in (" ", "A")

    work = Path(work_dir)
    work.mkdir(parents=True, exist_ok=True)
    out_path = work / f"{identifier.upper()}_{chain.id}_clean.pdb"
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
            f"{identifier} chain {chain.id} produced only {n_output_atoms} atoms / "
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
        "strip_hetatm": bool(strip_hetatm),
        "altloc_policy": "first_altloc_only",
        "chains_in_source": present_chains,
        "chains_removed": [c for c in present_chains if c != chain.id],
        "was_multi_chain": len(present_chains) > 1,
        "partner_chains_removed": partners,
        "warnings": [],
    }
    retained_het = applied["retained_hetero_residues"]
    if retained_het:
        applied["warnings"].append(
            f"{len(retained_het)} non-water hetero residue(s) retained "
            f"({sorted({r['residue_name'] for r in retained_het})}); these change "
            f"the detected cavity geometry"
        )
    if applied["was_multi_chain"]:
        # The induced-fit caveat (§3.4). Surfaced, not solved — and stated
        # precisely: the coordinates are still the ones the complex was solved
        # in, so this is the bound-state interface with the partner deleted,
        # not a relaxed apo structure.
        applied["warnings"].append(
            "The target chain was isolated computationally from a multi-chain "
            "complex. Its coordinates remain in the bound-state conformation, "
            "but the partner has been removed, exposing the interaction "
            "surface. The measured geometry therefore reflects the "
            "partner-free bound-state interface, not a relaxed unbound "
            "structure."
        )

    return PreparedTarget(
        identifier=identifier.upper(),
        chain_id=chain.id,
        model_index=model_index,
        path=out_path,
        source=fetched,
        applied=applied,
        residues=kept,
        name=_structure_name(fetched["path"]),
        structure_kind=fetched["structure_kind"],
        source_kind=source,
        site_residues=site_residues,
        partner_evidence=partner_evidence,
    )


def _resolve_site_residues(
    chain, requested: Iterable[int], identifier: str
) -> list[Residue]:
    """Match configured residue numbers against the chain's author numbering.

    Residue numbers are the depositor's labels, not sequential positions: a
    structure resolving residues 28-291 numbers its first residue 28. Biopython
    preserves that, so published residue tables line up — nothing is re-indexed.
    """
    wanted = []
    for item in requested or []:
        try:
            wanted.append(int(item))
        except (TypeError, ValueError) as exc:
            raise TargetError(
                f"binding-site residues must be plain integers in the author "
                f"numbering of chain {chain.id} (write 42, not 'K42'); got {item!r}"
            ) from exc
    if not wanted:
        return []

    by_number: dict[int, list] = {}
    insertion_coded = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        by_number.setdefault(int(residue.id[1]), []).append(residue)
        if residue.id[2].strip():
            insertion_coded.append(f"{residue.id[1]}{residue.id[2]}")

    missing = [n for n in wanted if n not in by_number]
    if missing:
        numbers = sorted(by_number)
        raise TargetError(
            f"{identifier} chain {chain.id}: binding-site residue(s) {missing[:8]} "
            f"are not present in this chain (it numbers {numbers[0]}-{numbers[-1]}, "
            f"{len(numbers)} residues resolved). Check the chain selector and that "
            f"the numbers use this structure's author numbering."
        )
    if insertion_coded:
        # Rare, and silently mismatching on one would point the tool at the
        # wrong residue. §3.6 says assert rather than assume.
        raise TargetError(
            f"{identifier} chain {chain.id} contains insertion-coded residues "
            f"({insertion_coded[:5]}), which this residue selector cannot address "
            f"unambiguously."
        )

    out = []
    for number in wanted:
        residue = by_number[number][0]
        out.append(
            Residue(
                chain_id=chain.id,
                residue_number=number,
                insertion_code="",
                residue_name=residue.get_resname().strip().upper(),
                record_type="ATOM",
            )
        )
    return out


def _partner_contacts(model, chain, partners: Sequence[str], site_residues) -> dict[str, Any]:
    """Which configured site residues actually touch the partner chain(s).

    Computed *before* the partner is stripped, which is the whole point of the
    ordering: afterwards there is nothing left to measure against. Disagreement
    here is reported rather than corrected — the configured list is the
    biologist's assertion, and quietly substituting our own would replace their
    judgement with a distance cutoff.
    """
    if not partners:
        return {"partner_chains": [], "computed": False}

    from Bio.PDB import NeighborSearch

    partner_atoms = [
        atom
        for cid in partners
        for residue in model[cid]
        if residue.id[0] == " " or residue.get_resname().strip().upper() not in WATER_RESNAMES
        for atom in residue
    ]
    if not partner_atoms:
        return {"partner_chains": list(partners), "computed": False,
                "reason": "partner chains contain no atoms"}

    search = NeighborSearch(partner_atoms)
    contacting: list[dict[str, Any]] = []
    for residue in chain:
        if residue.id[0] != " ":
            continue
        best = None
        for atom in residue:
            for neighbour in search.search(atom.coord, CONTACT_CUTOFF_A):
                distance = float(atom - neighbour)
                if best is None or distance < best:
                    best = distance
        if best is not None:
            contacting.append(
                {
                    "residue_number": int(residue.id[1]),
                    "residue_name": residue.get_resname().strip().upper(),
                    "min_distance_A": round(best, 2),
                }
            )

    contact_numbers = {c["residue_number"] for c in contacting}
    configured = {r.residue_number for r in site_residues}
    return {
        "partner_chains": list(partners),
        "computed": True,
        "cutoff_A": CONTACT_CUTOFF_A,
        "definition": (
            f"chain {chain.id} residues with any atom within {CONTACT_CUTOFF_A} A "
            f"of any atom of chain(s) {','.join(partners)}"
        ),
        "n_interface_residues": len(contacting),
        "interface_residues": sorted(contact_numbers),
        "interface_detail": contacting,
        "configured_site_residues": sorted(configured),
        "configured_in_interface": sorted(configured & contact_numbers),
        "configured_not_in_interface": sorted(configured - contact_numbers),
        "interface_not_configured": sorted(contact_numbers - configured),
    }


def _select_chain(model, chain_id: str | None, identifier: str):
    if chain_id:
        for chain in model:
            if chain.id == chain_id:
                return chain
        raise TargetError(
            f"{identifier}: chain {chain_id!r} not found; available: "
            f"{[c.id for c in model]}"
        )
    # A PDB file often contains several copies of the same protein; the first
    # chain with protein residues is the conventional choice.
    for chain in model:
        if any(r.get_resname().strip().upper() in AMINO_ACIDS for r in chain):
            return chain
    raise TargetError(f"{identifier}: no protein chain found in model")


def _structure_name(path: str | Path) -> str | None:
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if line.startswith("TITLE"):
            return line[10:].strip() or None
        if line.startswith("ATOM"):
            break
    return None
