"""fpocket execution and output parsing (spec §5.4).

fpocket is used purely as a measuring instrument: *where are the cavities and
how big are they?* Nothing here interprets its Druggability Score as evidence
of anything.

The parser is deliberately strict. `_info.txt` is a whitespace-sensitive block
format that a loose parser will happily misread into plausible-looking wrong
numbers, and a wrong pocket dimension propagates into every Tier 2 band. Any
non-blank line it cannot account for is an error, not a shrug.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

from ..errors import ExternalToolError

POCKET_HEADER = re.compile(r"^Pocket\s+(\d+)\s*:\s*$")
FIELD_LINE = re.compile(r"^\s*(?P<label>[A-Za-z][^:]*?)\s*:\s*(?P<value>\S.*?)\s*$")

#: Labels we depend on. Missing any of these is a hard failure — fpocket
#: changed its output format and every downstream number is suspect.
REQUIRED_LABELS = ("Score", "Number of Alpha Spheres", "Volume")


@dataclass
class AlphaSphere:
    index: int
    center_A: tuple[float, float, float]
    radius_A: float
    kind: str  # "apolar" | "polar" | "unknown"


@dataclass
class Residue:
    chain_id: str
    residue_number: int
    insertion_code: str
    residue_name: str | None = None
    record_type: str = "ATOM"

    def key(self) -> tuple[str, int, str]:
        """Identity for matching: name is validation, not the primary key."""
        return (self.chain_id.strip(), self.residue_number, self.insertion_code.strip())

    def to_dict(self) -> dict[str, Any]:
        return {
            "chain_id": self.chain_id,
            "residue_number": self.residue_number,
            "insertion_code": self.insertion_code,
            "residue_name": self.residue_name,
            "record_type": self.record_type,
        }


@dataclass
class Pocket:
    index: int
    metrics: dict[str, float | str]
    alpha_spheres: list[AlphaSphere] = field(default_factory=list)
    lining_residues: list[Residue] = field(default_factory=list)

    @property
    def score(self) -> float:
        return float(self.metrics["Score"])

    @property
    def druggability(self) -> float | None:
        value = self.metrics.get("Druggability Score")
        return float(value) if value is not None else None

    @property
    def volume_A3(self) -> float:
        return float(self.metrics["Volume"])

    @property
    def n_alpha_spheres_reported(self) -> int:
        return int(float(self.metrics["Number of Alpha Spheres"]))


def is_available(executable: str = "fpocket") -> bool:
    return shutil.which(executable) is not None


def run(pdb_path: str | Path, executable: str = "fpocket", extra_args: Iterable[str] = ()) -> dict[str, Any]:
    """Run fpocket on a prepared structure. Returns execution provenance."""
    path = Path(pdb_path)
    if not is_available(executable):
        raise ExternalToolError(
            f"{executable} not found on PATH. fpocket has no Windows build; "
            f"generate the target bundle on Linux (see .github/workflows/) and "
            f"score against the bundle instead."
        )
    command = [executable, "-f", str(path), *extra_args]
    proc = subprocess.run(command, capture_output=True, text=True)
    out_dir = path.parent / f"{path.stem}_out"
    if proc.returncode != 0 or not out_dir.is_dir():
        raise ExternalToolError(
            f"fpocket failed (exit {proc.returncode}) on {path}: "
            f"{(proc.stderr or proc.stdout).strip()[:400]}"
        )
    return {
        "status": "success",
        "version": _version(executable),
        "command": command,
        "exit_code": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "output_dir": str(out_dir),
    }


def load_pockets(out_dir: str | Path, stem: str) -> list[Pocket]:
    """Parse a complete `<stem>_out/` directory into Pocket objects.

    Cross-checks the three sources against each other: the pocket indices in
    `_info.txt`, the `pocket{N}_atm.pdb` files and the `pocket{N}_vert.pqr`
    files must agree exactly, and each pocket's reported alpha-sphere count must
    match the number of vertices actually present.
    """
    directory = Path(out_dir)
    info_path = directory / f"{stem}_info.txt"
    if not info_path.exists():
        raise ExternalToolError(f"fpocket info file not found: {info_path}")

    pockets = parse_info(info_path.read_text(encoding="utf-8"))
    info_indices = {p.index for p in pockets}

    pockets_dir = directory / "pockets"
    vert_indices = _indices_from(pockets_dir, "pocket*_vert.pqr", r"pocket(\d+)_vert\.pqr")
    atm_indices = _indices_from(pockets_dir, "pocket*_atm.pdb", r"pocket(\d+)_atm\.pdb")

    if info_indices != vert_indices or info_indices != atm_indices:
        raise ExternalToolError(
            "fpocket output is inconsistent: "
            f"info.txt has pockets {sorted(info_indices)}, "
            f"vert.pqr has {sorted(vert_indices)}, atm.pdb has {sorted(atm_indices)}"
        )

    for pocket in pockets:
        pocket.alpha_spheres = parse_pqr(
            (pockets_dir / f"pocket{pocket.index}_vert.pqr").read_text(encoding="utf-8")
        )
        pocket.lining_residues = parse_lining_residues(
            (pockets_dir / f"pocket{pocket.index}_atm.pdb").read_text(encoding="utf-8")
        )
        reported = pocket.n_alpha_spheres_reported
        if len(pocket.alpha_spheres) != reported:
            raise ExternalToolError(
                f"pocket {pocket.index}: info.txt reports {reported} alpha spheres "
                f"but pocket{pocket.index}_vert.pqr contains {len(pocket.alpha_spheres)}"
            )
        if not pocket.lining_residues:
            raise ExternalToolError(f"pocket {pocket.index} has no lining residues")
    return pockets


def parse_info(text: str) -> list[Pocket]:
    """Parse the blank-line-delimited `<stem>_info.txt` block format."""
    pockets: list[Pocket] = []
    index: int | None = None
    metrics: dict[str, float | str] = {}

    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip()
        if not line.strip():
            continue

        header = POCKET_HEADER.match(line.strip())
        if header:
            if index is not None:
                pockets.append(_finish_pocket(index, metrics))
            index, metrics = int(header.group(1)), {}
            continue

        match = FIELD_LINE.match(line)
        if not match:
            # Never silently skip: an unrecognised line means the format moved.
            raise ExternalToolError(
                f"unparsed line {lineno} in fpocket info output: {raw!r}"
            )
        if index is None:
            raise ExternalToolError(
                f"field line {lineno} appears before any 'Pocket N :' header: {raw!r}"
            )
        label = " ".join(match.group("label").split())
        metrics[label] = _coerce(match.group("value"))

    if index is not None:
        pockets.append(_finish_pocket(index, metrics))
    if not pockets:
        raise ExternalToolError("fpocket info output contained no pockets")
    return pockets


def parse_pqr(text: str) -> list[AlphaSphere]:
    """Alpha-sphere centres and radii from a `pocket{N}_vert.pqr` file.

    PQR is whitespace-delimited with a variable number of leading columns, so
    the numeric fields are taken from the right: `... x y z charge radius`.
    """
    spheres: list[AlphaSphere] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        fields = raw.split()
        if len(fields) < 7:
            raise ExternalToolError(f"malformed PQR line {lineno}: {raw!r}")
        try:
            radius = float(fields[-1])
            x, y, z = (float(v) for v in fields[-5:-2])
        except ValueError as exc:
            raise ExternalToolError(f"malformed PQR line {lineno}: {raw!r}") from exc
        kind = {"APOL": "apolar", "POL": "polar"}.get(fields[2].upper(), "unknown")
        spheres.append(
            AlphaSphere(index=len(spheres), center_A=(x, y, z), radius_A=radius, kind=kind)
        )
    if not spheres:
        raise ExternalToolError("PQR file contained no alpha spheres")
    return spheres


def parse_lining_residues(text: str) -> list[Residue]:
    """Unique residues lining a pocket, from `pocket{N}_atm.pdb`.

    Fixed-column PDB parsing: chain 22, residue number 23–26, insertion code 27
    (1-indexed, per the PDB format spec). Splitting on whitespace breaks when
    residue numbers reach five digits or run into the chain identifier.
    """
    seen: dict[tuple[str, int, str], Residue] = {}
    for lineno, raw in enumerate(text.splitlines(), start=1):
        if not raw.startswith(("ATOM", "HETATM")):
            continue
        if len(raw) < 27:
            raise ExternalToolError(f"truncated PDB line {lineno}: {raw!r}")
        try:
            number = int(raw[22:26])
        except ValueError as exc:
            raise ExternalToolError(f"malformed residue number on line {lineno}: {raw!r}") from exc
        residue = Residue(
            chain_id=raw[21].strip(),
            residue_number=number,
            insertion_code=raw[26].strip(),
            residue_name=raw[17:20].strip(),
            record_type="HETATM" if raw.startswith("HETATM") else "ATOM",
        )
        seen.setdefault(residue.key(), residue)
    return list(seen.values())


def _finish_pocket(index: int, metrics: dict[str, float | str]) -> Pocket:
    missing = [label for label in REQUIRED_LABELS if label not in metrics]
    if missing:
        raise ExternalToolError(
            f"pocket {index} is missing required fpocket field(s) {missing}; "
            f"found {sorted(metrics)}"
        )
    return Pocket(index=index, metrics=metrics)


def _indices_from(directory: Path, glob: str, pattern: str) -> set[int]:
    if not directory.is_dir():
        raise ExternalToolError(f"fpocket pockets directory not found: {directory}")
    regex = re.compile(pattern)
    return {
        int(m.group(1))
        for path in directory.glob(glob)
        if (m := regex.match(path.name))
    }


def _coerce(value: str) -> float | str:
    try:
        return float(value)
    except ValueError:
        return value


def _version(executable: str) -> str | None:
    """The version fpocket reports, via the one parser that knows its quirks.

    This string goes inside the bundle id, so "which fpocket produced these
    numbers" has to be the actual answer rather than the first line the tool
    happened to print.
    """
    from ..provenance import _cli_version

    return _cli_version([executable, "--version"])
