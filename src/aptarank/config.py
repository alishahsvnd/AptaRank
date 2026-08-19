"""Configuration loading, merging and validation.

One config object is resolved at the start of a run and copied verbatim into
the run artifact. Nothing downstream reads a tunable from anywhere else.

Layering, lowest precedence first:

    configs/default.yaml  ->  user config file  ->  --set key=value overrides
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
from pathlib import Path
from typing import Any, Iterable, Mapping

import yaml

from .errors import ConfigError

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parent.parent
DEFAULT_CONFIG_PATH = REPO_ROOT / "configs" / "default.yaml"

# Criterion definitions are structural, not tunable: they say what each number
# *means*. Only the weights (config) and the reference distribution (corpus)
# are adjustable.
#
#   mode "lower"     -> one-sided percentile, sign-flipped (lower value is better)
#   mode "higher"    -> one-sided percentile (higher value is better)
#   mode "two_sided" -> typicality; 1.0 at the corpus median, 0 at either tail
CRITERIA: dict[str, dict[str, str]] = {
    "mfe_norm": {
        "mode": "lower",
        "label": "Structure stability",
        "blurb": "length-normalised folding free energy",
    },
    "ensemble_defect": {
        "mode": "lower",
        "label": "Fold definition",
        "blurb": "how much the ensemble disagrees with the MFE structure",
    },
    "positional_entropy_mean": {
        "mode": "lower",
        "label": "Fold certainty",
        "blurb": "average per-position folding uncertainty",
    },
    "stem_fraction": {
        "mode": "two_sided",
        "label": "Structural composition",
        "blurb": "fraction of positions that are paired",
    },
    "gc_fraction": {
        "mode": "two_sided",
        "label": "Sequence composition",
        "blurb": "G+C content",
    },
}

_MISSING = object()

COMPOSITE_METHODS = ("corpus_weighted_mean", "batch_rank_aggregation")
RUN_MODES = ("standard", "fast")
PRIMARY_DESCRIPTORS = ("flexible", "extended")

#: Binding modes the comparison layer implements. `groove` is described in the
#: paper as future work and is named here so asking for it fails with that
#: answer rather than with "unknown value".
BINDING_MODES = ("pocket", "surface")
FUTURE_BINDING_MODES = ("groove",)
TARGET_SOURCES = ("pdb", "alphafold")
FOOTPRINT_MODELS = ("radius_of_gyration", "length")

#: Keys renamed by the refinements spec, old -> new. Accepted and rewritten so
#: existing configs, saved runs and CI invocations keep working; the rewrite is
#: recorded in `Config.sources` rather than applied silently.
LEGACY_KEYS: tuple[tuple[str, str], ...] = (
    ("tier2.a_per_nt", "tier2.geometry.a_per_nt_ss"),
    ("tier2.target.pdb_id", "tier2.target.id"),
    ("tier2.target.active_site_residues", "tier2.target.target_site_residues"),
)


class Config:
    """Immutable-by-convention view over the resolved configuration tree."""

    def __init__(self, data: Mapping[str, Any], sources: Iterable[str] = ()) -> None:
        self._data = copy.deepcopy(dict(data))
        self.sources = list(sources)

    # -- access ----------------------------------------------------------

    def get(self, dotted: str, default: Any = _MISSING) -> Any:
        node: Any = self._data
        for part in dotted.split("."):
            if not isinstance(node, Mapping) or part not in node:
                if default is _MISSING:
                    raise ConfigError(f"missing config key: {dotted}")
                return default
            node = node[part]
        return copy.deepcopy(node) if isinstance(node, (dict, list)) else node

    def __getitem__(self, key: str) -> Any:
        return self.get(key)

    def as_dict(self) -> dict[str, Any]:
        """The full resolved tree, for embedding in the run artifact."""
        return copy.deepcopy(self._data)

    def with_overrides(self, overrides: Mapping[str, Any]) -> "Config":
        merged = _deep_merge(self._data, overrides)
        return Config(merged, self.sources + ["<programmatic override>"])

    # -- derived ---------------------------------------------------------

    @property
    def is_fast(self) -> bool:
        return self.get("run.mode") == "fast"

    @property
    def shuffles_enabled(self) -> bool:
        return bool(self.get("tier1.shuffle.enabled")) and not self.is_fast

    def active_criteria(self) -> list[str]:
        """Criteria with a non-zero weight, in canonical order."""
        weights = self.get("tier1.weights")
        return [name for name in CRITERIA if float(weights.get(name, 0.0)) > 0.0]

    def scoring_signature(self) -> str:
        """Hash of everything that changes what a tier1_score *means*.

        Two runs may only be compared (evaluation E1, cross-target work) when
        their scoring signatures match. Deliberately excludes cosmetic keys
        such as output paths and worker counts.
        """
        relevant = {
            "criteria": {k: v["mode"] for k, v in CRITERIA.items()},
            "weights": self.get("tier1.weights"),
            "composite_method": self.get("tier1.composite.method"),
            "min_length": self.get("input.min_length"),
            "max_length": self.get("input.max_length"),
            "n_ensemble_samples": self.get("tier1.n_ensemble_samples"),
        }
        blob = json.dumps(relevant, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(blob.encode()).hexdigest()[:16]

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Config(sources={self.sources!r})"


def load_config(
    path: str | os.PathLike[str] | None = None,
    overrides: Mapping[str, Any] | None = None,
    cli_sets: Iterable[str] = (),
) -> Config:
    """Resolve the configuration from defaults, an optional file and overrides.

    `cli_sets` accepts `dotted.key=value` strings; values are parsed as YAML so
    `--set tier1.shuffle.n_shuffles=99` yields an int, not a string.
    """
    data = _read_yaml(DEFAULT_CONFIG_PATH)
    sources = [str(DEFAULT_CONFIG_PATH)]

    def layer(values: Mapping[str, Any], label: str) -> None:
        # Legacy names are rewritten per layer, so a renamed key is compared
        # against what the *user* set rather than against the default value the
        # new key always carries.
        nonlocal data
        migrated, notes = migrate_legacy_keys(values)
        data = _deep_merge(data, migrated)
        sources.append(label)
        sources.extend(notes)

    if path is not None:
        layer(_read_yaml(Path(path)), str(path))

    if overrides:
        layer(overrides, "<programmatic override>")

    for item in cli_sets:
        if "=" not in item:
            raise ConfigError(f"--set expects dotted.key=value, got: {item!r}")
        key, raw = item.split("=", 1)
        layer(_explode(key.strip(), yaml.safe_load(raw)), f"--set {item}")

    cfg = Config(data, sources)
    validate(cfg)
    return cfg


def migrate_legacy_keys(data: Mapping[str, Any]) -> tuple[dict[str, Any], list[str]]:
    """Rewrite pre-refinement key names, returning what was rewritten.

    A renamed key that silently did nothing would be the worst outcome here:
    `active_site_residues` still parsing but no longer selecting the site would
    point the tool at the wrong part of the protein and say nothing about it.
    """
    out = copy.deepcopy(dict(data))
    notes: list[str] = []
    for old, new in LEGACY_KEYS:
        value = _dig(out, old)
        if value is _MISSING or value in (None, [], {}):
            continue
        if _dig(out, new) not in (_MISSING, None, [], {}):
            raise ConfigError(
                f"both {old!r} (renamed) and {new!r} are set to non-empty values; "
                f"remove {old!r}"
            )
        out = _deep_merge(out, _explode(new, value))
        _drop(out, old)
        notes.append(f"<renamed {old} -> {new}>")
    return out, notes


#: The target-input contract (refinements §3.2), as the biologist writes it:
#:
#:     target_name: IGFBP3
#:     source: pdb
#:     id: 7WRQ
#:     chain: B
#:     binding_mode: surface
#:     partner_chain: C
#:     strip_hetatm: true
#:     target_site_residues: [7, 8, 9, ...]
#:
#: Left is the file's key, right is where it lands in the resolved config.
TARGET_SPEC_KEYS: dict[str, str] = {
    "target_name": "tier2.target.name",
    "name": "tier2.target.name",
    "source": "tier2.target.source",
    "id": "tier2.target.id",
    "chain": "tier2.target.chain",
    "binding_mode": "tier2.binding_mode",
    "partner_chain": "tier2.target.partner_chains",
    "partner_chains": "tier2.target.partner_chains",
    "strip_hetatm": "tier2.target.strip_hetatm",
    "target_site_residues": "tier2.target.target_site_residues",
    "model": "tier2.target.model",
    "retain_hetero_resnames": "tier2.target.retain_hetero_resnames",
    "allow_zero_overlap_fallback": "tier2.target.allow_zero_overlap_fallback",
}


def load_target_spec(path: str | os.PathLike[str]) -> dict[str, Any]:
    """Read a target description file into config overrides."""
    target_path = Path(path)
    if not target_path.is_file():
        raise ConfigError(f"target file not found: {target_path}")
    return parse_target_spec(target_path.read_text(encoding="utf-8"), str(target_path))


def parse_target_spec(text: str, label: str = "target description") -> dict[str, Any]:
    """Turn a target description into config overrides.

    An unknown key is an error rather than a shrug: a misspelled
    `target_site_residues` that was quietly ignored would point the tool at a
    site the biologist never chose, and nothing downstream would say so.
    """
    try:
        raw = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ConfigError(f"{label}: could not be read as a target file: {exc}") from exc
    if not isinstance(raw, Mapping):
        raise ConfigError(
            f"{label}: a target description is a list of `key: value` lines, e.g.\n"
            f"    target_name: IGFBP3\n    source: pdb\n    id: 7WRQ\n"
            f"    chain: B\n    binding_mode: surface"
        )

    unknown = [k for k in raw if k not in TARGET_SPEC_KEYS]
    if unknown:
        raise ConfigError(
            f"{label}: unrecognised key(s) {sorted(unknown)}. Supported keys: "
            f"{sorted(set(TARGET_SPEC_KEYS))}"
        )

    overrides: dict[str, Any] = {"tier2": {"enabled": True}}
    for key, value in raw.items():
        dotted = TARGET_SPEC_KEYS[key]
        if dotted == "tier2.target.partner_chains" and not isinstance(value, list):
            value = [value]
        if dotted == "tier2.target.target_site_residues":
            value = _residue_list(value, label)
        overrides = _deep_merge(overrides, _explode(dotted, value))
    # `_dig` returns a sentinel when the key is absent, and the sentinel is
    # truthy — so this must test for it explicitly.
    if _dig(overrides, "tier2.target.id") in (_MISSING, None, ""):
        raise ConfigError(
            f"{label}: a target description must set `id` (the PDB ID or UniProt "
            f"accession)"
        )
    return overrides


def _residue_list(value: Any, label: str) -> list[int]:
    if value in (None, ""):
        return []
    if not isinstance(value, (list, tuple)):
        raise ConfigError(
            f"{label}: target_site_residues must be a list, e.g. [42, 87, 119]"
        )
    out = []
    for item in value:
        try:
            out.append(int(item))
        except (TypeError, ValueError) as exc:
            raise ConfigError(
                f"{label}: target_site_residues must be plain integers in the "
                f"author numbering of the chosen chain — write 42, not {item!r}. "
                f"The amino-acid letter is redundant; the number and chain "
                f"select the residue."
            ) from exc
    return out


def validate(cfg: Config) -> None:
    """Fail fast and specifically on a config that cannot produce a valid run."""
    mode = cfg.get("run.mode")
    if mode not in RUN_MODES:
        raise ConfigError(f"run.mode must be one of {RUN_MODES}, got {mode!r}")

    method = cfg.get("tier1.composite.method")
    if method not in COMPOSITE_METHODS:
        raise ConfigError(
            f"tier1.composite.method must be one of {COMPOSITE_METHODS}, got {method!r}"
        )

    min_len, max_len = cfg.get("input.min_length"), cfg.get("input.max_length")
    if not (0 < min_len <= max_len):
        raise ConfigError(f"require 0 < min_length <= max_length, got {min_len}, {max_len}")

    weights = cfg.get("tier1.weights")
    unknown = set(weights) - set(CRITERIA)
    if unknown:
        raise ConfigError(f"unknown criteria in tier1.weights: {sorted(unknown)}")
    if any(float(w) < 0 for w in weights.values()):
        raise ConfigError("tier1.weights must be non-negative")
    if not cfg.active_criteria():
        raise ConfigError("tier1.weights leaves no criterion with non-zero weight")

    structural = cfg.get("tier1.shuffle.structural_criteria")
    unknown = set(structural) - set(CRITERIA)
    if unknown:
        raise ConfigError(
            f"unknown criteria in tier1.shuffle.structural_criteria: {sorted(unknown)}"
        )
    if not structural:
        raise ConfigError("tier1.shuffle.structural_criteria must not be empty")

    alpha = float(cfg.get("tier1.shuffle.alpha"))
    n_shuf = int(cfg.get("tier1.shuffle.n_shuffles"))
    if not 0 < alpha < 1:
        raise ConfigError(f"tier1.shuffle.alpha must be in (0,1), got {alpha}")
    if n_shuf < 1:
        raise ConfigError("tier1.shuffle.n_shuffles must be >= 1")
    # An exact Monte-Carlo p-value with M controls cannot go below 1/(M+1).
    if cfg.get("tier1.shuffle.enabled") and 1.0 / (n_shuf + 1) > alpha:
        raise ConfigError(
            f"tier1.shuffle.n_shuffles={n_shuf} can never reach alpha={alpha}: the "
            f"smallest attainable p-value is {1.0 / (n_shuf + 1):.4f}. "
            f"Use n_shuffles >= {int(round(1 / alpha)) - 1} or raise alpha."
        )

    if int(cfg.get("tier1.n_ensemble_samples")) < 1:
        raise ConfigError("tier1.n_ensemble_samples must be >= 1")

    descriptor = cfg.get("tier2.primary_descriptor")
    if descriptor not in PRIMARY_DESCRIPTORS:
        raise ConfigError(
            f"tier2.primary_descriptor must be one of {PRIMARY_DESCRIPTORS}, got {descriptor!r}"
        )

    mode = cfg.get("tier2.binding_mode")
    if mode in FUTURE_BINDING_MODES:
        raise ConfigError(
            f"binding mode {mode!r} is described as future work and is not "
            f"implemented. Available modes: {BINDING_MODES}."
        )
    if mode not in BINDING_MODES:
        raise ConfigError(
            f"tier2.binding_mode must be one of {BINDING_MODES}, got {mode!r}"
        )

    source = cfg.get("tier2.target.source")
    if source not in TARGET_SOURCES:
        raise ConfigError(
            f"tier2.target.source must be one of {TARGET_SOURCES}, got {source!r}"
        )

    for key in ("a_per_nt_ss", "a_per_bp_helix", "footprint_scale"):
        if float(cfg.get(f"tier2.geometry.{key}")) <= 0:
            raise ConfigError(f"tier2.geometry.{key} must be > 0")

    if float(cfg.get("tier2.surface.sigma_A2")) <= 0:
        raise ConfigError("tier2.surface.sigma_A2 must be > 0")

    footprint = cfg.get("tier2.surface.footprint_model")
    if footprint not in FOOTPRINT_MODELS:
        raise ConfigError(
            f"tier2.surface.footprint_model must be one of {FOOTPRINT_MODELS}, "
            f"got {footprint!r}"
        )

    weights = cfg.get("tier2.surface.weights")
    if any(float(w) < 0 for w in weights.values()):
        raise ConfigError("tier2.surface.weights must be non-negative")
    if sum(float(w) for w in weights.values()) <= 0:
        raise ConfigError("tier2.surface.weights leaves no signal with a positive weight")

    residues = cfg.get("tier2.target.target_site_residues")
    if not isinstance(residues, list):
        raise ConfigError(
            "tier2.target.target_site_residues must be a list of integer residue "
            f"numbers, got {type(residues).__name__}"
        )

    partners = cfg.get("tier2.target.partner_chains")
    chain = cfg.get("tier2.target.chain")
    if chain is not None and chain in (partners or []):
        raise ConfigError(
            f"tier2.target.chain {chain!r} is also listed in partner_chains; the "
            f"target chain cannot be its own binding partner"
        )

    bands = cfg.get("tier2.band_percentiles")
    if not 0 < bands["moderate"] < bands["strong"] < 1:
        raise ConfigError(
            "require 0 < tier2.band_percentiles.moderate < strong < 1, "
            f"got {bands}"
        )

    if float(cfg.get("tier2.sigma_A")) <= 0:
        raise ConfigError("tier2.sigma_A must be > 0")


# -- helpers -------------------------------------------------------------


def _read_yaml(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise ConfigError(f"config file not found: {path}")
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"config file must contain a mapping at top level: {path}")
    return data


def _deep_merge(base: Mapping[str, Any], override: Mapping[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(dict(base))
    for key, value in override.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), Mapping):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


def _dig(data: Mapping[str, Any], dotted: str) -> Any:
    node: Any = data
    for part in dotted.split("."):
        if not isinstance(node, Mapping) or part not in node:
            return _MISSING
        node = node[part]
    return node


def _drop(data: dict[str, Any], dotted: str) -> None:
    parts = dotted.split(".")
    node: Any = data
    for part in parts[:-1]:
        if not isinstance(node, dict) or part not in node:
            return
        node = node[part]
    if isinstance(node, dict):
        node.pop(parts[-1], None)


def _explode(dotted: str, value: Any) -> dict[str, Any]:
    parts = dotted.split(".")
    node: Any = value
    for part in reversed(parts):
        node = {part: node}
    return node
