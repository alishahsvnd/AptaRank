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

    if path is not None:
        user = _read_yaml(Path(path))
        data = _deep_merge(data, user)
        sources.append(str(path))

    if overrides:
        data = _deep_merge(data, overrides)
        sources.append("<programmatic override>")

    for item in cli_sets:
        if "=" not in item:
            raise ConfigError(f"--set expects dotted.key=value, got: {item!r}")
        key, raw = item.split("=", 1)
        data = _deep_merge(data, _explode(key.strip(), yaml.safe_load(raw)))
        sources.append(f"--set {item}")

    cfg = Config(data, sources)
    validate(cfg)
    return cfg


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


def _explode(dotted: str, value: Any) -> dict[str, Any]:
    parts = dotted.split(".")
    node: Any = value
    for part in reversed(parts):
        node = {part: node}
    return node
