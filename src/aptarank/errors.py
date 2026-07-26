"""Typed errors.

Split by who is responsible for the failure, so the CLI can print a message
that tells the user what to do rather than a traceback.
"""

from __future__ import annotations


class AptaRankError(Exception):
    """Base class for every error this package raises deliberately."""


class ConfigError(AptaRankError):
    """The configuration is malformed, incomplete or internally inconsistent."""


class InputError(AptaRankError):
    """The user-supplied candidate file could not be read or is unusable."""


class CorpusError(AptaRankError):
    """The reference corpus is missing, malformed, or incompatible with the cache."""


class FoldingError(AptaRankError):
    """ViennaRNA or forgi returned something we refuse to proceed on."""


class TargetError(AptaRankError):
    """Target structure fetch, preparation, or cavity detection failed."""


class ExternalToolError(AptaRankError):
    """An external CLI (fpocket, pdb2pqr, apbs) is missing or failed."""


class ArtifactError(AptaRankError):
    """A run artifact could not be written, read, or failed validation."""
