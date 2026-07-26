"""AptaRank — two-tier evaluation and ranking of generated RNA aptamer candidates.

Tier 1 (target-agnostic intrinsic quality) produces the ranking.
Tier 2 (target-aware geometric plausibility) annotates it and never reorders it.

Nothing in this package predicts binding. Target-aware output is plausibility
evidence only.
"""

__version__ = "0.1.0"
ARTIFACT_SCHEMA_VERSION = "1.0.0"

TIER2_CAVEAT = (
    "Tier 2 reflects geometric plausibility only — it is not evidence of binding."
)
