"""Evaluation experiments E1–E5 (spec §8).

Everything here reads run artifacts or the corpus and produces the paper's
quantitative claims. There is no separate analysis pipeline: the numbers in the
paper come from the same code path as the numbers in the demo.

One methodological rule governs the whole module: **composite weights are never
fitted to these groups.** Tuning weights to maximise separation and then
reporting that separation would inflate the result by construction.
"""

from .groups import ComparisonGroups, build_groups  # noqa: F401
