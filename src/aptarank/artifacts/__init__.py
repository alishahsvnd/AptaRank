"""Run artifacts: the only interface between the pipeline and everything else.

The dashboard, the evaluation scripts and the paper figures all read artifacts
and hold no computation of their own. That decoupling is what lets the paper's
figures be generated from the same runs the live demo uses.
"""

from .io import build_artifact, read_artifact, write_artifact  # noqa: F401
