"""Tier 2 — target-aware geometric plausibility.

Tier 2 annotates the Tier 1 ranking. It never reorders it (spec §6).

Nothing here is a binding prediction. The output is: *is this candidate's loop
geometry plausible against a cavity detected on this target, compared to
shuffled controls?* That is one signal among several, and the code, the field
names and the UI copy all have to keep saying so.

Structure: the heavy external tools (fpocket, PDB2PQR, APBS) run once per
target on Linux and emit an immutable, checksummed **target bundle**. Scoring
consumes the bundle and needs no external tools at all.
"""
