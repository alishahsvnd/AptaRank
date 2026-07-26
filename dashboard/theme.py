"""Colour and chart styling for the dashboard.

Two deliberate choices, both about honesty rather than taste:

* **Bands use one ordinal blue ramp, not a traffic light.** "Strong" means the
  candidate's loop geometry agrees with this target's cavity better than most
  shuffled controls. It does not mean "good candidate", and it certainly does
  not mean "binds". Green-amber-red would say exactly that, so bands are steps
  of a single hue — ordered, but carrying no verdict.
* **"Not evaluated" is neutral grey, never the weak colour.** A candidate below
  the Tier 2 cut has no target-aware evidence at all; painting it as a weak
  match would fabricate evidence.

Palettes are validated ordinal ramps (single hue, monotone lightness, visible
step gaps, light end clear of the surface) in both light and dark mode.
"""

from __future__ import annotations

from typing import Any

BAND_ORDER = ["strong", "moderate", "weak", "not_evaluated"]

LIGHT = {
    "surface": "#fcfcfb",
    "text_primary": "#0b0b0b",
    "text_secondary": "#52514e",
    "muted": "#898781",
    "grid": "#e1e0d9",
    "axis": "#c3c2b7",
    "series": "#2a78d6",
    "bands": {
        "strong": "#104281",
        "moderate": "#2a78d6",
        "weak": "#86b6ef",
        "not_evaluated": "#898781",
    },
}

DARK = {
    "surface": "#1a1a19",
    "text_primary": "#ffffff",
    "text_secondary": "#c3c2b7",
    "muted": "#898781",
    "grid": "#2c2c2a",
    "axis": "#383835",
    "series": "#3987e5",
    "bands": {
        "strong": "#9ec5f4",
        "moderate": "#3987e5",
        "weak": "#184f95",
        "not_evaluated": "#898781",
    },
}

BAND_LABEL = {
    "strong": "strong",
    "moderate": "moderate",
    "weak": "weak",
    "not_evaluated": "not evaluated",
}


def palette(theme_type: str | None) -> dict[str, Any]:
    return DARK if (theme_type or "light").lower() == "dark" else LIGHT


def band_scale(colors: dict[str, Any]) -> tuple[list[str], list[str]]:
    return BAND_ORDER, [colors["bands"][b] for b in BAND_ORDER]


def chart_config(colors: dict[str, Any]) -> dict[str, Any]:
    """Recessive grid and axes; the data carries the ink."""
    return {
        "background": "transparent",
        "axis": {
            "domainColor": colors["axis"],
            "gridColor": colors["grid"],
            "gridWidth": 1,
            "labelColor": colors["muted"],
            "tickColor": colors["axis"],
            "titleColor": colors["text_secondary"],
            "titleFontWeight": "normal",
            "labelFont": "system-ui, -apple-system, Segoe UI, sans-serif",
            "titleFont": "system-ui, -apple-system, Segoe UI, sans-serif",
        },
        "legend": {
            "labelColor": colors["text_secondary"],
            "titleColor": colors["text_secondary"],
            "titleFontWeight": "normal",
            "symbolStrokeWidth": 0,
        },
        "view": {"stroke": "transparent"},
    }


CSS = """
<style>
  .block-container { padding-top: 2.2rem; max-width: 1500px; }
  h1, h2, h3 { letter-spacing: -0.01em; }
  .apt-caveat {
      font-size: 0.86rem; color: var(--text-color); opacity: 0.72;
      border-left: 3px solid #2a78d6; padding: 0.35rem 0 0.35rem 0.7rem;
      margin: 0.4rem 0 0.9rem 0;
  }
  .apt-chip {
      display: inline-block; padding: 0.12rem 0.55rem; margin: 0.15rem 0.3rem 0.15rem 0;
      border-radius: 999px; font-size: 0.78rem; line-height: 1.5;
      border: 1px solid rgba(128,128,128,0.35);
  }
  .apt-chip-caution { border-color: #d03b3b; }
  .apt-chip-positive { border-color: #2a78d6; }
  .apt-seq {
      font-family: ui-monospace, "Cascadia Mono", Consolas, monospace;
      font-size: 0.78rem; line-height: 1.45; word-break: break-all;
      letter-spacing: 0.02em;
  }
  .apt-dev-banner {
      background: rgba(208,59,59,0.10); border: 1px solid #d03b3b;
      border-radius: 6px; padding: 0.55rem 0.8rem; margin-bottom: 0.9rem;
      font-size: 0.88rem;
  }
  div[data-testid="stMetricValue"] { font-size: 1.35rem; }
</style>
"""
