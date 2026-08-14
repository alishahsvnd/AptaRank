"""The dashboard's theme is configuration, and it has to actually be loaded.

Two failures this pins. The first is that Streamlit resolves
`.streamlit/config.toml` relative to the *working directory*, so a service
started from anywhere else silently ignores it and follows the viewer's system
preference instead — which is exactly what the deployment did. The second is
drift between the two places colour is declared: Streamlit's own widgets come
from config.toml, the Altair charts from dashboard/theme.py, and the config
comment promises they are the same palette.
"""

from __future__ import annotations

import tomllib
from pathlib import Path

from dashboard.theme import LIGHT

REPO_ROOT = Path(__file__).resolve().parent.parent
CONFIG = REPO_ROOT / ".streamlit" / "config.toml"


def theme() -> dict:
    return tomllib.loads(CONFIG.read_text(encoding="utf-8"))["theme"]


def test_the_dashboard_defaults_to_light_not_to_the_system_preference():
    assert theme()["base"] == "light"


def test_the_light_theme_is_fully_specified():
    """A partially specified theme lets the browser preference fill the gaps."""
    for key in ("primaryColor", "backgroundColor", "secondaryBackgroundColor",
                "textColor"):
        assert theme().get(key), f"{key} is not set, so the viewer's system decides it"


def test_streamlit_widgets_and_altair_charts_share_one_palette():
    """The config comment promises this; nothing else enforces it."""
    assert theme()["primaryColor"] == LIGHT["series"]
    assert theme()["backgroundColor"] == LIGHT["surface"]
    assert theme()["textColor"] == LIGHT["text_primary"]


def test_the_service_starts_from_the_checkout_so_the_config_is_found():
    """Streamlit reads .streamlit/config.toml from the working directory."""
    script = (REPO_ROOT / "deploy" / "aptarank.sh").read_text(encoding="utf-8")
    start = script.split("start)", 1)[1].split("stop)", 1)[0]
    assert 'cd "$APP_DIR"' in start
    assert start.index('cd "$APP_DIR"') < start.index("streamlit run")


def test_the_local_launcher_does_the_same():
    launcher = (REPO_ROOT / "scripts" / "start.py").read_text(encoding="utf-8")
    assert "cwd=str(REPO_ROOT)" in launcher


# -- user-facing copy ----------------------------------------------------


def test_a_binding_mode_has_one_name_on_every_screen():
    """§1.1: the drift this rule exists to stop is a mode called one thing while
    you choose it and another while you read the result."""
    from dashboard.inputs import (
        BINDING_MODE_DESCRIPTION, BINDING_MODE_LABEL, BINDING_MODE_PREMISE,
    )

    assert BINDING_MODE_LABEL["pocket"] == "Pocket/groove recognition"
    assert BINDING_MODE_LABEL["surface"] == "Surface-patch recognition"
    # Every mode is named, described for the chooser, and described for the reader.
    assert set(BINDING_MODE_LABEL) == set(BINDING_MODE_PREMISE) == set(BINDING_MODE_DESCRIPTION)


def test_the_chooser_and_the_result_page_say_different_amounts():
    """Choosing needs 'which kind of site is this?'; reading a band needs the
    mechanism the comparison assumed. Same name, different depth."""
    from dashboard.inputs import BINDING_MODE_DESCRIPTION, BINDING_MODE_PREMISE

    for mode in ("pocket", "surface"):
        assert len(BINDING_MODE_DESCRIPTION[mode]) > len(BINDING_MODE_PREMISE[mode])
        assert "Common for" in BINDING_MODE_PREMISE[mode]

    assert "van der Waals" in BINDING_MODE_DESCRIPTION["pocket"]
    assert "exosite" in BINDING_MODE_DESCRIPTION["surface"]


def test_the_shuffled_control_explains_what_it_tests():
    from dashboard.views import SHUFFLE_HELP

    assert "letter composition" in SHUFFLE_HELP
    assert "scramble the order" in SHUFFLE_HELP
    # It must describe what passing means, not merely assert the word.
    assert "how the sequence is arranged" in SHUFFLE_HELP


def test_a_library_can_describe_itself_in_its_manifest(tmp_path):
    """So one library's wording never has to be hard-coded into the UI."""
    import json
    from dashboard.inputs import inspect_library

    path = tmp_path / "lib.csv"
    path.write_text(
        "id,sequence,target_name\n"
        + "".join(f"r{i},{'ACGU' * 8},T{i % 7}\n" for i in range(200)),
        encoding="utf-8",
    )
    path.with_suffix(".manifest.json").write_text(
        json.dumps({
            "name": "UTexas Aptamer Database subset",
            "description": "336 experimentally validated RNA aptamers spanning "
                           "multiple target types.",
            "source": "s", "curator": "c", "curated_date": "d",
        }),
        encoding="utf-8",
    )
    library = inspect_library(path)
    assert library.name == "UTexas Aptamer Database subset"
    assert library.describe().startswith("336 experimentally validated RNA aptamers")


def test_the_preselected_library_is_the_one_labelled_default(tmp_path):
    """The marker follows whichever library is actually preselected."""
    import json
    from dashboard.inputs import inspect_library
    from dashboard.newrun import _library_label

    path = tmp_path / "lib.csv"
    path.write_text(
        "id,sequence,target_name\n"
        + "".join(f"r{i},{'ACGU' * 8},T{i % 7}\n" for i in range(200)),
        encoding="utf-8",
    )
    path.with_suffix(".manifest.json").write_text(
        json.dumps({"name": "Some library", "description": "200 aptamers.",
                    "source": "s", "curator": "c", "curated_date": "d"}),
        encoding="utf-8",
    )
    library = inspect_library(path)
    assert "(default)" in _library_label(library, default=True)
    assert "(default)" not in _library_label(library)


def test_the_multi_chain_caveat_says_the_same_thing_before_and_after_the_run():
    """The New Analysis warning and the prepared target's own warning describe
    one fact, so they must not drift apart — and neither may claim the
    coordinates relaxed, because deleting the partner does not move an atom."""
    from pathlib import Path

    from dashboard.inputs import TargetEvidence, TargetRequest, review

    prepared = TargetEvidence(path=Path("x"), pdb_id="7WRQ", chain="B",
                              synthetic=False, binding_mode="surface",
                              was_multi_chain=True)
    verdict = review({"ok": True, "n_submitted": 5, "n_valid": 5, "n_rejected": 0},
                     None, TargetRequest(kind="prepared", prepared=prepared,
                                         label="7WRQ", binding_mode="surface"),
                     preset="standard")
    before = next(w for w in verdict["warnings"] if "multi-chain" in w)

    source = (REPO_ROOT / "src/aptarank/tier2/target.py").read_text(encoding="utf-8")
    after = source.split("was_multi_chain\"]:", 1)[1].split("return PreparedTarget", 1)[0]

    for text in (before, after):
        assert "bound-state" in text
        assert "unbound" in text          # only ever to say it is *not* that
        assert "relax" in text
    assert "computational removal" in before
    assert "isolated computationally" in after
