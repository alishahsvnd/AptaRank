"""The dashboard must actually render.

Everything else here tests the functions the dashboard calls. This runs the app
itself, because the failure mode these panels have is a `KeyError` on an
artifact field that moved — which no unit test of a helper would catch, and
which a user would meet as a stack trace instead of a result.
"""

from __future__ import annotations

import sys

import pytest

pytest.importorskip("streamlit", reason="the dashboard is an optional extra")
from streamlit.testing.v1 import AppTest  # noqa: E402

from .conftest import make_surface_bundle, make_synthetic_bundle  # noqa: E402

REPO_ROOT = __import__("pathlib").Path(__file__).resolve().parent.parent
APP = REPO_ROOT / "dashboard" / "streamlit_app.py"


def run_app(runs_dir, view: str | None = None, artifact: str | None = None) -> AppTest:
    """Start the app, then navigate it the way a user would.

    The view is selected through the sidebar radio rather than by writing
    session state, because a value written before the first run is not applied —
    which quietly turned an earlier version of these tests into assertions about
    the sidebar.
    """
    # The app reads --runs-dir from sys.argv, and AppTest runs the script in
    # this process. Setting an `args` attribute on the AppTest does nothing —
    # it just adds an attribute, and the app quietly falls back to the default
    # runs directory, which is how an earlier version of this file managed to
    # assert against an empty Results page.
    original = sys.argv
    sys.argv = ["streamlit_app.py", "--runs-dir", str(runs_dir)]
    try:
        app = AppTest.from_file(str(APP), default_timeout=120)
        app.run()
        assert not app.exception, app.exception

        if artifact:
            app.session_state["artifact_path"] = artifact
        if view:
            app.sidebar.radio[0].set_value(view).run()
        elif artifact:
            app.run()
        return app
    finally:
        sys.argv = original


def rendered_text(app: AppTest) -> str:
    parts = [str(m.value) for m in app.markdown]
    parts += [str(c.value) for c in app.caption]
    parts += [str(w.value) for w in app.warning]
    parts += [str(s.value) for s in app.success]
    parts += [str(i.value) for i in app.info]
    return " ".join(parts)


@pytest.fixture(scope="module")
def artifact_run(tmp_path_factory, request):
    """One real surface-mode run, rendered by the real app."""
    from aptarank.config import load_config
    from aptarank.pipeline import run_pipeline
    from aptarank.tier2 import bundle as bundle_mod

    tmp = tmp_path_factory.mktemp("render")
    bundle_path = bundle_mod.write(make_surface_bundle(tmp), tmp / "targets")
    cfg = load_config(
        overrides={
            "corpus": {"path": str(request.path.parent / "fixtures" / "mini_corpus.csv"),
                       "is_placeholder": True, "allow_placeholder": True,
                       "cache_dir": str(tmp / "cache")},
            "tier1": {"n_ensemble_samples": 10, "shuffle": {"n_shuffles": 20},
                      "parallel": {"workers": 1}},
            "tier2": {"enabled": True, "binding_mode": "surface", "n_candidates": 6,
                      "bundle_path": str(bundle_path),
                      "calibration": {"bank_size": 40, "cache_dir": str(tmp / "bank")}},
            "output": {"dir": str(tmp / "runs"), "n_diagrams": 2},
        }
    )
    out = run_pipeline(cfg, request.path.parent / "fixtures" / "mini_candidates.csv")
    runs_dir = tmp / "runs"
    return runs_dir, str(out.path)


def test_the_app_starts(tmp_path):
    app = run_app(tmp_path)
    assert not app.exception


def test_the_results_view_renders_a_surface_mode_run(artifact_run):
    runs_dir, artifact = artifact_run
    app = run_app(runs_dir, view="Results", artifact=artifact)
    assert not app.exception

    text = rendered_text(app)
    # The results screen itself, not just the sidebar: the target panel and the
    # scatter heading are only drawn in this view.
    assert "Surface-patch recognition" in text
    assert "Binding-site area" in text or "binding-site residues" in text
    # The user-facing vocabulary, on the screen the user reads (§1.1).
    assert "aptamer-likeness" in text.lower()
    assert "compatibility" in text.lower()


def test_the_sidebar_states_what_tier_2_is_and_is_not(tmp_path):
    app = run_app(tmp_path)
    captions = " ".join(str(c.value) for c in app.caption)
    assert "aptamer-likeness" in captions
    assert "aptamer-target compatibility" in captions
    assert "not a binding prediction" in captions


def test_a_development_run_is_labelled_before_the_results_are_shown(artifact_run):
    """The placeholder corpus makes this run ineligible; the page must say so."""
    runs_dir, artifact = artifact_run
    app = run_app(runs_dir, view="Results", artifact=artifact)
    markdown = rendered_text(app)
    assert "Development run" in markdown
    assert "publication_eligible" in markdown
    assert "synthetic example data" in markdown


def test_the_new_analysis_view_renders(tmp_path):
    app = run_app(tmp_path, view="New analysis")
    assert not app.exception
    body = rendered_text(app)
    assert "Protein target" in body
    assert "Reference library" in body
    # The renamed section and its new description (§1.3).
    assert "binding site" in body
    assert "does not predict binding" in body


def test_the_recent_and_progress_views_render(tmp_path):
    for view in ("Recent analyses", "Progress"):
        app = run_app(tmp_path, view=view)
        assert not app.exception, f"{view} raised {app.exception}"


def test_a_pocket_mode_artifact_also_renders(tmp_path_factory, request):
    """Both modes go through the same panels; neither may KeyError on the other."""
    from aptarank.config import load_config
    from aptarank.pipeline import run_pipeline
    from aptarank.tier2 import bundle as bundle_mod

    tmp = tmp_path_factory.mktemp("render_pocket")
    bundle_path = bundle_mod.write(make_synthetic_bundle(tmp), tmp / "targets")
    cfg = load_config(
        overrides={
            "corpus": {"path": str(request.path.parent / "fixtures" / "mini_corpus.csv"),
                       "is_placeholder": True, "allow_placeholder": True,
                       "cache_dir": str(tmp / "cache")},
            "tier1": {"n_ensemble_samples": 10, "shuffle": {"n_shuffles": 20},
                      "parallel": {"workers": 1}},
            "tier2": {"enabled": True, "n_candidates": 6,
                      "bundle_path": str(bundle_path),
                      "calibration": {"bank_size": 40, "cache_dir": str(tmp / "bank")}},
            "output": {"dir": str(tmp / "runs"), "n_diagrams": 2},
        }
    )
    out = run_pipeline(cfg, request.path.parent / "fixtures" / "mini_candidates.csv")
    app = run_app(tmp / "runs", view="Results", artifact=str(out.path))
    assert not app.exception
