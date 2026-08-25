"""Reading the window cache back, with the screen applied.

Every analysis starts here, so the screening decision is made once and in one
place. Loading is by study because the whole cache is two million rows wide
enough that holding all of it plus a model is wasteful when no analysis needs
more than one cohort at a time.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from . import config, windows as windowmod

ACTION_MODEL_DEFAULT = "oref_6h75"


def _files(study: str | None = None):
    pat = f"{study}_*.parquet" if study else "*.parquet"
    return sorted(config.WINDOW_CACHE.glob(pat))


def load(study: str | None = None, screened: bool = True, strict: bool = False,
         model: str | None = None, columns: list[str] | None = None) -> pd.DataFrame:
    """Windows for one study or all of them.

    `model` names the insulin model whose action columns become the plain `a_pre`,
    `a_in` and their routine-adjusted twins. Loop subjects carry several; everyone
    else carries the oref model only, so asking for a Loop preset outside Loop
    falls back rather than failing.
    """
    frames = []
    for f in _files(study):
        d = pd.read_parquet(f, columns=columns)
        if d.empty:
            continue
        frames.append(d)
    if not frames:
        return pd.DataFrame()
    w = pd.concat(frames, ignore_index=True)
    if screened:
        w = w[windowmod.screen(w, strict=strict)].copy()
    return attach_action(w, model)


def attach_action(w: pd.DataFrame, model: str | None = None) -> pd.DataFrame:
    """Give the chosen insulin model's action columns their plain names."""
    if w.empty:
        return w
    model = model or ACTION_MODEL_DEFAULT
    if f"a_pre_{model}" not in w.columns:
        model = ACTION_MODEL_DEFAULT
    for pref in ("a_pre", "a_in", "a_net_pre", "a_net_in", "iob0"):
        w[pref] = w[f"{pref}_{model}"]
    w["a_tot"] = w.a_pre + w.a_in
    w["a_net"] = w.a_net_pre + w.a_net_in
    w["action_model"] = model
    return w


def per_subject_models(choice: pd.DataFrame | None = None) -> dict[str, str]:
    """Which insulin model each Loop subject was judged to be running."""
    if choice is None:
        f = config.RESULTS / "inv009_loop_model_choice.parquet"
        if not f.exists():
            return {}
        choice = pd.read_parquet(f)
    return dict(zip(choice.subject_id, choice.model))


def cohort_summary(strict: bool = False) -> pd.DataFrame:
    """How much usable data each cohort actually contributes."""
    rows = []
    for study in config.COHORTS:
        w = load(study, screened=False)
        if w.empty:
            continue
        ok = windowmod.screen(w, strict=strict)
        per = w[ok].groupby("subject_id").size()
        rows.append(dict(
            study=study,
            label=config.COHORTS[study]["label"],
            closed_loop=config.COHORTS[study]["closed_loop"],
            logs_carbs=config.COHORTS[study]["carbs"],
            subjects=w.subject_id.nunique(),
            candidate_windows=len(w),
            screened_windows=int(ok.sum()),
            subjects_with_enough=int((per >= config.MIN_WINDOWS).sum()),
            median_windows=float(per.median()) if len(per) else np.nan,
            median_tdd=float(w.groupby("subject_id").tdd_u.first().median()),
            median_age=float(w.groupby("subject_id").age.first().median()),
        ))
    return pd.DataFrame(rows)
