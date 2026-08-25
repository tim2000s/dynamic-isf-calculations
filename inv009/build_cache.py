"""Build the per-subject window cache.

One parquet of windows per subject. Grids are not kept: they are a few seconds
to rebuild and would be four gigabytes to store, while the windows are the thing
every analysis reads and are small.

Runs a pool of seven, which is half the cores, so the machine stays usable while
this runs. Everything a worker needs travels in the job tuple, because macOS
spawns rather than forks and a worker re-imports this module rather than
inheriting anything set up in main.

    python3 -m inv009.build_cache --probe            # one subject, report memory
    python3 -m inv009.build_cache --study Loop
    python3 -m inv009.build_cache                    # everything
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import resource
import sys
import time
import traceback

import numpy as np
import pandas as pd

from . import config, db, grid as gridmod, tdd as tddmod, windows as windowmod


def analyse(job) -> tuple[str, int, str]:
    """One subject, start to finish. Returns (subject_id, n_windows, note)."""
    subject_id, study, age = job
    try:
        out = config.WINDOW_CACHE / f"{subject_id.replace(':', '_')}.parquet"
        streams = db.streams(subject_id)
        g = gridmod.build_grid(streams)
        if g is None:
            return subject_id, 0, "no usable span"
        subj = tddmod.subject_level(g)
        tw = tddmod.windowed(g)
        models = config.CACHED_MODELS if study == "Loop" else ("oref_6h75",)
        w = windowmod.extract_windows(g, tw, subject_id, study, age, subj, models)
        if w.empty:
            return subject_id, 0, "no candidate windows"
        w.to_parquet(out, index=False)
        return subject_id, len(w), ""
    except Exception:
        return subject_id, 0, traceback.format_exc().splitlines()[-1]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--study", default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--probe", action="store_true",
                    help="run one subject in this process and report peak memory")
    ap.add_argument("--workers", type=int, default=config.WORKERS)
    args = ap.parse_args()
    config.ensure_dirs()

    subs = db.subjects(args.study)
    subs = subs[subs.study_name.isin(config.COHORTS)]
    if args.limit:
        subs = subs.head(args.limit)
    jobs = [(r.subject_id, r.study_name, r.age_years) for r in subs.itertuples()]

    if args.probe:
        # Measure one worker before launching seven. Twelve workers each holding
        # a full history is how this machine got wedged in INV-008.
        job = next(j for j in jobs if j[1] == "Loop")
        t = time.time()
        sid, n, note = analyse(job)
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
        print(f"probe {sid}: {n} windows in {time.time() - t:.1f}s, peak RSS {rss:.2f} GB")
        print(f"  {args.workers} workers would need about {rss * args.workers:.1f} GB")
        return 0

    print(f"{len(jobs)} subjects, {args.workers} workers")
    t0 = time.time()
    done = total = 0
    notes: dict[str, int] = {}
    with mp.Pool(args.workers, maxtasksperchild=8) as pool:
        for sid, n, note in pool.imap_unordered(analyse, jobs, chunksize=4):
            done += 1
            total += n
            if note:
                notes[note] = notes.get(note, 0) + 1
            if done % 100 == 0 or done == len(jobs):
                print(f"  {done}/{len(jobs)}  {total:,} windows  {time.time() - t0:.0f}s",
                      flush=True)
    print(f"done: {total:,} windows from {done} subjects in {time.time() - t0:.0f}s")
    for note, k in sorted(notes.items(), key=lambda kv: -kv[1]):
        print(f"  {k:4d} x {note}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
