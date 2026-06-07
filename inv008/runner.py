"""OREF-INV-008 parallel orchestrator.

Usage:
    python -m inv008.runner --stage 1 --platforms v6 v7
    python -m inv008.runner --stage 2 --platforms v5 v6 v7
    python -m inv008.runner --stage all
    python -m inv008.runner --stage 2 --users U073 U074 --workers 4 --force

Design (Mac mini M4 Pro, 10P+4E, 64 GB):
  * one user = one task; multiprocessing.Pool(workers, maxtasksperchild=4)
  * workers open their own Postgres connections / read their own files — the parent
    never holds bulk data
  * resumable: users whose output parquet already exists are skipped unless --force
  * every status line is timestamped; a full log is written to inv008_cache/logs/
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from datetime import datetime

from inv008 import config


def log(fh, msg: str) -> None:
    line = f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {msg}"
    print(line, flush=True)
    fh.write(line + "\n")
    fh.flush()


def _stage1_task(args):
    from inv008 import stage1_tdd
    try:
        return stage1_tdd.run_user(args)
    except Exception as e:  # never kill the pool on one bad user
        return {"user": args[0], "status": "error", "reason": f"{type(e).__name__}: {e}"}


def _stage2_task(args):
    from inv008 import stage2_replay
    try:
        return stage2_replay.run_user(args)
    except Exception as e:
        return {"user": args[0], "status": "error", "reason": f"{type(e).__name__}: {e}"}


def run_stage(stage: int, platforms: tuple[str, ...], users: list[str] | None,
              workers: int, force: bool, fh) -> list[dict]:
    if stage == 1:
        from inv008 import stage1_tdd
        tasks = stage1_tdd.user_list(tuple(p for p in platforms if p in ("v6", "v7")))
        out_dir, task_fn = config.TDD_DIR, _stage1_task
    else:
        from inv008 import stage2_replay
        tasks = stage2_replay.user_list(platforms)
        out_dir, task_fn = config.REPLAY_DIR, _stage2_task

    if users:
        tasks = [t for t in tasks if t[0] in set(users)]
    if not force:
        done = {p.stem for p in out_dir.glob("*.parquet")}
        skipped = [t for t in tasks if t[0] in done]
        tasks = [t for t in tasks if t[0] not in done]
        if skipped:
            log(fh, f"stage {stage}: {len(skipped)} users already done (use --force to redo)")

    log(fh, f"stage {stage}: {len(tasks)} users, {workers} workers")
    if not tasks:
        return []

    results = []
    t0 = time.time()
    ctx = mp.get_context("spawn")  # macOS-safe; no inherited DB handles
    with ctx.Pool(processes=workers, maxtasksperchild=config.MAXTASKSPERCHILD) as pool:
        for i, res in enumerate(pool.imap_unordered(task_fn, tasks), 1):
            results.append(res)
            tag = res["status"].upper()
            extra = res.get("reason") or ", ".join(
                f"{k}={v}" for k, v in res.items() if k not in ("user", "status") and v is not None)
            log(fh, f"[{i}/{len(tasks)}] {res['user']} {tag} {extra}")
    dt = time.time() - t0
    ok = sum(1 for r in results if r["status"] == "ok")
    log(fh, f"stage {stage} finished: {ok}/{len(tasks)} ok in {dt/60:.1f} min")
    return results


def main() -> None:
    ap = argparse.ArgumentParser(description="OREF-INV-008 ISF replay orchestrator")
    ap.add_argument("--stage", choices=["1", "2", "all"], default="all")
    ap.add_argument("--platforms", nargs="+", default=list(config.PLATFORMS),
                    choices=list(config.PLATFORMS))
    ap.add_argument("--users", nargs="+", default=None)
    ap.add_argument("--workers", type=int, default=config.DEFAULT_WORKERS)
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    config.ensure_dirs()
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_path = config.LOG_DIR / f"run_{run_id}.log"
    manifest = {"run_id": run_id, "argv": sys.argv[1:], "source_commit": config.SOURCE_COMMIT,
                "stages": {}}

    with open(log_path, "w") as fh:
        log(fh, f"INV-008 run {run_id} | platforms={args.platforms} workers={args.workers}")
        stages = [1, 2] if args.stage == "all" else [int(args.stage)]
        for s in stages:
            res = run_stage(s, tuple(args.platforms), args.users, args.workers, args.force, fh)
            manifest["stages"][str(s)] = res
        mpath = config.LOG_DIR / f"manifest_{run_id}.json"
        mpath.write_text(json.dumps(manifest, indent=1))
        log(fh, f"manifest → {mpath}")


if __name__ == "__main__":
    main()
