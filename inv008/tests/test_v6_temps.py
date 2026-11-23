"""The AndroidAPS export temp-basal reader (inv008.sources.load_v6_temps)."""

import json

from inv008.sources import load_v6_temps


def test_both_layouts_duplicates_and_deletions(tmp_path):
    v1 = {
        "date": 1_000_000,
        "durationInMinutes": 30,
        "isAbsolute": False,
        "percentRate": 70,
        "absoluteRate": 0,
        "isValid": True,
        "isDeletion": False,
        "queuedOn": 1,
    }
    later = {**v1, "percentRate": 50, "queuedOn": 2}  # same start, queued later: wins
    gone = {**v1, "date": 2_000_000, "isDeletion": True}
    v2 = {
        "timestamp": 3_000_000,
        "duration": 600_000,
        "isAbsolute": True,
        "rate": 0.4,
        "isValid": "TemporaryBasal(id=1, isValid=true, rate=0.4)",
        "dateCreated": 5,
    }
    v2_bad = {
        **v2,
        "timestamp": 4_000_000,
        "isValid": "TemporaryBasal(id=2, isValid=false)",
    }
    for i, recs in enumerate(([v1, gone], [later, v2, v2_bad], [v1])):
        d = tmp_path / f"upload-num{i}"
        d.mkdir()
        (d / "TemporaryBasals.json").write_text(json.dumps(recs))
    got = load_v6_temps(tmp_path)
    assert [(e.ts, e.duration_min, e.rate_u_h, e.percent) for e in got] == [
        (1000.0, 30.0, None, 50.0),
        (3000.0, 10.0, 0.4, None),
    ]


def test_repeated_boluses_count_once():
    from inv008.sources import _bolus_arrays

    ts, u = _bolus_arrays([20.0, 10.0, 20.0, 20.0], [1.0, 0.5, 1.0, 0.3])
    assert ts.tolist() == [10.0, 20.0, 20.0] and u.tolist() == [0.5, 0.3, 1.0]
    ts, u = _bolus_arrays([], [])
    assert len(ts) == 0 and len(u) == 0


def test_later_deletion_and_history_rows(tmp_path):
    kept = {
        "date": 1_000_000,
        "durationInMinutes": 30,
        "isAbsolute": True,
        "absoluteRate": 1.2,
        "percentRate": 0,
        "isValid": True,
        "isDeletion": False,
        "queuedOn": 1,
    }
    deleted_later = {**kept, "isDeletion": True, "queuedOn": 2}
    history = {
        "timestamp": 5_000_000,
        "duration": 600_000,
        "isAbsolute": True,
        "rate": 3.0,
        "isValid": "TemporaryBasal(isValid=true)",
        "dateCreated": 9,
        "referenceId": 7,
    }
    current = {**history, "rate": 0.5, "referenceId": None, "dateCreated": 8}
    for i, recs in enumerate(([kept], [deleted_later, history, current])):
        d = tmp_path / f"upload-num{i}"
        d.mkdir()
        (d / "TemporaryBasals.json").write_text(json.dumps(recs))
    assert [(e.ts, e.rate_u_h) for e in load_v6_temps(tmp_path)] == [(5000.0, 0.5)]


def test_reader_passes_temps_through(tmp_path, monkeypatch):
    from inv008 import config, sources

    up = tmp_path / "raw" / "direct-sharing-396" / "upload-num1"
    up.mkdir(parents=True)
    (up / "Treatments.json").write_text(
        json.dumps([{"date": 1000, "insulin": 1.0}, {"date": 1000, "insulin": 1.0}])
    )
    (up / "TemporaryBasals.json").write_text(
        json.dumps(
            [
                {
                    "date": 2000,
                    "durationInMinutes": 30,
                    "isAbsolute": False,
                    "percentRate": 50,
                    "absoluteRate": 0,
                }
            ]
        )
    )
    monkeypatch.setattr(config, "NS_SAMPLES", tmp_path)
    _, u, temps = sources.load_v6_events("raw")
    assert u.tolist() == [1.0] and [e.percent for e in temps] == [50.0]


def test_timestamps_with_milliseconds_and_offset_parse():
    from inv008.sources import parse_ts

    assert parse_ts("2017-02-14T13:06:39.000-06:00") == parse_ts("2017-02-14T19:06:39Z")
    assert parse_ts("2017-02-14T19:06:39") == parse_ts("2017-02-14T19:06:39+00:00")
    assert parse_ts("not a time") is None


def test_v7_reader_counts_a_bolus_in_two_files_once(tmp_path, monkeypatch):
    from inv008 import config, sources

    ds = tmp_path / "raw" / "direct-sharing-31"
    ds.mkdir(parents=True)
    rec = {"created_at": "2017-02-14T13:06:39.000-06:00", "insulin": 1.5}
    for i in (1, 2):
        (ds / f"x_treatments_{i}.json").write_text(json.dumps([rec]))
    monkeypatch.setattr(config, "NS_SAMPLES", tmp_path)
    ts, u, _ = sources.load_v7_events("raw")
    assert u.tolist() == [1.5] and ts.tolist() == [
        sources.parse_ts("2017-02-14T19:06:39Z")
    ]
