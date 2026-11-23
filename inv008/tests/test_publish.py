"""Per-person results are removed before publishing (inv008.publish)."""

from inv008.publish import strip_json, strip_md


def test_json_person_records_dropped():
    assert strip_json([{"user_id": "U900"}]) is None
    assert strip_json(
        {"n": 3, "per_person": [{"user_id": "U900"}], "per_user": {}}
    ) == {"n": 3}
    assert strip_json({"n": 3, "mixed_unit_users": ["U901", "U902"]}) == {"n": 3}


def test_markdown_person_rows_and_empty_tables_removed():
    md = "\n".join(
        [
            "# T",
            "",
            "| a | b |",
            "|---|---|",
            "| all | 1 |",
            "",
            "| user | x |",
            "|---|---|",
            "| U900 | 2 |",
            "| U901 | 3 |",
            "",
            "end",
        ]
    )
    assert strip_md(md) == "\n".join(
        ["# T", "", "| a | b |", "|---|---|", "| all | 1 |", "", "", "end"]
    )


def test_publish_refuses_an_id_the_filters_leave(tmp_path, monkeypatch):
    import pytest

    import inv008.publish as pub

    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "x.md").write_text(
        "Two people switched units: U900 and U901.\n"
    )
    monkeypatch.setattr(pub.config, "ROOT", tmp_path)
    monkeypatch.setattr(pub, "REPO", tmp_path / "repo")
    with pytest.raises(SystemExit):
        pub.publish("x.md")
