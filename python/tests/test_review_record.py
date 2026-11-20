"""Every review id the repository cites resolves in the phase 3 review record (council round 10, Mg4).

Ids are cited as "r{N} {Id}" or "council round N, {Id}" in the plan, the CHANGELOG, the README and the
package source; each must be a row id (findings or limits) in round N's section of
`docs/reviews/2026-11-phase3-council.md`. The verbatim rulings and reports under `docs/reviews/phase3/`
carry no home or scratch paths.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
SUMMARY = ROOT / "docs" / "reviews" / "2026-11-phase3-council.md"
RECORD = ROOT / "docs" / "reviews" / "phase3"

pytestmark = pytest.mark.skipif(not SUMMARY.exists(), reason="repository docs/ not present")

CITE_SHORT = re.compile(r"\br(\d+) ([A-Z][A-Za-z]*\d+)")
CITE_LONG = re.compile(r"council round (\d+)[^)]*?, ([A-Z][A-Za-z]*\d+)", re.IGNORECASE)


def _cited_in() -> list[Path]:
    return [ROOT / "docs" / "phase3-plan.md", ROOT / "CHANGELOG.md", ROOT / "README.md",
            *sorted((ROOT / "python" / "dadl").rglob("*.py"))]


def citations() -> set[tuple[int, str, str]]:
    """(round, id, file) for every citation."""
    out = set()
    for path in _cited_in():
        text = path.read_text(encoding="utf-8")
        for rx in (CITE_SHORT, CITE_LONG):
            for m in rx.finditer(text):
                out.add((int(m.group(1)), m.group(2), str(path.relative_to(ROOT))))
    return out


def row_ids() -> dict[int, set[str]]:
    """Round -> the ids in the first column of every table row of its section ("S6 + D8" gives both)."""
    ids: dict[int, set[str]] = {}
    current = None
    for line in SUMMARY.read_text(encoding="utf-8").splitlines():
        head = re.match(r"^## Round (\d+)\b", line)
        if head:
            current = int(head.group(1))
            ids[current] = set()
        elif line.startswith("## "):
            current = None
        elif current is not None and line.startswith("| ") and not line.startswith("| ---"):
            first = line.split("|")[1]
            ids[current].update(re.findall(r"[A-Z][A-Za-z]*\d+", first))
    return ids


def test_the_summary_has_ten_rounds():
    assert sorted(row_ids()) == list(range(1, 11))


def test_every_cited_id_resolves_in_its_round():
    ids = row_ids()
    missing = sorted(f"{file}: r{n} {i}" for n, i, file in citations() if i not in ids.get(n, set()))
    assert not missing, "unresolved review ids:\n" + "\n".join(missing)


def test_the_resolver_finds_citations_and_rejects_an_unknown_id():
    # The test of the test: the two forms are found, and an id absent from a round is not resolved.
    found = {(n, i) for n, i, _ in citations()}
    assert (8, "Pl6") in found and (9, "Rd1") in found and (10, "Rg1") in found
    ids = row_ids()
    assert "S6" in ids[1] and "D8" in ids[1]
    assert "Rg9" not in ids[10] and "Rd1" not in ids[8]
    assert CITE_LONG.search("(council round 9, Rd1b: an echo)").group(2) == "Rd1"


def test_the_record_has_no_home_or_scratch_paths():
    files = sorted(RECORD.glob("*")) + [SUMMARY]
    assert len(files) == 41 and all(f.suffix == ".md" for f in files)
    for f in files:
        text = f.read_text(encoding="utf-8")
        assert "/Users/" not in text and "/private/tmp" not in text, f.name
