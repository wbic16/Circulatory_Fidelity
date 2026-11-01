#!/usr/bin/env python3
"""
Tests for cf_phext. Runnable two ways:

    python3 test_cf_phext.py     (plain runner, prints PASS/FAIL per test)
    pytest test_cf_phext.py      (same functions under pytest)

Covers: parser agreement with the reference phext parser, matrix D+R
reconstruction, diagnose round-trip, CSV round-trip for all 8 upstream
simulation files (field-text exact, incl. quoted fields containing
commas), whole-file UTF-8 validity, and an SQ host/select round-trip
(SKIP with a printed reason only if SQ genuinely cannot run).
"""

from __future__ import annotations

import importlib.util
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cf_phext as px  # noqa: E402

def _find_sim_dir():
    """Simulations dir in either layout: integration bundle
    (../upstream/simulations) or upstream repo (../simulations,
    when this file sits in python/)."""
    for cand in [HERE.parent / "upstream" / "simulations",
                 HERE.parent / "simulations",
                 HERE / "simulations"]:
        if cand.is_dir():
            return cand
    return HERE.parent / "upstream" / "simulations"


SIM_DIR = _find_sim_dir()
CSV_FILES = sorted(p.name for p in SIM_DIR.glob("*.csv")) if SIM_DIR.is_dir() else []
TMP = HERE / "output" / "test_tmp"


def _tmp(name):
    TMP.mkdir(parents=True, exist_ok=True)
    p = TMP / name
    p.unlink(missing_ok=True)
    return p


def _reference_parser():
    """Locate parse_phext.py (the reference dead-reckoning parser)."""
    for cand in [
        HERE.parent.parent / "phext" / "parse_phext.py",
        Path.home() / "workspace" / "phext" / "parse_phext.py",
    ]:
        if cand.is_file():
            spec = importlib.util.spec_from_file_location("parse_phext", cand)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod
    return None


def test_parser_agreement():
    ref = _reference_parser()
    if ref is None:
        return "SKIP: reference parse_phext.py not found"
    path = _tmp("agreement.phext")
    rng = np.random.default_rng(3)
    A = rng.standard_normal((3, 3))
    px.matrix_to_phext(A @ A.T + 3 * np.eye(3), path, kind="precision")
    ref_scrolls = dict()
    for coord_s, text in ref.parse(str(path)):
        ref_scrolls[px.parse_coord(coord_s)] = text
    mine = px.parse_file(path)
    assert ref_scrolls == mine, (
        f"parser disagreement: ref {len(ref_scrolls)} scrolls vs mine {len(mine)}; "
        f"first diff: {next(((c, ref_scrolls.get(c), mine.get(c)) for c in set(ref_scrolls) | set(mine) if ref_scrolls.get(c) != mine.get(c)), None)}"
    )


def test_matrix_reconstruction():
    rng = np.random.default_rng(42)
    A = rng.standard_normal((6, 6))
    M = A @ A.T + 6 * np.eye(6)
    for kind, Mat in (("precision", M), ("correlation", np.corrcoef(rng.standard_normal((6, 50))))):
        path = _tmp(f"matrix_{kind}.phext")
        vol = px.matrix_to_phext(Mat, path, kind=kind)
        D, R = px.decomposition_from_phext(path, vol)
        assert np.max(np.abs(D - np.diag(np.diag(Mat)))) < 1e-12
        assert np.max(np.abs(R - (Mat - np.diag(np.diag(Mat))))) < 1e-12
        back = px.matrix_from_phext(path, vol)
        err = float(np.max(np.abs(back - Mat)))
        assert err < 1e-12, f"{kind}: reconstruction error {err}"


def test_diagnose_roundtrip():
    cf = px.load_cf()
    rng = np.random.default_rng(2026)
    n = 1500
    x = rng.standard_normal(n)
    y = 0.8 * x + 0.6 * rng.standard_normal(n)
    expected = cf.diagnose(x, y, model_type="filtering")
    path = _tmp("diagnose.phext")
    info = px.diagnose_to_phext(x, y, path, model_type="filtering")
    v = info["volume"]
    ic = float(px.read_scroll(path, px.coord_str(px.addr(v, 2, 1, 1, 1))))
    se = float(px.read_scroll(path, px.coord_str(px.addr(v, 2, 1, 1, 2))))
    nn = int(px.read_scroll(path, px.coord_str(px.addr(v, 2, 1, 1, 4))))
    verdict = px.read_scroll(path, px.coord_str(px.addr(v, 2, 1, 1, 6)))
    assert abs(ic - expected.ic) < 1e-12
    assert abs(se - expected.se) < 1e-12
    assert nn == expected.n
    assert verdict == expected.recommendation
    # The lanes must hold the precision decomposition of the pair.
    D, R = px.decomposition_from_phext(path, v)
    assert np.max(np.abs((D + R) - info["precision"])) < 1e-12


def test_windowed_max_addressable():
    cf = px.load_cf()
    rng = np.random.default_rng(5)
    T = 400
    z = rng.standard_normal(T)
    noise = rng.standard_normal(T)
    x = 0.15 * z + np.sqrt(1 - 0.15 ** 2) * noise
    x[300:] = 0.9 * z[300:] + np.sqrt(1 - 0.81) * noise[300:]
    result = cf.windowed_ic(z, x, window_size=50)
    path = _tmp("windowed.phext")
    vol, max_coord = px.windowed_to_phext(result, path)
    assert abs(float(px.read_scroll(path, max_coord)) - float(result["ic_max"])) < 1e-12
    recorded = px.read_scroll(path, px.coord_str(px.addr(vol, 2, 1, 1, 8)))
    assert recorded == max_coord


def _norm_trailing_newline(s: str) -> str:
    return s[:-1] if s.endswith("\n") else s


def test_csv_roundtrip_all():
    assert len(CSV_FILES) == 8, f"expected 8 simulation CSVs, found {CSV_FILES}"
    path = _tmp("all_csv.phext")
    volumes = {}
    for name in CSV_FILES:
        volumes[name] = px.import_csv_to_phext(SIM_DIR / name, path)
    for name in CSV_FILES:
        exported = px.export_phext_to_csv(path, volumes[name])
        original = (SIM_DIR / name).read_text(encoding="utf-8")
        assert _norm_trailing_newline(exported) == _norm_trailing_newline(original), \
            f"CSV round-trip mismatch: {name}"


def test_csv_quoted_field_with_comma():
    """threshold_calibration.csv's ic_bin column contains commas inside
    quotes, e.g. "(-0.000437, 0.0437]". Per-field verbatim storage must
    keep that field as ONE scroll and reproduce the quoting on export."""
    import csv as _csv
    import io as _io
    src = SIM_DIR / "threshold_calibration.csv"
    rows = list(_csv.reader(_io.StringIO(src.read_text(encoding="utf-8"), newline="")))
    target = next(r[3] for r in rows[1:] if "," in r[3])
    assert target.startswith("(") and target.endswith("]")
    path = _tmp("quoted.phext")
    vol = px.import_csv_to_phext(src, path)
    row_idx = next(i for i, r in enumerate(rows[1:], start=0) if "," in r[3])
    stored = px.read_scroll(path, px.coord_str(px.addr(vol, 6, row_idx + 2, 1, 4)))
    assert stored == target, f"quoted field mangled: {stored!r} != {target!r}"
    exported = px.export_phext_to_csv(path, vol)
    assert f'"{target}"' in exported, "export did not re-quote the comma field"


def test_utf8_whole_file():
    path = _tmp("utf8.phext")
    rng = np.random.default_rng(11)
    x = rng.standard_normal(500)
    y = 0.5 * x + np.sqrt(0.75) * rng.standard_normal(500)
    px.diagnose_to_phext(x, y, path)
    raw = Path(path).read_bytes()
    raw.decode("utf-8")  # strict: must not raise
    assert px.parse_file(path)[px.parse_coord("1.1.1/1.1.1/1.1.1")] == px.MANIFEST_TEXT


def test_sq_roundtrip():
    sq = shutil.which("sq") or str(Path.home() / ".cargo" / "bin" / "sq")
    if not Path(sq).exists():
        return "SKIP: sq binary not found"
    work = TMP / "sqtest"
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True)
    src = _tmp("sq_src.phext")
    rng = np.random.default_rng(21)
    A = rng.standard_normal((3, 3))
    M = A @ A.T + 3 * np.eye(3)
    vol = px.matrix_to_phext(M, src, kind="precision")
    phext_file = work / "sq_demo.phext"
    shutil.copyfile(src, phext_file)
    # Pick an R-lane coordinate: entry (0,1) value scroll.
    r_coord = px.coord_str(px.addr(vol, px.BOOK_R, 1, 2, 1))
    expected = px.read_scroll(src, r_coord)
    env = dict(os.environ)
    env["PATH"] = str(Path(sq).parent) + os.pathsep + env.get("PATH", "")

    def run(*args, timeout=15):
        return subprocess.run([sq, *args], cwd=work, env=env,
                              capture_output=True, text=True, timeout=timeout)

    share = subprocess.Popen([sq, "share", phext_file.name], cwd=work, env=env,
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:
        up = False
        for _ in range(40):
            try:
                r = run("status", timeout=5)
                if r.returncode == 0 and "Hosting" in r.stdout:
                    up = True
                    break
            except subprocess.TimeoutExpired:
                pass
            time.sleep(0.5)
        if not up:
            return "SKIP: sq daemon did not come up in its temp dir"
        toc_out = run("toc").stdout
        assert "1.1.1/1.1.1/1.1.1" in toc_out, f"sq toc missing manifest: {toc_out[:200]}"
        sel = run("select", r_coord).stdout.strip()
        # sq prints "coord: text"
        got = sel.split(": ", 1)[1] if ": " in sel else sel
        assert got == expected, f"sq select {r_coord}: {got!r} != reader {expected!r}"
    finally:
        try:
            run("shutdown", timeout=10)
        except Exception:
            pass
        try:
            share.wait(timeout=10)
        except Exception:
            share.kill()


ALL_TESTS = [
    test_parser_agreement,
    test_matrix_reconstruction,
    test_diagnose_roundtrip,
    test_windowed_max_addressable,
    test_csv_roundtrip_all,
    test_csv_quoted_field_with_comma,
    test_utf8_whole_file,
    test_sq_roundtrip,
]


def main():
    failures = 0
    skips = 0
    for t in ALL_TESTS:
        try:
            result = t()
            if isinstance(result, str) and result.startswith("SKIP"):
                skips += 1
                print(f"SKIP {t.__name__}: {result[5:]}")
            else:
                print(f"PASS {t.__name__}")
        except Exception as e:  # noqa: BLE001
            failures += 1
            print(f"FAIL {t.__name__}: {type(e).__name__}: {e}")
    passed = len(ALL_TESTS) - failures - skips
    print(f"\n{passed}/{len(ALL_TESTS)} tests passed"
          + (f", {skips} skipped" if skips else ""))
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
