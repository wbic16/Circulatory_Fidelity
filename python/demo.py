#!/usr/bin/env python3
"""
End-to-end demo of the Circulatory Fidelity x phext integration.

  (a) synthetic filtering model with known rho ~ 0.8 -> diagnose ->
      phext -> verdict read back by coordinate alone
  (b) 4x4 precision matrix -> phext (D lane + R lane) -> reconstruct
  (c) windowed series with one planted high-coupling regime -> phext ->
      IC_max found at an addressable window scroll
  (d) all 8 upstream simulations/*.csv -> one phext (one volume per
      dataset) -> exported back and compared byte-for-byte
  (e) worked query on the imported validation data (see
      worked_example()): by coordinate query alone, find where the
      SVF mse_ratio first reaches >= 2 and how the HLM pooling data
      reads at comparable IC — the Dependency Asymmetry, addressed.

Run:  python3 python/demo.py   (from the repo root; outputs under python/output/)
      python3 demo.py          (from inside python/, same outputs)
Pure Python + numpy/scipy via the upstream module. No SQ required.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cf_phext as px  # noqa: E402

OUT = HERE / "output"
SIM_DIR = next(
    (c for c in [HERE.parent / "upstream" / "simulations",
                 HERE.parent / "simulations"] if c.is_dir()),
    HERE.parent / "upstream" / "simulations",
)
CSV_FILES = [
    "copula_validation.csv",
    "dsprites_proximal_dominance.csv",
    "hlm_validation.csv",
    "ic_psis_comprehensive_validation.csv",
    "svf_validation.csv",
    "three_layer_validation.csv",
    "threshold_calibration.csv",
    "trigger_experiment.csv",
]


def demo_diagnose():
    print("== (a) diagnose: synthetic filtering model, rho ~ 0.8 ==")
    rng = np.random.default_rng(20261101)
    n = 2000
    x = rng.standard_normal(n)
    y = 0.8 * x + 0.6 * rng.standard_normal(n)
    path = OUT / "diagnose_demo.phext"
    path.unlink(missing_ok=True)
    info = px.diagnose_to_phext(x, y, path, model_type="filtering")
    d = info["diagnostic"]
    verdict = px.read_scroll(path, info["verdict_coord"])
    ic_back = float(px.read_scroll(path, info["ic_coord"]))
    print(f"  upstream diagnose(): IC = {d.ic:.4f} +/- {d.se:.4f}, n = {d.n}")
    print(f"  verdict read back @ {info['verdict_coord']}: {verdict}")
    print(f"  IC read back @ {info['ic_coord']}: {ic_back:.4f}")
    assert verdict == d.recommendation
    assert abs(ic_back - d.ic) < 1e-12
    return path, info


def demo_matrix():
    print("== (b) matrix: 4x4 precision matrix -> D lane + R lane -> back ==")
    rng = np.random.default_rng(7)
    A = rng.standard_normal((4, 4))
    M = A @ A.T + 4.0 * np.eye(4)
    path = OUT / "matrix_demo.phext"
    path.unlink(missing_ok=True)
    vol = px.matrix_to_phext(M, path, kind="precision")
    D, R = px.decomposition_from_phext(path, vol)
    back = px.matrix_from_phext(path, vol)
    err = float(np.max(np.abs(back - M)))
    print(f"  volume {vol}; max |reconstruction error| = {err:.3e}")
    print(f"  ||D||_F = {np.linalg.norm(D):.4f}  ||R||_F = {np.linalg.norm(R):.4f}"
          f"  (MFVI keeps D, discards R)")
    max_r_coord = px.read_scroll(path, px.coord_str(px.addr(vol, px.BOOK_SUMMARY, 1, 1, 4)))
    print(f"  largest |R| entry lives at {max_r_coord}"
          f" = {px.read_scroll(path, max_r_coord)}")
    assert err < 1e-12
    return path


def demo_windowed(path, info):
    print("== (c) windowed: planted high-coupling regime ==")
    cf = px.load_cf()
    rng = np.random.default_rng(99)
    T = 600
    z = rng.standard_normal(T)
    noise = rng.standard_normal(T)
    x = np.empty(T)
    x[:400] = 0.1 * z[:400] + np.sqrt(1 - 0.01) * noise[:400]     # low coupling
    x[400:] = 0.9 * z[400:] + np.sqrt(1 - 0.81) * noise[400:]     # high coupling
    result = cf.windowed_ic(z, x, window_size=50)
    vol, max_coord = px.windowed_to_phext(result, path)  # appends a volume
    ic_at_max = float(px.read_scroll(path, max_coord))
    verdict_coord = px.coord_str(px.addr(vol, px.BOOK_SUMMARY, 1, 1, 6))
    print(f"  appended as volume {vol} of {path.name}")
    print(f"  IC_max = {result['ic_max']:.4f} across {result['n_windows']} windows")
    print(f"  IC_max lives at {max_coord} -> read back {ic_at_max:.4f}")
    print(f"  verdict @ {verdict_coord}: {px.read_scroll(path, verdict_coord)}")
    assert abs(ic_at_max - float(result["ic_max"])) < 1e-12
    return vol, max_coord, result


def demo_csv_roundtrip():
    print("== (d) CSV: import all 8 simulation datasets, export back ==")
    path = OUT / "simulations.phext"
    path.unlink(missing_ok=True)
    exp_dir = OUT / "exported"
    exp_dir.mkdir(parents=True, exist_ok=True)
    all_ok = True
    for name in CSV_FILES:
        vol = px.import_csv_to_phext(SIM_DIR / name, path)
        text = px.export_phext_to_csv(path, vol, exp_dir / name)
        orig = (SIM_DIR / name).read_text(encoding="utf-8")
        ok = text == orig
        all_ok &= ok
        print(f"  volume {vol}: {name}: {'byte-exact' if ok else 'MISMATCH'}")
    assert all_ok
    return path


def _dataset_volume(doc, name):
    for v in doc.volumes():
        if doc.volume_kind(v) == "dataset" and doc.get_a(v, 1, 1, 1, 3) == name:
            return v
    raise KeyError(name)


def _column(doc, volume, col_idx, n_rows):
    return [doc.get_a(volume, px.BOOK_DATASET, r + 2, 1, col_idx + 1) for r in range(n_rows)]


def worked_example(path):
    """A question the phext answers by coordinate query alone.

    SVF (filtering): as IC grows, mse_ratio grows — MFVI fails.
    HLM (pooling): the same IC axis means the opposite thing — at high
    IC the groups are well separated and pooling matters least.
    We find, by scanning dataset-lane scrolls, the first SVF row whose
    mse_ratio >= 2 and report its address; then the mean mse_ratio of
    HLM rows in the top IC decile vs the bottom decile.
    """
    print("== (e) worked example: Dependency Asymmetry by coordinate query ==")
    doc = px.PhextDoc.load(path)
    # SVF: header coupling,ic,mf_mse,oracle_mse,mse_ratio,rep
    # Aggregate by the coupling column (400 reps per level): single
    # replicates are noisy, level means are the honest summary.
    v = _dataset_volume(doc, "svf_validation.csv")
    n_rows = int(doc.get_a(v, 1, 2, 1, 1))
    coupling = [float(t) for t in _column(doc, v, 0, n_rows)]
    mse = [float(t) for t in _column(doc, v, 4, n_rows)]
    levels = sorted(set(coupling))
    means = {c: float(np.mean([m for cc, m in zip(coupling, mse) if cc == c]))
             for c in levels}
    cross = next(c for c in levels if means[c] >= 2.0)
    first_row = coupling.index(cross)
    row_coord = px.coord_str(px.addr(v, px.BOOK_DATASET, first_row + 2, 1, 5))
    print(f"  SVF volume {v}: mean mse_ratio by coupling level rises from "
          f"{means[levels[0]]:.3f} (coupling {levels[0]}) to "
          f"{means[levels[-1]]:.2f} (coupling {levels[-1]})")
    print(f"    it first reaches >= 2 at coupling = {cross:.3f} "
          f"(mean mse_ratio = {means[cross]:.3f} over 400 reps);")
    print(f"    that level's first replicate mse_ratio scroll is at {row_coord}"
          f" = {mse[first_row]:.3f}")
    first = first_row
    # HLM: header tau,icc,reliability,ic,no_pool_mse,partial_pool_mse,mse_ratio,rep
    vh = _dataset_volume(doc, "hlm_validation.csv")
    nh = int(doc.get_a(vh, 1, 2, 1, 1))
    ic_h = np.array([float(t) for t in _column(doc, vh, 3, nh)])
    mse_h = np.array([float(t) for t in _column(doc, vh, 6, nh)])
    lo = mse_h[ic_h <= np.quantile(ic_h, 0.1)].mean()
    hi = mse_h[ic_h >= np.quantile(ic_h, 0.9)].mean()
    print(f"  HLM volume {vh}: mean no-pool/partial-pool mse_ratio,"
          f" bottom IC decile = {lo:.2f}, top IC decile = {hi:.2f}")
    print("    filtering reads high IC as failure; pooling reads the same")
    print("    axis in reverse — both answers fetched by coordinate, no")
    print("    matrix was parsed to get them.")
    return {"svf_cross_coupling": cross, "svf_coord": row_coord,
            "hlm_lo": float(lo), "hlm_hi": float(hi)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    path_a, info = demo_diagnose()
    demo_matrix()
    demo_windowed(path_a, info)
    path_d = demo_csv_roundtrip()
    worked_example(path_d)
    print("\nDemo complete: all sections passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
