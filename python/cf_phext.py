#!/usr/bin/env python3
"""
cf_phext — a phext integration for Circulatory Fidelity (CF).

Circulatory Fidelity splits a precision matrix by *position*:

    Λ = D + R      D: diagonal entries (i == j), retained by MFVI
                   R: off-diagonal entries (i != j), discarded by MFVI

Phext gives structure an *address*: a text stream cut by ASCII control
delimiters into scrolls at 9-dimensional coordinates

    Library.Shelf.Series / Collection.Volume.Book / Chapter.Section.Scroll

This module stores CF diagnostics in phext so that the coordinate space
*is* the coupling graph: the nodal part D and the relational part R live
in distinct Book lanes, every matrix entry is one record of per-field
scrolls, and a single relational entry can be fetched by coordinate
without parsing a matrix. CSV datasets are stored field-verbatim, one
scroll per field, so import/export round-trips are text-exact.

Coordinate scheme (also written into every phext as the manifest scroll
at 1.1.1/1.1.1/1.1.1 — see MANIFEST_TEXT):

    Library.Shelf.Series/Collection = 1.1.1/1 (fixed for study phexts)
    Volume  = one diagnostic run / one matrix / one windowed result /
              one imported dataset, incrementing per volume
    Book    = artifact-class lane:
                1 manifest/meta   2 summary+verdict   3 nodal D
                4 relational R    5 windowed series   6 dataset rows
                7 provenance
    Chapter.Section.Scroll = the record address inside a lane
              (per-lane field orders are documented in the manifest)

Dependencies: Python standard library + numpy. The upstream
circulatory_fidelity module (numpy/scipy) is imported by file path only
by the diagnose/windowed helpers; the phext reader/writer itself needs
nothing beyond the standard library.

License: MIT (same as upstream Circulatory Fidelity).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

__version__ = "1.1.0"

# Phext storage schema version written by this module. Schema 1 files
# (written by cf_phext 1.0.0) have matrix-entry fields 1-6 only and no
# Book 7 Scroll 6; schema 2 adds fields 7-9 (conditional coupling) and
# stamps the version at Book 7 Scroll 6 of every volume. Readers treat
# a missing Book 7 Scroll 6 as schema "1" (see coupling_from_phext).
SCHEMA_VERSION = "2"

# ---------------------------------------------------------------------------
# Phext core: delimiters, coordinates, serialize/parse
# ---------------------------------------------------------------------------
# Delimiter hierarchy, highest -> lowest (must agree with the reference
# parser, libphext / phext.io):
#   0x01 Library, 0x1F Shelf, 0x1E Series, 0x1D Collection, 0x1C Volume,
#   0x1A Book, 0x19 Chapter, 0x18 Section, 0x17 Scroll
# Hitting a delimiter increments that level and resets all lower levels
# to 1. Coordinates are 1-indexed; the first scroll is 1.1.1/1.1.1/1.1.1.
LEVEL_BYTES = [0x01, 0x1F, 0x1E, 0x1D, 0x1C, 0x1A, 0x19, 0x18, 0x17]
LEVEL_INDEX = {b: i for i, b in enumerate(LEVEL_BYTES)}
_DELIM_RE = re.compile(rb"[\x01\x1f\x1e\x1d\x1c\x1a\x19\x18\x17]")
_DELIM_CHARS = {chr(b) for b in LEVEL_BYTES}

# Lane (Book) numbers
BOOK_MANIFEST = 1
BOOK_SUMMARY = 2
BOOK_D = 3
BOOK_R = 4
BOOK_WINDOWED = 5
BOOK_DATASET = 6
BOOK_PROVENANCE = 7


def coord_str(coord) -> str:
    """(9 ints) -> '1.1.1/1.1.1/1.1.1'."""
    c = [int(v) for v in coord]
    return f"{c[0]}.{c[1]}.{c[2]}/{c[3]}.{c[4]}.{c[5]}/{c[6]}.{c[7]}.{c[8]}"


def parse_coord(text: str):
    """'1.1.1/1.2.3/4.5.6' -> tuple of 9 ints."""
    parts = text.strip().split("/")
    if len(parts) != 3:
        raise ValueError(f"bad coordinate: {text!r}")
    out = []
    for part in parts:
        nums = part.split(".")
        if len(nums) != 3:
            raise ValueError(f"bad coordinate: {text!r}")
        out.extend(int(n) for n in nums)
    if any(v < 1 for v in out):
        raise ValueError(f"coordinates are 1-indexed: {text!r}")
    return tuple(out)


def addr(volume: int, book: int, chapter: int, section: int, scroll: int):
    """Study-phext address: fixed 1.1.1/1 prefix + (volume, book, chapter.section.scroll)."""
    return (1, 1, 1, 1, int(volume), int(book), int(chapter), int(section), int(scroll))


def serialize(scrolls: dict) -> bytes:
    """dict {coord tuple -> text} -> phext bytes.

    Scrolls are emitted in ascending coordinate order. Between the
    current coordinate and the next target, delimiter bytes advance the
    coordinate one step at a time (each delimiter flushes one scroll
    slot; slots not present in the dict become empty scrolls, exactly as
    the reference parser reads them back).
    """
    items = sorted(scrolls.items())
    out = bytearray()
    cur = [1] * 9
    for coord, text in items:
        tgt = list(coord)
        if any(t < 1 for t in tgt):
            raise ValueError(f"coordinates are 1-indexed: {coord}")
        for ch in text:
            if ch in _DELIM_CHARS:
                raise ValueError(
                    f"scroll text at {coord_str(coord)} contains a phext "
                    f"delimiter character U+{ord(ch):04X}; phext scrolls "
                    "cannot contain raw delimiter bytes"
                )
        while cur != tgt:
            # Highest level that still differs; sorted order guarantees
            # cur[h] < tgt[h] there.
            h = next(i for i in range(9) if cur[i] != tgt[i])
            if cur[h] > tgt[h]:
                raise ValueError(
                    f"coordinates out of order: cannot go from "
                    f"{coord_str(cur)} to {coord_str(tgt)}"
                )
            out.append(LEVEL_BYTES[h])
            cur[h] += 1
            for j in range(h + 1, 9):
                cur[j] = 1
        out += text.encode("utf-8")
    return bytes(out)


def parse_bytes(data: bytes) -> dict:
    """phext bytes -> dict {coord tuple -> text}.

    Semantics identical to the reference dead-reckoning parser
    (parse_phext.py): walk the stream, a delimiter flushes the current
    scroll and advances its level (resetting lower levels to 1); a
    trailing empty scroll produced by a final delimiter is dropped,
    interior empty scrolls are kept. Delimiter bytes are single-byte
    ASCII controls and never occur inside a multi-byte UTF-8 sequence,
    so splitting on them is byte-safe.
    """
    scrolls = {}
    coord = [1] * 9
    pos = 0
    for m in _DELIM_RE.finditer(data):
        scrolls[tuple(coord)] = data[pos:m.start()].decode("utf-8", errors="replace")
        lvl = LEVEL_INDEX[m.group()[0]]
        coord[lvl] += 1
        for j in range(lvl + 1, 9):
            coord[j] = 1
        pos = m.end()
    tail = data[pos:].decode("utf-8", errors="replace")
    if tail != "" or not scrolls:
        scrolls[tuple(coord)] = tail
    return scrolls


def parse_file(path) -> dict:
    return parse_bytes(Path(path).read_bytes())


class PhextDoc:
    """An in-memory phext: dict of coord tuple -> scroll text."""

    def __init__(self, scrolls=None):
        self.scrolls = dict(scrolls) if scrolls else {}

    @classmethod
    def load(cls, path):
        return cls(parse_file(path))

    def save(self, path):
        Path(path).write_bytes(serialize(self.scrolls))

    # -- scroll access ------------------------------------------------------
    def set(self, coord, text: str):
        self.scrolls[tuple(coord)] = str(text)

    def get(self, coord, default=None):
        return self.scrolls.get(tuple(coord), default)

    def set_a(self, volume, book, chapter, section, scroll, text):
        self.set(addr(volume, book, chapter, section, scroll), text)

    def get_a(self, volume, book, chapter, section, scroll, default=None):
        return self.get(addr(volume, book, chapter, section, scroll), default)

    # -- volume bookkeeping ---------------------------------------------------
    def volumes(self):
        return sorted({c[4] for c in self.scrolls})

    def next_volume(self):
        vs = self.volumes()
        return (max(vs) + 1) if vs else 1

    def clear_volume(self, volume):
        self.scrolls = {c: t for c, t in self.scrolls.items() if c[4] != volume}

    def volume_kind(self, volume):
        return self.get_a(volume, BOOK_MANIFEST, 1, 1, 2)


# ---------------------------------------------------------------------------
# Manifest (written at 1.1.1/1.1.1/1.1.1 in every phext this module writes)
# ---------------------------------------------------------------------------
MANIFEST_TEXT = """\
Circulatory Fidelity x Phext — study phext manifest
====================================================
This phext stores Circulatory Fidelity (CF) diagnostics in phext form.
Phext: a text stream cut by ASCII control delimiters into scrolls at 9D
coordinates Library.Shelf.Series/Collection.Volume.Book/Chapter.Section.Scroll
(1-indexed; first scroll 1.1.1/1.1.1/1.1.1). Delimiters, highest to lowest:
0x01 Library, 0x1F Shelf, 0x1E Series, 0x1D Collection, 0x1C Volume,
0x1A Book, 0x19 Chapter, 0x18 Section, 0x17 Scroll. A delimiter increments
its level and resets all lower levels to 1. Coordinates are internal to
this file: the same coordinate in a different phext is a different scroll.

Why: CF decomposes a precision matrix by position, Lambda = D + R
(D diagonal i==j, retained by MFVI; R off-diagonal i!=j, discarded by
MFVI). Phext gives that split an address: D and R occupy distinct lanes
below, so one relational entry can be fetched by coordinate alone.

Coordinate map (Library.Shelf.Series/Collection fixed at 1.1.1/1):
  Volume  = one diagnostic run / matrix / windowed result / dataset.
  Book    = artifact-class lane:
    1 manifest/meta   2 summary+verdict   3 nodal D lane
    4 relational R lane   5 windowed series   6 dataset rows
    7 provenance

Book 1 (manifest/meta), Chapter 1, Section 1:
  Scroll 1 = for Volume 1: this manifest. For Volume >1: volume header.
  Scroll 2 = volume kind: diagnose | matrix | windowed | dataset
  Scroll 3 = volume label / source name
  Chapter 2, Section 1: Scroll 1.. = kind-specific counts
    (dataset: 1=n_data_rows, 2=n_columns; matrix: 1=dimension;
     diagnose: 1=n; windowed: 1=n_windows)

Book 2 (summary+verdict) layout is PER VOLUME KIND: read the kind
from Book 1 Chapter 1 Section 1 Scroll 2 first, then apply that
kind's layout below.
  diagnose volume, Chapter 1, Section 1, Scroll:
    1=ic 2=se 3=method 4=n 5=model_type 6=recommendation(verdict)
    7=rho (signed Pearson) 8=ic_band 9=mse_ratio_predicted
  diagnose volume, Chapter 2/3/4, Section 1, Scroll 1..4 =
    correlation / precision / covariance 2x2 entries, row-major.
  matrix volume, Chapter 1, Section 1, Scroll:
    1=kind 2=dimension 3=max_abs_offdiagonal 4=coordinate of that entry's
    value scroll 5=frobenius_norm_R 6=frobenius_norm_D
  windowed volume, Chapter 1, Section 1, Scroll:
    1=ic_max 2=ic_mean 3=ic_std 4=n_windows 5=window_size
    6=recommendation(verdict) 7=max_window_index (0-based)
    8=coordinate of the max window's IC scroll (Book 5)
    9=se_per_window (may be absent)

Matrix entry records (Books 3 and 4), per-field Scroll:
  1=value 2=ic_marginal 3=se 4=i (0-based) 5=j (0-based) 6=kind
  7=ic_conditional 8=partial_corr (signed) 9=marginal_corr (signed)
  (Fields 7-9 exist only in schema 2 files; see 'Schema versions'.)
  D lane (Book 3): Chapter=i+1, Section=1, diagonal entries only (i==j).
  R lane (Book 4): Chapter=i+1, Section=j+1, off-diagonal only (i!=j).
  MARGINAL vs CONDITIONAL coupling (schema 2): two different relational
  properties, stored separately, never mixed in one field.
  Field 2 'ic_marginal' and field 9 'marginal_corr' describe MARGINAL
  dependence: field 9 is the signed correlation corr_ij of the pair
  (for kind=correlation, the entry itself; for kind=covariance, the
  normalized covariance; for kind=precision, corr_ij of the implied
  covariance inv(Lambda)); field 2 is its magnitude. Marginal
  correlation can be entirely indirect: in the chain A-B-C with
  corr(A,B)=corr(B,C)=0.6, corr(A,C)=0.36 although A and C are
  conditionally independent given B. Upstream CF's pairwise IC is a
  marginal (copula) quantity by design; field 2 has stored exactly
  that since schema 1, and keeps that meaning.
  Field 7 'ic_conditional' and field 8 'partial_corr' describe
  CONDITIONAL dependence: the partial correlation of the pair given
  all other variables, partial_ij = -Lambda_ij/sqrt(Lambda_ii *
  Lambda_jj) from the precision matrix Lambda (for kind=precision,
  the stored matrix itself; otherwise inv of the stored matrix).
  In a Gaussian graphical model the precision matrix encodes
  conditional independence: Lambda_ij = 0 iff the pair is
  conditionally independent given the rest. For the chain above,
  partial(A,C) = 0 while marginal(A,C) = 0.36. On the diagonal both
  couplings are 1.0 by convention. Neither marginal nor conditional
  dependence is causation.
  'se' field: Fisher SE 1/sqrt(n-3) when a sample size n is known,
  otherwise the text 'nan'.
  Floats are stored with repr(), which round-trips float64 exactly.

Schema versions: Book 7 Scroll 6 (below) records the schema version
of each volume ('2' for files written by cf_phext >= 1.1.0). Files
without fields 7-9 and without Book 7 Scroll 6 are schema 1; their
conditional coupling is not stored but can be derived at read time
from the stored value fields (see coupling_from_phext in cf_phext.py;
the signed marginal correlation is NOT stored in schema 1 and cannot
be recovered from the entry record alone).

Book 5 (windowed series): Chapter=window_index+1, Section=1, Scroll:
  1=ic 2=window_center 3=is_max ('true'/'false')

Book 6 (dataset rows): Chapter 1 = header row (Scroll=column#+1, the
  header field text verbatim). Data row r (0-based) at Chapter=r+2,
  Section=1, Scroll=column#+1, the ORIGINAL field text verbatim
  (post-CSV-parse field value). Export re-serializes with the csv
  module (QUOTE_MINIMAL, LF line endings), reproducing the source
  files byte-for-byte for the upstream simulations/*.csv files.

Book 7 (provenance), Chapter 1, Section 1, Scroll:
  1=upstream CF version (or 'n/a') 2=UTC timestamp (ISO 8601)
  3=source description 4=cf_phext version 5=detail (model_type/method,
  matrix kind, or SHA-256 of an imported CSV)
  6=schema version ('2'; ABSENT in schema 1 files — absence means '1')

Thresholds (filtering): IC <0.25 safe, 0.25-0.35 weak, 0.35-0.55
moderate, 0.55-0.70 strong, >0.70 structured inference required.
Pooling inverts the reading (see upstream diagnose()).
"""


def _ensure_manifest(doc: PhextDoc):
    key = addr(1, BOOK_MANIFEST, 1, 1, 1)
    if doc.get(key) in (None, ""):
        doc.set(key, MANIFEST_TEXT)


def _stamp_meta(doc, volume, kind, label, counts, provenance_detail, cf_version, source):
    if volume != 1:
        doc.set_a(volume, BOOK_MANIFEST, 1, 1, 1,
                  f"Volume {volume}: kind={kind}, label={label}")
    doc.set_a(volume, BOOK_MANIFEST, 1, 1, 2, kind)
    doc.set_a(volume, BOOK_MANIFEST, 1, 1, 3, label)
    for k, v in enumerate(counts, start=1):
        doc.set_a(volume, BOOK_MANIFEST, 2, 1, k, str(v))
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 1, cf_version)
    # Timestamp: current UTC, unless CF_PHEXT_TIMESTAMP is set (any
    # string), which makes generated phexts byte-reproducible.
    ts = os.environ.get("CF_PHEXT_TIMESTAMP") or \
        datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 2, ts)
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 3, source)
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 4, __version__)
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 5, provenance_detail)
    doc.set_a(volume, BOOK_PROVENANCE, 1, 1, 6, SCHEMA_VERSION)


def _open_doc(path, volume):
    """Load existing phext (append) or start a new one; resolve volume."""
    p = Path(path)
    if p.exists() and p.stat().st_size > 0:
        doc = PhextDoc.load(p)
        if volume is None:
            volume = doc.next_volume()
        else:
            doc.clear_volume(volume)
    else:
        doc = PhextDoc()
        if volume is None:
            volume = 1
    _ensure_manifest(doc)
    if volume == 1:
        # Volume 1's meta scroll 1 is the manifest by definition.
        doc.set(addr(1, BOOK_MANIFEST, 1, 1, 1), MANIFEST_TEXT)
    return doc, volume


def _f(v) -> str:
    """Float -> repr text (round-trips float64 exactly)."""
    return repr(float(v))


def ic_band(ic: float, model_type: str = "filtering") -> str:
    """Threshold band for an IC value. The filtering scale is upstream's
    interpretive scale; pooling inverts the reading (upstream diagnose()
    recommendation text is the authoritative verdict in both cases)."""
    if model_type == "pooling":
        if ic > 0.70:
            return "pooling: very strong separation"
        if ic > 0.55:
            return "pooling: strong separation"
        if ic > 0.35:
            return "pooling: moderate homogeneity"
        return "pooling: weak separation"
    if ic < 0.25:
        return "safe"
    if ic < 0.35:
        return "weak"
    if ic < 0.55:
        return "moderate"
    if ic < 0.70:
        return "strong"
    return "structured required"


# ---------------------------------------------------------------------------
# Upstream CF module loading (by file path; no installation required)
# ---------------------------------------------------------------------------
def load_cf(cf_module_path=None):
    """Import upstream circulatory_fidelity.py by path.

    Resolution order: explicit argument, $CF_PHEXT_UPSTREAM, a sibling
    circulatory_fidelity.py (repo layout: python/cf_phext.py next to
    python/circulatory_fidelity.py), then the phext-integration layout
    (../upstream/python/circulatory_fidelity.py relative to this file).
    """
    candidates = []
    if cf_module_path:
        candidates.append(Path(cf_module_path))
    if os.environ.get("CF_PHEXT_UPSTREAM"):
        candidates.append(Path(os.environ["CF_PHEXT_UPSTREAM"]))
    here = Path(__file__).resolve()
    candidates.append(here.parent / "circulatory_fidelity.py")
    candidates.append(here.parent.parent / "upstream" / "python" / "circulatory_fidelity.py")
    candidates.append(here.parent.parent.parent / "upstream" / "python" / "circulatory_fidelity.py")
    for cand in candidates:
        if cand.is_file():
            spec = importlib.util.spec_from_file_location("circulatory_fidelity", cand)
            mod = importlib.util.module_from_spec(spec)
            # Register before exec: upstream uses @dataclass, which looks
            # the module up in sys.modules during class creation.
            sys.modules["circulatory_fidelity"] = mod
            spec.loader.exec_module(mod)
            return mod
    raise FileNotFoundError(
        "circulatory_fidelity.py not found; pass cf_module_path or set "
        f"CF_PHEXT_UPSTREAM. Tried: {[str(c) for c in candidates]}"
    )


# ---------------------------------------------------------------------------
# Matrix decomposition -> phext (Books 3 + 4)
# ---------------------------------------------------------------------------
def _associated_corr(matrix: np.ndarray, kind: str) -> np.ndarray:
    """Marginal correlation matrix associated with `matrix`.

    This is the correlation whose magnitude field 2 (ic_marginal) has
    always stored: for kind=correlation the matrix itself, for
    kind=covariance the normalized covariance, for kind=precision the
    correlation of the implied covariance inv(Lambda). Marginal
    correlation includes indirect dependence (see partial_correlations
    for the conditional counterpart)."""
    M = np.asarray(matrix, dtype=np.float64)
    if kind == "correlation":
        return M.copy()
    if kind == "covariance":
        d = np.sqrt(np.diag(M))
        return M / np.outer(d, d)
    if kind == "precision":
        cov = np.linalg.inv(M)
        d = np.sqrt(np.diag(cov))
        return cov / np.outer(d, d)
    raise ValueError(f"unknown matrix kind: {kind!r}")


def _partial_from_precision(P: np.ndarray) -> np.ndarray:
    """Signed partial correlations from a precision matrix P:
    partial_ij = -P_ij / sqrt(P_ii * P_jj) for i != j; 1.0 on the
    diagonal by convention. Pure normalization — NO matrix inverse."""
    d = np.sqrt(np.diag(P))
    partial = -P / np.outer(d, d)
    np.fill_diagonal(partial, 1.0)
    return partial


def partial_correlations(matrix: np.ndarray, kind: str = "precision") -> np.ndarray:
    """Signed partial (conditional) correlation matrix for `matrix`.

    partial_ij is the correlation of variables i and j conditioned on
    all remaining variables. In a Gaussian graphical model it is read
    off the precision matrix Lambda: partial_ij =
    -Lambda_ij / sqrt(Lambda_ii * Lambda_jj).

    Cost note: for kind='precision' the input IS Lambda, so this is a
    normalization only and performs NO matrix inversion. For
    kind='correlation'/'covariance', Lambda = inv(matrix) is required
    and exactly one inverse is computed. The dual holds for marginal
    correlations: direct for correlation/covariance inputs, one
    inverse for precision inputs."""
    M = np.asarray(matrix, dtype=np.float64)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("matrix must be square 2D")
    if kind == "precision":
        P = M
    elif kind in ("correlation", "covariance"):
        P = np.linalg.inv(M)
    else:
        raise ValueError(f"unknown matrix kind: {kind!r}")
    return _partial_from_precision(P)


def conditional_ic(matrix: np.ndarray, kind: str = "precision") -> np.ndarray:
    """|partial correlation| matrix — the conditional counterpart of
    CF's (marginal) Inference Coupling. See partial_correlations for
    the definition and cost notes; for kind='precision' no inverse is
    performed."""
    return np.abs(partial_correlations(matrix, kind))


# Backwards-compatible alias matching the module's naming style.
marginal_correlations = _associated_corr


def _corr_and_precision(M: np.ndarray, kind: str):
    """(marginal corr, precision) for a stored matrix, computing the
    dense inverse EXACTLY ONCE regardless of kind:
      precision:   corr needs inv(M); the precision is M itself.
      covariance:  corr is a normalization of M; precision is inv(M).
      correlation: corr is M itself;      precision is inv(M)."""
    if kind == "precision":
        cov = np.linalg.inv(M)
        d = np.sqrt(np.diag(cov))
        return cov / np.outer(d, d), M
    if kind == "covariance":
        d = np.sqrt(np.diag(M))
        return M / np.outer(d, d), np.linalg.inv(M)
    if kind == "correlation":
        return M.copy(), np.linalg.inv(M)
    raise ValueError(f"unknown matrix kind: {kind!r}")


def _write_matrix_lanes(doc, volume, matrix, kind, n=None):
    M = np.asarray(matrix, dtype=np.float64)
    if M.ndim != 2 or M.shape[0] != M.shape[1]:
        raise ValueError("matrix must be square 2D")
    dim = M.shape[0]
    # One inverse total (see _corr_and_precision): schema 1 computed the
    # marginal correlations with the same single inverse for
    # kind=precision and none for the other kinds; fields 1-6 below are
    # byte-identical to schema 1 for identical inputs.
    corr, prec = _corr_and_precision(M, kind)
    partial = _partial_from_precision(prec)
    se_text = _f(1.0 / np.sqrt(n - 3)) if (n is not None and n > 3) else "nan"
    D = np.diag(np.diag(M))
    R = M - D
    for i in range(dim):
        for j in range(dim):
            book = BOOK_D if i == j else BOOK_R
            section = 1 if i == j else j + 1
            base = (volume, book, i + 1, section)
            doc.set_a(*base, 1, _f(M[i, j]))
            doc.set_a(*base, 2, _f(abs(corr[i, j])))
            doc.set_a(*base, 3, se_text)
            doc.set_a(*base, 4, str(i))
            doc.set_a(*base, 5, str(j))
            doc.set_a(*base, 6, kind)
            # Schema 2 (additive): conditional coupling, explicitly.
            doc.set_a(*base, 7, _f(abs(partial[i, j])))
            doc.set_a(*base, 8, _f(partial[i, j]))
            doc.set_a(*base, 9, _f(corr[i, j]))
    # Summary scrolls for the matrix volume.
    offdiag_abs = np.abs(R)
    if dim > 1:
        k = int(np.argmax(offdiag_abs))
        mi, mj = divmod(k, dim)
        max_coord = coord_str(addr(volume, BOOK_R, mi + 1, mj + 1, 1))
        max_val = float(offdiag_abs[mi, mj])
    else:
        max_coord, max_val = "", 0.0
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 1, kind)
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 2, str(dim))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 3, _f(max_val))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 4, max_coord)
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 5, _f(np.linalg.norm(R)))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 6, _f(np.linalg.norm(D)))
    return D, R


def matrix_to_phext(matrix, path, kind="precision", n=None, volume=None,
                    label=None):
    """Write a square matrix to `path` as a new (or specified) volume,
    decomposed into the nodal D lane (Book 3) and relational R lane
    (Book 4). Returns the volume number used."""
    M = np.asarray(matrix, dtype=np.float64)
    doc, volume = _open_doc(path, volume)
    _write_matrix_lanes(doc, volume, M, kind, n=n)
    _stamp_meta(doc, volume, "matrix", label or f"{kind} matrix {M.shape[0]}x{M.shape[1]}",
                [M.shape[0]], f"kind={kind}",
                "n/a (matrix input)", "matrix_to_phext")
    doc.save(path)
    return volume


def decomposition_from_phext(path, volume=None):
    """Read back (D, R) matrices for a matrix/diagnose volume."""
    doc = PhextDoc.load(path)
    if volume is None:
        kinds = [(v, doc.volume_kind(v)) for v in doc.volumes()]
        mats = [v for v, k in kinds if k in ("matrix", "diagnose")]
        if not mats:
            raise ValueError("no matrix/diagnose volume in this phext")
        volume = mats[-1]
    dim_text = doc.get_a(volume, BOOK_SUMMARY, 1, 1, 2)
    if doc.volume_kind(volume) == "diagnose":
        dim = 2
    elif dim_text is not None:
        dim = int(dim_text)
    else:
        raise ValueError(f"volume {volume} has no matrix dimension record")
    D = np.zeros((dim, dim))
    R = np.zeros((dim, dim))
    for i in range(dim):
        t = doc.get_a(volume, BOOK_D, i + 1, 1, 1)
        if t is not None:
            D[i, i] = float(t)
        for j in range(dim):
            if i == j:
                continue
            t = doc.get_a(volume, BOOK_R, i + 1, j + 1, 1)
            if t is not None:
                R[i, j] = float(t)
    return D, R


def matrix_from_phext(path, volume=None) -> np.ndarray:
    """Reconstruct the full matrix D + R from a phext volume."""
    D, R = decomposition_from_phext(path, volume)
    return D + R


def schema_version_of(doc: PhextDoc, volume) -> str:
    """Schema version of a volume: Book 7 Scroll 6, or '1' if absent
    (schema 1 files predate the stamp)."""
    return doc.get_a(volume, BOOK_PROVENANCE, 1, 1, 6) or "1"


def coupling_from_phext(path, volume, i, j) -> dict:
    """Read one matrix entry's coupling record from a phext volume.

    Returns a dict with keys:
      value           the stored matrix entry (field 1)
      ic_marginal     |marginal corr_ij| (field 2; CF's pairwise IC
                      for this pair — a marginal quantity by design)
      ic_conditional  |partial corr_ij| (field 7, schema 2; derived
                      from the stored values for schema 1 files)
      partial_corr    signed partial correlation (field 8 / derived)
      marginal_corr   signed marginal correlation (field 9, schema 2;
                      None for schema 1 — the sign was not stored)
      se              Fisher SE stored in field 3 (float; nan if the
                      writer knew no sample size)
      kind            'precision' | 'correlation' | 'covariance'
      schema_version  '2' if the volume carries fields 7-9 / the
                      Book 7 Scroll 6 stamp, else '1'
      notes           list of derivation/limitation notes (schema 1)

    Schema-1 derivation: the reconstructed matrix Mat = D + R is the
    stored matrix. If kind is 'precision', Mat IS the precision matrix
    and partial_ij = -Mat_ij/sqrt(Mat_ii*Mat_jj) directly. If kind is
    'correlation'/'covariance', the precision is inv(Mat) (one
    read-time inverse). The signed marginal correlation cannot be
    derived from a schema-1 entry record (only its magnitude was
    stored), so marginal_corr is None there rather than a guess."""
    doc = PhextDoc.load(path)
    if volume is None:
        raise ValueError("coupling_from_phext requires an explicit volume")
    i, j = int(i), int(j)
    book = BOOK_D if i == j else BOOK_R
    section = 1 if i == j else j + 1
    base = (volume, book, i + 1, section)

    def _get(scroll):
        return doc.get_a(*base, scroll)

    kind = _get(6)
    if kind is None:
        raise ValueError(
            f"no matrix entry record at volume {volume}, i={i}, j={j}")
    value = float(_get(1))
    ic_marginal = float(_get(2))
    se_text = _get(3)
    se = float(se_text) if se_text is not None else float("nan")
    schema = schema_version_of(doc, volume)
    notes = []
    if _get(7) is not None:
        ic_conditional = float(_get(7))
        partial_corr = float(_get(8))
        marginal_corr = float(_get(9))
    else:
        # Schema 1: derive conditional coupling from stored values.
        Mat = matrix_from_phext(path, volume)
        if kind == "precision":
            P = Mat
        else:
            P = np.linalg.inv(Mat)
            notes.append(
                "ic_conditional derived at read time via one inverse of "
                "the reconstructed stored matrix (schema 1 stores no "
                "conditional fields)")
        if i == j:
            partial_corr = 1.0
        else:
            partial_corr = float(-P[i, j] / np.sqrt(P[i, i] * P[j, j]))
        ic_conditional = abs(partial_corr)
        marginal_corr = None
        notes.append(
            "marginal_corr (signed) is not stored in schema 1; only "
            "its magnitude (ic_marginal) was recorded, so the sign is "
            "reported as None rather than reconstructed")
        if schema == "1":
            notes.append("schema 1 file: fields 7-9 absent; conditional "
                         "values derived from the stored value fields")
    return {
        "value": value,
        "ic_marginal": ic_marginal,
        "ic_conditional": ic_conditional,
        "partial_corr": partial_corr,
        "marginal_corr": marginal_corr,
        "marginal_corr_signed": marginal_corr,
        "se": se,
        "kind": kind,
        "schema_version": schema,
        "notes": notes,
    }


# ---------------------------------------------------------------------------
# diagnose() -> phext
# ---------------------------------------------------------------------------
def diagnose_to_phext(x, y, path, model_type="filtering", method="copula",
                      cf_module_path=None, volume=None, label=None):
    """Run upstream diagnose(x, y) and write the result as a phext volume:
    summary+verdict (Book 2), the 2x2 precision decomposition in the
    D/R lanes (Books 3/4), correlation/precision/covariance entries in
    Book 2 Chapters 2-4, manifest + provenance. Returns a dict with the
    ICDiagnostic, the volume used, and key coordinates."""
    cf = load_cf(cf_module_path)
    x = np.asarray(x, dtype=np.float64).flatten()
    y = np.asarray(y, dtype=np.float64).flatten()
    diag = cf.diagnose(x, y, model_type=model_type, method=method)

    rho = float(np.corrcoef(x, y)[0, 1])
    corr = np.array([[1.0, rho], [rho, 1.0]])
    cov = np.cov(x, y)
    prec = np.linalg.inv(corr) if abs(rho) < 1.0 else np.linalg.pinv(corr)

    doc, volume = _open_doc(path, volume)
    # D/R lanes carry the precision decomposition Lambda = D + R.
    _write_matrix_lanes(doc, volume, prec, "precision", n=diag.n)
    # Summary + verdict.
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 1, _f(diag.ic))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 2, _f(diag.se))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 3, diag.method)
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 4, str(diag.n))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 5, model_type)
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 6, diag.recommendation)
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 7, _f(rho))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 8, ic_band(diag.ic, model_type))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 9, _f(cf.mse_ratio_predicted(diag.ic)))
    for ch, M in ((2, corr), (3, prec), (4, cov)):
        for k in range(4):
            i, j = divmod(k, 2)
            doc.set_a(volume, BOOK_SUMMARY, ch, 1, k + 1, _f(M[i, j]))
    _stamp_meta(doc, volume, "diagnose", label or f"diagnose ({model_type})",
                [diag.n], f"model_type={model_type}, method={method}",
                getattr(cf, "__version__", "unknown"),
                "diagnose_to_phext via upstream circulatory_fidelity.diagnose")
    doc.save(path)
    return {
        "diagnostic": diag,
        "volume": volume,
        "ic_coord": coord_str(addr(volume, BOOK_SUMMARY, 1, 1, 1)),
        "verdict_coord": coord_str(addr(volume, BOOK_SUMMARY, 1, 1, 6)),
        "rho": rho,
        "correlation": corr,
        "precision": prec,
        "covariance": cov,
    }


# ---------------------------------------------------------------------------
# windowed_ic() result -> phext
# ---------------------------------------------------------------------------
def windowed_to_phext(result, path, volume=None, label=None):
    """Write a windowed_ic() result dict as a phext volume. Each window
    is addressed in Book 5 (Chapter = window#+1); the Book 2 verdict
    scrolls record IC_max and the coordinate of the window scroll that
    holds it — the Maximal Coupling Rule, made addressable.
    Returns (volume, max_window_coord_str)."""
    ic_series = [float(v) for v in result["ic_series"]]
    centers = [int(v) for v in result["window_centers"]]
    if not ic_series:
        raise ValueError("windowed result has no windows")
    max_idx = int(np.argmax(ic_series))
    doc, volume = _open_doc(path, volume)
    for w, (ic, center) in enumerate(zip(ic_series, centers)):
        doc.set_a(volume, BOOK_WINDOWED, w + 1, 1, 1, _f(ic))
        doc.set_a(volume, BOOK_WINDOWED, w + 1, 1, 2, str(center))
        doc.set_a(volume, BOOK_WINDOWED, w + 1, 1, 3,
                  "true" if w == max_idx else "false")
    max_coord = coord_str(addr(volume, BOOK_WINDOWED, max_idx + 1, 1, 1))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 1, _f(result["ic_max"]))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 2, _f(result["ic_mean"]))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 3, _f(result["ic_std"]))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 4, str(int(result["n_windows"])))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 5, str(int(result["window_size"])))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 6, str(result["recommendation"]))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 7, str(max_idx))
    doc.set_a(volume, BOOK_SUMMARY, 1, 1, 8, max_coord)
    if result.get("se_per_window") is not None:
        doc.set_a(volume, BOOK_SUMMARY, 1, 1, 9, _f(result["se_per_window"]))
    _stamp_meta(doc, volume, "windowed", label or "windowed_ic result",
                [int(result["n_windows"])],
                f"window_size={int(result['window_size'])}",
                "n/a (windowed result dict)", "windowed_to_phext")
    doc.save(path)
    return volume, max_coord


# ---------------------------------------------------------------------------
# CSV <-> phext (field-text verbatim, round-trip exact)
# ---------------------------------------------------------------------------
def import_csv_to_phext(csv_path, phext_path, volume=None):
    """Import a CSV as a new phext volume (Book 6 dataset lane).
    Header row at Chapter 1; data row r at Chapter r+2; each field is
    one scroll holding the parsed field text verbatim. Returns volume."""
    csv_path = Path(csv_path)
    raw = csv_path.read_bytes()
    rows = list(csv.reader(io.StringIO(raw.decode("utf-8"), newline="")))
    if not rows:
        raise ValueError(f"empty CSV: {csv_path}")
    header, data = rows[0], rows[1:]
    doc, volume = _open_doc(phext_path, volume)
    for c, field in enumerate(header):
        doc.set_a(volume, BOOK_DATASET, 1, 1, c + 1, field)
    for r, row in enumerate(data):
        for c, field in enumerate(row):
            doc.set_a(volume, BOOK_DATASET, r + 2, 1, c + 1, field)
    _stamp_meta(doc, volume, "dataset", csv_path.name,
                [len(data), len(header)],
                f"sha256={hashlib.sha256(raw).hexdigest()}",
                "n/a (dataset import)", f"import_csv_to_phext from {csv_path}")
    doc.save(phext_path)
    return volume


def export_phext_to_csv(phext_path, volume, out_path=None):
    """Export a dataset volume back to CSV text. Re-serializes with the
    csv module (QUOTE_MINIMAL, LF endings). Writes to out_path if given.
    Returns the CSV text."""
    doc = PhextDoc.load(phext_path)
    if doc.volume_kind(volume) != "dataset":
        raise ValueError(f"volume {volume} is not a dataset volume")
    counts = [doc.get_a(volume, BOOK_MANIFEST, 2, 1, k) for k in (1, 2)]
    n_rows, n_cols = int(counts[0]), int(counts[1])
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, lineterminator="\n")
    header = [doc.get_a(volume, BOOK_DATASET, 1, 1, c + 1, "") for c in range(n_cols)]
    writer.writerow(header)
    for r in range(n_rows):
        writer.writerow(
            [doc.get_a(volume, BOOK_DATASET, r + 2, 1, c + 1, "") for c in range(n_cols)]
        )
    text = buf.getvalue()
    if out_path:
        Path(out_path).write_text(text, encoding="utf-8", newline="")
    return text


# ---------------------------------------------------------------------------
# Single-coordinate query / TOC
# ---------------------------------------------------------------------------
def read_scroll(path, coord_str_) -> str:
    """Fetch one scroll by coordinate string, e.g. '1.1.1/1.1.4/2.1.1'.
    Raises KeyError if the coordinate holds no scroll."""
    doc = PhextDoc.load(path)
    coord = parse_coord(coord_str_)
    if coord not in doc.scrolls:
        raise KeyError(f"no scroll at {coord_str_}")
    return doc.scrolls[coord]


def toc(path):
    """Sorted list of (coord_str, text) for non-empty scrolls."""
    doc = PhextDoc.load(path)
    return [(coord_str(c), t) for c, t in sorted(doc.scrolls.items()) if t != ""]


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def _cmd_diagnose_demo(args):
    rng = np.random.default_rng(20261101)
    n = 2000
    x = rng.standard_normal(n)
    y = 0.8 * x + 0.6 * rng.standard_normal(n)
    out = args.out or "cf_demo.phext"
    if Path(out).exists():
        Path(out).unlink()
    info = diagnose_to_phext(x, y, out, model_type="filtering")
    d = info["diagnostic"]
    print(f"IC = {d.ic:.4f} +/- {d.se:.4f}  (n={d.n}, method={d.method})")
    print(f"verdict @ {info['verdict_coord']}: {d.recommendation}")
    print(f"phext written: {out} (volume {info['volume']})")
    return 0


def _cmd_matrix_demo(args):
    rng = np.random.default_rng(7)
    A = rng.standard_normal((4, 4))
    M = A @ A.T + 4.0 * np.eye(4)
    out = args.out or "cf_matrix_demo.phext"
    if Path(out).exists():
        Path(out).unlink()
    vol = matrix_to_phext(M, out, kind="precision")
    back = matrix_from_phext(out, vol)
    err = float(np.max(np.abs(back - M)))
    print(f"4x4 precision matrix -> phext volume {vol} -> reconstructed")
    print(f"max |reconstruction error| = {err:.3e}")
    print(f"phext written: {out}")
    return 0 if err < 1e-12 else 1


def _cmd_import(args):
    vol = import_csv_to_phext(args.csv, args.phext)
    print(f"imported {args.csv} as volume {vol} of {args.phext}")
    return 0


def _cmd_export(args):
    text = export_phext_to_csv(args.phext, args.volume, args.out)
    if not args.out:
        sys.stdout.write(text)
    else:
        print(f"exported volume {args.volume} of {args.phext} -> {args.out}")
    return 0


def _cmd_query(args):
    try:
        print(read_scroll(args.phext, args.coord))
    except KeyError as e:
        print(str(e), file=sys.stderr)
        return 1
    return 0


def _cmd_toc(args):
    for c, t in toc(args.phext):
        first = t.splitlines()[0] if t.splitlines() else ""
        print(f"{c}: {first[:100]}")
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(prog="cf_phext", description=__doc__.splitlines()[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("diagnose-demo"); p.add_argument("--out"); p.set_defaults(fn=_cmd_diagnose_demo)
    p = sub.add_parser("matrix-demo"); p.add_argument("--out"); p.set_defaults(fn=_cmd_matrix_demo)
    p = sub.add_parser("import"); p.add_argument("csv"); p.add_argument("phext"); p.set_defaults(fn=_cmd_import)
    p = sub.add_parser("export"); p.add_argument("phext"); p.add_argument("volume", type=int)
    p.add_argument("--out"); p.set_defaults(fn=_cmd_export)
    p = sub.add_parser("query"); p.add_argument("phext"); p.add_argument("coord"); p.set_defaults(fn=_cmd_query)
    p = sub.add_parser("toc"); p.add_argument("phext"); p.set_defaults(fn=_cmd_toc)
    args = ap.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    raise SystemExit(main())
