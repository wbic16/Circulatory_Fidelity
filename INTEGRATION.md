# INTEGRATION.md — Circulatory Fidelity × Phext: design

## 1. First-principles derivation

Circulatory Fidelity's core object is a split **by matrix position**:

```
Λ = D + R     D[i,j] = Λ[i,j] if i == j else 0      (nodal, kept by MFVI)
              R[i,j] = Λ[i,j] if i != j else 0      (relational, dropped by MFVI)
```

Everything CF predicts follows from how much structure sits in R
relative to D. But in every existing representation — a NumPy array, a
CSV, a Julia `Matrix` — position is *implicit*: to read one relational
entry you load and index the whole matrix, and D-vs-R is a convention
in the reader's head, not a property of the stored artifact.

Phext is a text encoding in which **position is explicit**: a stream
cut by ASCII control delimiters into scrolls, each living at a
9-dimensional coordinate
`Library.Shelf.Series/Collection.Volume.Book/Chapter.Section.Scroll`
(1-indexed; delimiters, highest→lowest: `0x01, 0x1F, 0x1E, 0x1D, 0x1C,
0x1A, 0x19, 0x18, 0x17`; hitting a delimiter increments that level and
resets all lower levels to 1). Coordinates are internal to a phext
file — the same coordinate in another phext is a different scroll.

Put the two together and the design is forced, not decorative:

- CF splits structure by position ⇒ give each position-class an
  **address range**. D and R become distinct *lanes* (Books 3 and 4),
  so "the part MFVI discards" is a place, not a computation.
- A matrix entry is a small record (value, its IC, its SE, its i, j)
  ⇒ store it as **per-field scrolls** at the coordinate derived from
  (i, j). One relational entry is then one coordinate fetch — by the
  module's `read_scroll`, or by SQ's `select` — with no matrix parse.
- CF's workflow objects (a `diagnose()` verdict, a `windowed_ic()`
  result, a validation table) are also records ⇒ each gets a Volume,
  and its summary/verdict gets fixed, documented scrolls.
- The Maximal Coupling Rule says a windowed diagnosis turns on the
  single worst window ⇒ the verdict scroll stores **the coordinate of
  that window's IC scroll**. The rule becomes an address.
- Validation CSVs are flat rows whose relational identity was
  collapsed into columns ⇒ store each field verbatim as its own
  scroll (row → Chapter, column → Scroll), so the table round-trips
  byte-exact and any field is addressable (this is also how
  Julia-produced tables enter the same phext: via CSV).

Two phext properties matter and are stated from the format's
definition, not by analogy: a scroll's content is text and the
delimiters are ASCII control bytes (themselves valid UTF-8), so a
conforming phext is valid UTF-8 as a whole byte stream; and arbitrary
binary, if ever needed, would ride inside scrolls as text (e.g.
base64). This integration stores only text.

## 2. Coordinate map

`Library.Shelf.Series/Collection` is fixed at `1.1.1/1` for study
phexts. `Volume` = one diagnostic run / matrix / windowed result /
imported dataset, incrementing per volume. `Book` = artifact lane:

| Book | Lane | Contents |
|-----:|------|----------|
| 1 | manifest/meta | global manifest at Vol 1 (see below); per-volume kind, label, counts |
| 2 | summary+verdict | diagnose/matrix/windowed summary fields (below) |
| 3 | nodal D | diagonal entries only |
| 4 | relational R | off-diagonal entries only |
| 5 | windowed series | one Chapter per window |
| 6 | dataset rows | one Chapter per row, one Scroll per column |
| 7 | provenance | upstream version, UTC timestamp, source, detail |

Every phext written by `cf_phext.py` carries the full map as its
manifest scroll at `1.1.1/1.1.1/1.1.1`, so a file is self-describing
without this document.

**Book 1** — Ch.1 §1: S1 = manifest (Vol 1) / volume header (Vol >1);
S2 = kind (`diagnose|matrix|windowed|dataset`); S3 = label/source
name. Ch.2 §1: kind-specific counts (dataset: rows, columns; matrix:
dimension; diagnose: n; windowed: n_windows).

**Book 2, diagnose volume** — Ch.1 §1 scrolls: 1=ic, 2=se, 3=method,
4=n, 5=model_type, 6=recommendation (the verdict), 7=rho (signed),
8=ic_band, 9=mse_ratio_predicted. Ch.2/3/4 §1 S1–4: the 2×2
correlation / precision / covariance entries, row-major.

**Book 2, matrix volume** — Ch.1 §1: 1=kind, 2=dimension,
3=max |off-diagonal|, 4=coordinate of that entry's value scroll,
5=‖R‖_F, 6=‖D‖_F.

**Book 2, windowed volume** — Ch.1 §1: 1=ic_max, 2=ic_mean, 3=ic_std,
4=n_windows, 5=window_size, 6=recommendation, 7=max window index
(0-based), 8=coordinate of the max window's IC scroll, 9=se_per_window
(when upstream provides it).

**Matrix entry records (Books 3/4)** — Chapter = i+1; Section = 1 in
the D lane, j+1 in the R lane; per-field Scroll: 1=value,
2=ic_marginal, 3=se, 4=i (0-based), 5=j (0-based), 6=kind, and —
**schema 2 (additive)** — 7=ic_conditional, 8=partial_corr (signed),
9=marginal_corr (signed). Fields 1–6 are byte-identical to schema 1
for identical inputs; field 2 is now *documented* as what it always
stored, the marginal coupling (§8). `se` is the Fisher SE 1/√(n−3)
when a sample size is known, else the text `nan`. Floats are stored
via `repr()`, which round-trips float64 exactly — that is what makes
the 1e-12 reconstruction tolerances in the tests honest (the seeded
4×4 demo reconstructs with error exactly 0.0).

**Book 5** — Chapter = window#+1, §1: 1=ic, 2=window_center,
3=is_max (`true`/`false`).

**Book 6** — Chapter 1 = header row (Scroll = column#+1, header text
verbatim). Data row r (0-based) at Chapter r+2, §1, Scroll =
column#+1, holding the parsed field text verbatim. Export
re-serializes with Python's `csv` module (QUOTE_MINIMAL, LF endings),
which reproduces all 8 upstream `simulations/*.csv` files
byte-for-byte — including the quoted, comma-containing `ic_bin`
fields of `threshold_calibration.csv`, which are stored as single
scrolls and re-quoted on export (tested explicitly).

**Book 7** — Ch.1 §1: 1=upstream CF version (or `n/a`), 2=UTC
timestamp, 3=source description, 4=cf_phext version, 5=detail
(model_type/method, matrix kind, or the SHA-256 of an imported CSV),
6=**schema version** (`2`; absent in schema 1 files, whose readers
treat absence as `1`).

For a `diagnose` volume, the D/R lanes hold the **precision**
decomposition of the pair (Λ from the inverse correlation matrix) —
the decomposition CF is about — while the correlation matrix itself
is in Book 2 Ch.2 for convenience.

## 3. How each upstream API maps

| Upstream (`circulatory_fidelity.py` v1.1.1) | Phext bridge (`cf_phext.py`) |
|---|---|
| `diagnose(x, y, model_type, method)` → `ICDiagnostic` | `diagnose_to_phext(...)` → Book 2 verdict + precision D/R lanes + provenance |
| `inference_coupling`, `ic_gaussian` | values land in Book 2 S1/S7 and entry `ic` fields; thresholds via `ic_band()` mirroring upstream's scale |
| precision/correlation/covariance matrix (user-supplied) | `matrix_to_phext(...)` / `matrix_from_phext(...)` / `decomposition_from_phext(...)` |
| `windowed_ic(z, x, ...)` → dict | `windowed_to_phext(result, ...)` → Book 5 windows + Book 2 verdict naming the max window's coordinate |
| `windowed_ic_envelope(...)` | each size's result is one `windowed_to_phext` volume (not wrapped further) |
| `simulations/*.csv` validation tables | `import_csv_to_phext` / `export_phext_to_csv` (byte-exact) |
| `balance_factor`, `ic_from_icc`, `mse_ratio_predicted`, `ic_threshold` | scalar helpers; `mse_ratio_predicted` is recorded in diagnose volumes (Book 2 S9); the rest are one-scroll writes if a user wants them — deliberately not wrapped, to avoid API surface the data doesn't need |

Module location: `cf_phext.load_cf()` imports the upstream module by
file path — explicit argument, `$CF_PHEXT_UPSTREAM`, a sibling
`circulatory_fidelity.py` (the proposed repo layout), or
`../upstream/python/` relative to the integration layout. No
installation and no vendored copy.

## 4. SQ usage (optional)

SQ v0.6.0 hosts a phext for coordinate queries over its CLI/REST:

```sh
cd <dir containing run.phext>
sq share run.phext          # starts the daemon (foreground)
sq toc                      # table of contents
sq select 1.1.1/1.1.4/2.4.1 # one relational entry, by address
sq shutdown
```

The integration test does exactly this against a generated matrix
phext and asserts SQ's `select` returns the same text as the built-in
reader. SQ is a bonus receipt only: the core path (module, demo,
tests) is pure Python and never requires it. Operational note,
evidenced with SQ v0.6.0: do not replace a hosted `.phext` file while
its daemon is running — shut the daemon down first.

**What is SQ-verified, exactly (2026-11-01, SQ v0.6.0):** for
coordinates whose components are all ≤ 999, `sq select` agrees with
the built-in reader — in the integration test (small matrix phext,
clean sequential run: PASS) and in manual probes against the full
simulations phext (values spot-checked across all 8 volumes).
**What is NOT SQ-verified:** coordinates with a component ≥ 1000 —
libphext saturates components at 999 (§5), so dataset rows ≥ 998 are
not SQ-addressable in this schema; and scripted SQ daemon startup on
a shared machine proved intermittent (the test reports SKIP with the
daemon's log when the daemon does not start; see §5). No claim in
this document rests on SQ behaving beyond that envelope.

## 5. Limitations (only those derived or evidenced)

- **Delimiter bytes cannot appear in scroll text.** Derived from the
  format definition: the 9 delimiter bytes are the structure. The
  writer raises `ValueError` if asked to store one. (None occur in
  CF's numeric/text data or the upstream CSVs — checked in tests.)
- **Whole-file rewrite on save.** This implementation serializes the
  full document per write (as does SQ on insert). Evidenced cost: the
  8-file simulations phext (33,965 data rows, 249,832 scrolls,
  3.6 MB) imports + exports in seconds; a workload of many small
  appends to one huge phext would feel the rewrite. A
  coordinate-addressed single-scroll update path is SQ's job, not
  this module's.
- **Provenance timestamps make hashes run-specific.** Each volume
  records a UTC timestamp, so two runs of the demo produce phexts
  with different SHA-256 hashes. Set `CF_PHEXT_TIMESTAMP` to any
  fixed string to make generation byte-reproducible.
- **CSV round-trip exactness is a property of these files' dialect**
  (QUOTE_MINIMAL-style quoting, LF endings, trailing newline), which
  the tests verify for all 8 upstream files. A CSV produced with a
  different dialect (CRLF, quote-all) would round-trip field-exact
  but not necessarily byte-exact. *Hypothesis (labelled):* normalizing
  dialect on import would close this; not implemented, not needed for
  the upstream corpus.
- **The `diagnose` precision lane uses Λ = corr⁻¹.** For a pair this
  is the natural 2×2 precision; for multivariate diagnoses a user
  should store their model precision via `matrix_to_phext` instead.
- **SQ/libphext caps coordinate components at 999 (evidenced).**
  On a synthetic phext with scrolls at Chapters 999/1000/1001, SQ's
  TOC lists one scroll (at 999, holding the last-written value) and
  selects at ≥1000 return empty. This schema's dataset lane puts row
  r at Chapter r+2, so SQ interop covers only the first ~997 data
  rows of each dataset volume. Neither `cf_phext` nor the reference
  dead-reckoning parser imposes this cap, and every non-SQ test
  covers the full range. Repacking rows to respect 999 would change
  existing coordinates; it is a question for the maintainer (§8),
  not a silent change.
- **Schema 1 files do not store the signed marginal correlation.**
  Only its magnitude (field 2) was recorded, so
  `coupling_from_phext` returns `marginal_corr=None` for schema 1
  rather than reconstructing a sign the file does not contain.
  (The conditional values *are* derivable for schema 1, from the
  stored matrix values — see §8.)
- **Point queries are full parses in `cf_phext` (measured, §9).**
  Addressability is not yet random access; an index changes that
  (prototype results in §9), but `cf_phext` does not ship one.
- **No native Julia phext implementation yet.** Julia results enter
  via CSV today (PHEXT.md shows how). The format spec needed to write
  one is §1–2 of this document plus the in-file manifest.

## 6. PR notes (what would land upstream)

Additive only — **no upstream files modified**:

- `python/cf_phext.py` (this module; its upstream-module lookup
  finds the sibling `circulatory_fidelity.py` first)
- `python/test_cf_phext.py` (the SQ test self-skips, with its
  daemon's log, if `sq` is absent or cannot start; the
  reference-parser agreement test self-skips if `parse_phext.py`
  isn't available outside the phext project)
- `python/demo.py` (the end-to-end demo, repo layout)
- `python/fixtures/schema1_chain.phext` (a schema-1 file generated
  by cf_phext 1.0.0, used by the backward-compatibility test)
- `verify.sh` at repo root (upstream suite + integration tests +
  demo, repo-relative, with an SQ preflight line)
- `PHEXT.md` at repo root (short README-style section + quickstart)
- `INTEGRATION.md` at repo root (this document)

`CITATION.cff`, `LICENSE`, the paper, Julia sources, notebooks and
figures are untouched, and CF's authorship/citation are credited, not
diluted: PHEXT.md closes by citing the CF paper (DOI
10.5281/zenodo.18121821). A pre-assembled tree with exactly this
layout is in the integration bundle's `repo/` directory.

## 7. Verification

```sh
bash verify.sh
```

runs upstream's own `test_circulatory_fidelity.py`, then
`test_cf_phext.py` (14 tests), then `demo.py`, then prints SHA-256
hashes of the generated phexts. Fresh-copy cost: the three commands
above, nothing to install beyond CF's existing numpy/scipy. The repo
copy has its own `verify.sh` with the same contents in repo layout;
both print an SQ preflight line and count PASS / SKIP / FAIL
separately (a skip is not a pass).

## 8. Marginal vs conditional coupling (schema 2)

Schema 1 stored, in each matrix-entry record, a *value* from one
matrix alongside an *IC* computed from a different (associated)
matrix, without saying so. For kind=precision the value is an entry
of Λ — a **conditional**-dependence quantity (in a Gaussian
graphical model, Λ_ij = 0 iff variables i, j are conditionally
independent given all the rest) — while the `ic` field was
|corr_ij| of the implied covariance Λ⁻¹ — a **marginal** quantity
that includes dependence mediated by other variables. Neither is
wrong; they answer different questions, and a record that mixes them
unlabelled invites misreading.

The canonical example is a three-variable Gaussian chain A–B–C with

```
Σ = [[1, .6, .36], [.6, 1, .6], [.36, .6, 1]]   (marginal correlations)
Λ = Σ⁻¹ = [[1.5625, -.9375, 0], [-.9375, 2.125, -.9375], [0, -.9375, 1.5625]]
```

- corr(A, C) = **0.36** — yet A and C are conditionally independent
  given B: Λ_AC = 0, so the partial correlation
  ρ_AC·B = −Λ_AC/√(Λ_AA·Λ_CC) = **0**. The 0.36 is entirely the
  indirect path through B (0.6 × 0.6).
- ρ_AB·C = 0.9375/√(1.5625·2.125) = **0.5144957554…**, which differs
  from the marginal corr(A, B) = 0.6.
- In two variables the two notions coincide (there is nothing else
  to condition on): for ρ = 0.8 both are 0.8. The distinction only
  exists in dimension ≥ 3 — which is exactly where a stored matrix
  stops being "a pair".

**Schema 2 (additive, backward-compatible)** makes both explicit
per entry record: field 2 keeps its schema-1 value and is now
documented as `ic_marginal`; field 7 adds `ic_conditional`
= |partial_ij|; field 8 adds the signed `partial_corr`; field 9 adds
the signed `marginal_corr`. Every volume stamps its schema version
at Book 7 Scroll 6 (`2`); a missing stamp means schema 1.
`coupling_from_phext(path, volume, i, j)` returns both couplings for
either schema — for schema 1 files it *derives* the conditional
values from the stored matrix values (directly from Λ for
kind=precision; via one read-time inverse otherwise) and returns
`marginal_corr=None`, because schema 1 never stored the sign and we
do not invent data. A schema-1 fixture generated by the original
code (`fixtures/schema1_chain.phext`, the chain's Λ) is in the test
suite as the compatibility evidence.

Why field 2 stays marginal: **upstream CF's pairwise IC is a
marginal quantity by design** — `inference_coupling(x, y)` is a
copula correlation of a pair of sample vectors; there is no
conditioning set. Redefining field 2 to a conditional quantity
would silently change CF's diagnostic, which this integration must
not do. Conditional coupling is offered as additional, labelled
information.

Two cautions, stated plainly: partial correlation is still a
*correlation* — neither marginal nor conditional dependence implies
causation, and a zero partial correlation is conditional
independence only under the Gaussian (graphical-model) reading.
And a question we explicitly leave to Aaron: **whether conditional
IC deserves to be a first-class CF diagnostic** (e.g. an
`inference_coupling_conditional` in upstream, with its own
thresholds and validation) is a scientific decision for the
maintainer. This PR only makes the stored semantics explicit; it
does not propose new CF theory.

**Cost (measured, `perf_coupling.py`; medians, numpy 1.26.4, this
machine).** One of the two paths needs a dense inverse; which one
depends on the input kind, symmetrically. For kind=precision,
marginal = inverse + normalize, conditional = normalize only (the
`partial_correlations(Λ, 'precision')` path provably performs no
inverse — there is a regression test that monkeypatches
`np.linalg.inv` to raise). For kind=correlation, the costs swap.

| d | precision: marginal | precision: conditional | correlation: marginal | correlation: conditional |
|--:|--------------------:|-----------------------:|----------------------:|-------------------------:|
| 10 | 0.014 ms | 0.009 ms | 0.001 ms | 0.015 ms |
| 50 | 0.154 ms | 0.013 ms | 0.001 ms | 0.157 ms |
| 100 | 0.996 ms | 0.024 ms | 0.002 ms | 1.004 ms |
| 200 | 7.384 ms | 0.088 ms | 0.008 ms | 7.385 ms |
| 400 | 57.424 ms | 0.336 ms | 0.059 ms | 56.835 ms |

Honest reading: below d ≈ 50 every path is under ~0.16 ms and the
distinction is cost-irrelevant; at d = 400 the inverse-carrying path
costs ~57 ms against ~0.3 ms (precision inputs) — material if a
caller needs *only* conditional coupling of a large precision
matrix, which is why the no-inverse entry points exist. Inside
`_write_matrix_lanes` both couplings are written and the inverse is
computed exactly once per call, as in schema 1.

## 9. Storage performance evidence (summary)

Full method and numbers: `BENCHMARK.md` in the integration bundle
(workspace-only; not part of the PR). Headlines from the final run,
on the simulations phext (249,832 scroll slots, identical seeded
255-query set, all answers verified): serialized size is a wash
(phext 3,634,939 B vs CSV 3,630,327 B, 1.001×); the current
`read_scroll` costs a full parse per query (median 385 ms) versus
CSV parse-per-query 12 ms; a prototype coordinate→(offset, length)
index built in 1.97 s (~64 MB traced) answers the same queries in a
median of 8.6 µs; full scans still favour bulk parsing (CSV 0.15 s,
phext load+walk 0.45 s); a one-field update costs a whole-file
rewrite in both (CSV 10 ms for its 0.76 MB file; phext 1.12 s for
3.63 MB). Conclusion recorded there: addressability is not yet
random access in `cf_phext`; the indexed-reader follow-up is
justified by measurement, and no speed claim is made for the
current implementation.
