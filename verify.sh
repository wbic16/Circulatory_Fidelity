#!/bin/sh
# verify.sh — Circulatory Fidelity x phext verification, REPO layout.
# Run from anywhere:  bash verify.sh        (repo root)
# Runs upstream's own test suite (python/test_circulatory_fidelity.py),
# the phext integration tests (python/test_cf_phext.py), and the demo
# (python/demo.py) — all repo-relative. Prints PASS/FAIL per check and
# SHA-256 hashes of the generated phexts. Non-zero exit on any failure.
# Tests that cannot run in this environment report SKIP with a reason
# (counted separately; a skip is not a pass).
set -u
cd "$(dirname "$0")"
PY="${PYTHON:-python3}"
FAIL=0

echo "== preflight: SQ availability =="
SQBIN="${SQ_BIN:-}"
if [ -z "$SQBIN" ]; then SQBIN="$(command -v sq || true)"; fi
if [ -z "$SQBIN" ] && [ -x "$HOME/.cargo/bin/sq" ]; then SQBIN="$HOME/.cargo/bin/sq"; fi
if [ -n "$SQBIN" ]; then
  SQPROBE_DIR="$(mktemp -d)"
  SQPROBE="$(cd "$SQPROBE_DIR" && "$SQBIN" --version 2>&1 | head -1)"
  rm -rf "$SQPROBE_DIR"
  echo "sq found at $SQBIN (probe in empty dir: ${SQPROBE:-no output})"
else
  echo "sq not found — SQ interop test will be SKIPPED (SQ interoperability UNVERIFIED in this environment)"
fi

step() {
  desc="$1"; shift
  echo ""
  echo "== $desc =="
  if "$@"; then echo "PASS: $desc"; else echo "FAIL: $desc"; FAIL=1; fi
}

step "upstream test suite (python/test_circulatory_fidelity.py)" \
  sh -c "cd python && $PY test_circulatory_fidelity.py"

step "integration tests (python/test_cf_phext.py)" \
  sh -c "cd python && $PY test_cf_phext.py"

step "end-to-end demo (python/demo.py)" \
  sh -c "cd python && $PY demo.py"

echo ""
echo "== SHA-256 of generated demo phexts =="
if command -v sha256sum >/dev/null 2>&1; then
  sha256sum python/output/*.phext || FAIL=1
else
  echo "sha256sum not available"; FAIL=1
fi

echo ""
if [ "$FAIL" -eq 0 ]; then
  echo "------------------------------------------------------------"
  echo "All checks passed."
else
  echo "VERIFICATION FAILED — see FAIL lines above."
fi
exit $FAIL
