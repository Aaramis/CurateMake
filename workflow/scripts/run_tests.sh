#!/usr/bin/env bash
# Run every test file in the repository and fail if any of them fails.
#
# The test files are script-style: each prints "PASSED n/m" and exits non-zero on
# failure, so they run under a plain interpreter and need no test framework. pytest
# cannot collect them as a directory, because they call sys.exit() at module level.
# They are all offline: the CoL client is faked, so nothing here touches ChecklistBank
# and no result can drift with a CoL release.
#
# Usage: bash workflow/scripts/run_tests.sh [interpreter]
set -uo pipefail
PY="${1:-python3}"
cd "$(cd "$(dirname "$0")/../.." && pwd)"

TESTS=(
  workflow/scripts/test_col_taxonomic_matcher.py
  workflow/scripts/test_sativa_taxonomy_report.py
  workflow/scripts/utils/test_retry_policy.py
  workflow/scripts/utils/test_taxonomy_normalize.py
  workflow/scripts/MSA/ehmm/test_build_ehmm_msa.py
  workflow/scripts/Sativa/test_mis_confidence_filter.py
  workflow/scripts/Sativa/test_auto_sativa.py
  workflow/scripts/Sativa/test_classify_modification.py
  workflow/scripts/Sativa/test_vote_sativa.py
)

rc=0
for t in "${TESTS[@]}"; do
  out="$("$PY" "$t" 2>&1)"; code=$?
  line="$(printf '%s\n' "$out" | grep -E '^PASSED|^FAILED' | tail -1)"
  printf '%-44s %s\n' "$(basename "$t")" "${line:-aucune ligne de resultat}"
  if [ "$code" -ne 0 ]; then
    rc=1
    printf '%s\n' "$out" | tail -25
  fi
done

if [ "$rc" -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit "$rc"
