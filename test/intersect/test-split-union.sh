#!/usr/bin/env bash
# Regression tests for the union used by intersect -split -f.
set -eu
BT=${BT:-../../bin/bedtools}
[ -x "$BT" ] || { echo "bedtools executable not found: $BT" >&2; exit 2; }
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
failures=0
tests=0
printf 'chr1\t0\t100\tquery\t0\t+\t0\t100\t0\t1\t100\t0\n' > "$tmp/a"

check_case() {
    local name=$1 fraction=$2 expected=$3 mode
    for mode in normal reversed sorted; do
        local options=()
        case "$mode" in
            normal) cp "$tmp/b" "$tmp/input" ;;
            reversed) awk '{ lines[NR]=$0 } END { for (i=NR;i>0;i--) print lines[i] }' "$tmp/b" > "$tmp/input" ;;
            sorted) LC_ALL=C sort -k1,1 -k2,2n "$tmp/b" > "$tmp/input"; options=(-sorted) ;;
        esac
        "$BT" intersect -a "$tmp/a" -b "$tmp/input" -split -f "$fraction" -u "${options[@]}" > "$tmp/observed"
        if [ "$expected" = 1 ]; then cp "$tmp/a" "$tmp/expected"; else : > "$tmp/expected"; fi
        tests=$((tests + 1))
        if diff -u "$tmp/expected" "$tmp/observed"; then
            echo "PASS: $name ($mode)"
        else
            echo "FAIL: $name ($mode)"
            failures=$((failures + 1))
        fi
    done
}

printf 'chr1\t0\t100\n' > "$tmp/b"
check_case single 0.9 1
printf 'chr1\t0\t100\nchr1\t40\t60\n' > "$tmp/b"
check_case nested 0.9 1
printf 'chr1\t0\t60\nchr1\t40\t100\n' > "$tmp/b"
check_case partial 0.9 1
printf 'chr1\t0\t30\nchr1\t20\t60\nchr1\t50\t100\n' > "$tmp/b"
check_case chain 0.9 1
printf 'chr1\t0\t100\nchr1\t0\t100\n' > "$tmp/b"
check_case duplicate 1.0 1
printf 'chr1\t0\t40\nchr1\t40\t100\n' > "$tmp/b"
check_case bookended 1.0 1
printf 'chr1\t0\t40\nchr1\t60\t100\n' > "$tmp/b"
check_case disjoint-covered 0.8 1
check_case disjoint-gap 0.81 0
printf 'chr1\t0\t40\nchr1\t10\t20\nchr1\t60\t100\n' > "$tmp/b"
check_case nested-disjoint 0.8 1
printf 'chr1\t120\t140\n' > "$tmp/b"
check_case no-overlap 0.1 0
printf 'chr1\t0\t100\tquery\t0\t+\t0\t100\t0\t2\t40,40\t0,60\n' > "$tmp/a"
printf 'chr1\t0\t100\nchr1\t10\t20\n' > "$tmp/b"
check_case split-query 1.0 1
printf 'chr1\t40\t60\n' > "$tmp/b"
check_case intron-only 0.1 0
printf 'chr1\t1000\t1100\tquery\t0\t+\t1000\t1100\t0\t1\t100\t0\n' > "$tmp/a"
printf 'chr1\t1000\t1100\nchr1\t1040\t1060\n' > "$tmp/b"
check_case coordinate-offset 0.9 1
printf 'split-union: %s checks, %s failures\n' "$tests" "$failures"
[ "$failures" -eq 0 ]
