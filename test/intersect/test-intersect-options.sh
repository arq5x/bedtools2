. ./helpers-intersect.sh options
set -e

# Output modes and overlap rules of bedtools intersect.
#
# Expected values come from two sources, never from running bedtools:
#   hand       small inputs worked out by hand (literal expected files).
#   reference  an awk model of the documented semantics (ref.awk below),
#              applied to deterministic generated inputs (fixed LCG).
#
# Hit order inside one A record, which the output shows:
#   unsorted mode: bins are visited from the smallest level (16 kb bins) to
#                  the largest, by ascending bin number within a level, and
#                  in B file order within a bin;
#   -sorted mode : B file order.
# Fixtures that keep every B inside one 16 kb bin therefore give B file
# order in both modes, and are checked in both modes with one expectation.

EXP=$TMPD/exp
OBS=$TMPD/obs
N=0

nextid() { N=$((N + 1)); ID=$(printf '%03d' $N); }

# tbegin <mode> <args...>: begin the next id and log it for the targets table.
tbegin()
{
	local id mode=$1
	shift
	nextid
	id=$ID
	begin $id
	if [ -n "${TARGETS_OUT-}" ]; then
		printf '%s\t%s\t%s\t%s\n' $id "$* $mode" "$CLASS" "$DERIV" >> "$TARGETS_OUT"
	fi
}

# ck <expected file> <options...>: unsorted run then -sorted run on $A, $B.
ck()
{
	local exp=$1 m
	shift
	for m in unsorted -sorted; do
		tbegin $m "$@"
		if [ $m = unsorted ]; then
			$BT intersect "$@" -a $A -b $B > $OBS
		else
			$BT intersect "$@" -sorted -a $A -b $B > $OBS
		fi
		check $OBS $exp
	done
}

# cku / cks: only the unsorted run / only the -sorted run.
cku()
{
	local exp=$1
	shift
	tbegin unsorted "$@"
	$BT intersect "$@" -a $A -b $B > $OBS
	check $OBS $exp
}

cks()
{
	local exp=$1
	shift
	tbegin -sorted "$@"
	$BT intersect "$@" -sorted -a $A -b $B > $OBS
	check $OBS $exp
}

# ---------------------------------------------------------------------------
# Reference model (documented semantics of the overlap rules and modes).
#   hit: same chromosome; positive-length overlap; -s: both strands known
#   and equal; -S: both known and different; with -f/-F the overlap divided
#   by the length of A / of B must be >= f / F, -r sets F=f, -e accepts
#   either fraction instead of both.
#   modes: see the print section.  Variables: opts (the bedtools options),
#   order=bin (unsorted hit order) or file (-sorted hit order).
# ---------------------------------------------------------------------------
cat > $TMPD/ref.awk <<'EOF'
function bin(s, e,    i, sh, a, b) {
	e = e - 1
	sh = 16384
	for (i = 0; i < 8; i++) {
		a = int(s / sh); b = int(e / sh)
		if (a == b) return i * 1000000 + a
		sh = sh * 8
	}
	return -1
}
function nullrec(nf) {
	if (nf == 3) return ".\t-1\t-1"
	if (nf == 4) return ".\t-1\t-1\t."
	if (nf == 5) return ".\t-1\t-1\t.\t-1"
	return ".\t-1\t-1\t.\t-1\t."
}
BEGIN {
	FS = OFS = "\t"
	n = split(opts, o, " ")
	for (i = 1; i <= n; i++) {
		t = o[i]
		if (t == "-f") { f = o[i + 1] + 0; i++ }
		else if (t == "-F") { F = o[i + 1] + 0; i++ }
		else if (t == "-r") r = 1
		else if (t == "-e") e = 1
		else if (t == "-s") s = 1
		else if (t == "-S") S = 1
		else if (t == "-wa") wa = 1
		else if (t == "-wb") wb = 1
		else if (t == "-loj") loj = 1
		else if (t == "-wo") wo = 1
		else if (t == "-wao") { wao = 1; wo = 1 }
		else if (t == "-u") u = 1
		else if (t == "-c") c = 1
		else if (t == "-v") v = 1
	}
	if (r) F = f
}
END {
	if (shape) { c = 0; for (l in lvl) c++; print c }
}
FNR == NR {
	nb++
	bl[nb] = $0; bs[nb] = $2 + 0; be[nb] = $3 + 0; bst[nb] = $6
	bk[nb] = bin($2 + 0, $3 + 0)
	lvl[int(bk[nb] / 1000000)] = 1
	k = ++cnt[$1]
	idx[$1, k] = nb
	bnf = NF
	next
}
{
	as = $2 + 0; ae = $3 + 0; alen = ae - as
	nh = 0
	for (k = 1; k <= cnt[$1]; k++) {
		j = idx[$1, k]
		lo = (as > bs[j]) ? as : bs[j]
		hi = (ae < be[j]) ? ae : be[j]
		ov = hi - lo
		if (ov <= 0) continue
		if (s && !($6 != "." && $6 == bst[j])) continue
		if (S && !($6 != "." && bst[j] != "." && $6 != bst[j])) continue
		if (f > 0 || F > 0) {
			fa = ov / alen; fb = ov / (be[j] - bs[j])
			if (e) { if (!(fa >= f || fb >= F)) continue }
			else if (!(fa >= f && fb >= F)) continue
		}
		nh++; hit[nh] = j; hov[nh] = ov; hlo[nh] = lo; hhi[nh] = hi
	}
	if (order == "bin") {
		for (p = 2; p <= nh; p++) {
			j = hit[p]; q = hov[p]; l = hlo[p]; h = hhi[p]
			m = p - 1
			while (m >= 1 && (bk[hit[m]] > bk[j] || (bk[hit[m]] == bk[j] && hit[m] > j))) {
				hit[m + 1] = hit[m]; hov[m + 1] = hov[m]; hlo[m + 1] = hlo[m]; hhi[m + 1] = hhi[m]
				m--
			}
			hit[m + 1] = j; hov[m + 1] = q; hlo[m + 1] = l; hhi[m + 1] = h
		}
	}
	if (u) { if (nh > 0) print $0; next }
	if (v) { if (nh == 0) print $0; next }
	if (c) { print $0, nh; next }
	if (nh == 0) {
		if (wao) print $0, nullrec(bnf), 0
		else if (loj) print $0, nullrec(bnf)
		next
	}
	for (p = 1; p <= nh; p++) {
		if (!wa && !wb && !wo && !loj) {
			line = $1 OFS hlo[p] OFS hhi[p]
			for (q = 4; q <= NF; q++) line = line OFS $q
			print line
		} else if ((wa && wb) || loj) print $0, bl[hit[p]]
		else if (wa) print $0
		else if (wb) {
			line = $1 OFS hlo[p] OFS hhi[p]
			for (q = 4; q <= NF; q++) line = line OFS $q
			print line, bl[hit[p]]
		} else print $0, bl[hit[p]], hov[p]
	}
}
EOF

cat > $TMPD/gen.awk <<'EOF'
# Deterministic BED6 generator: Park-Miller LCG, positive lengths.
function rnd() { x = (x * 16807) % 2147483647; return x / 2147483647 }
BEGIN {
	x = seed
	for (i = 0; i < 5; i++) rnd()
	for (i = 1; i <= n; i++) {
		c = 1 + int(rnd() * nchr)
		u = rnd()
		len = 1 + int(u * u * u * maxlen)
		st = int(rnd() * (maxc - len + 1))
		sd = rnd()
		strand = (sd < 0.4) ? "+" : ((sd < 0.8) ? "-" : ".")
		printf "chr%d\t%d\t%d\tn%d\t%d\t%s\n", c, st, st + len, i, int(rnd() * 1000), strand
	}
}
EOF

# grp: stdin to stdout with the lines of each A record (keyed by field 4, the
# unique A name, in order of first appearance) sorted among themselves, so
# the order of hits within one A record is not compared.
grp()
{
	awk 'BEGIN { FS = OFS = "\t" } { if (!($4 in g)) g[$4] = ++n; print g[$4], $0 }' \
		| sort -t "$(printf '\t')" -k1,1n -k2 | cut -f2-
}

# refck <A> <B> <order> <bedtools options...>: compare one run with ref.awk.
# order is unsorted or sorted (-sorted, B file order). Unsorted hit order
# across bin levels is not documented, so unsorted runs compare the hits of
# each A record order-insensitively (A record order is still compared);
# sorted runs compare the exact order.
refck()
{
	local a=$1 b=$2 ord=$3
	shift 3
	tbegin $ord "$@"
	if [ $ord = unsorted ]; then
		$BT intersect "$@" -a $a -b $b | grp > $OBS
		awk -v opts="$*" -v order=file -f $TMPD/ref.awk $b $a | grp > $EXP
	else
		$BT intersect "$@" -sorted -a $a -b $b > $OBS
		awk -v opts="$*" -v order=file -f $TMPD/ref.awk $b $a > $EXP
	fi
	check $OBS $EXP
}

# ===========================================================================
## fixture F1: five A records, nine B records, all inside one 16 kb bin
CLASS="hand F1: BED6, one bin, abutting, duplicate B, absent chromosomes"
DERIV=hand
A=$TMPD/f1a.bed
B=$TMPD/f1b.bed
rows > $A <<'EOF'
chr1 100 200 a1 1 +
chr1 300 400 a2 2 -
chr1 500 600 a3 3 .
chr1 700 800 a4 4 +
chr2 10 20 a5 5 +
EOF
rows > $B <<'EOF'
chr1 100 200 b3 12 -
chr1 150 250 b1 10 +
chr1 150 250 b1 10 +
chr1 190 195 b2 11 -
chr1 300 400 b4 13 +
chr1 350 700 b5 14 .
chr1 590 650 b7 16 -
chr1 600 650 b6 15 +
chr3 1 2 b8 17 +
EOF
A1="chr1 100 200 a1 1 +"
A2="chr1 300 400 a2 2 -"
A3="chr1 500 600 a3 3 ."
A4="chr1 700 800 a4 4 +"
A5="chr2 10 20 a5 5 +"
B3="chr1 100 200 b3 12 -"
B1="chr1 150 250 b1 10 +"
B2="chr1 190 195 b2 11 -"
B4="chr1 300 400 b4 13 +"
B5="chr1 350 700 b5 14 ."
B7="chr1 590 650 b7 16 -"
NUL=". -1 -1 . -1 ."

# Default output: the overlapping part of A; hits in B order; A3 hit by B5 and
# B7 but not by B6 (abuts at 600); A4 abuts B5 at 700.
rows > $EXP <<'EOF'
chr1 100 200 a1 1 +
chr1 150 200 a1 1 +
chr1 150 200 a1 1 +
chr1 190 195 a1 1 +
chr1 300 400 a2 2 -
chr1 350 400 a2 2 -
chr1 500 600 a3 3 .
chr1 590 600 a3 3 .
EOF
ck $EXP

# -wa repeats the whole A record once per hit.
{ for i in 1 2 3 4; do echo "$A1"; done; echo "$A2"; echo "$A2"; echo "$A3"; echo "$A3"; } | rows > $EXP
ck $EXP -wa

# -wb: overlapping part of A followed by the whole B record.
{
echo "chr1 100 200 a1 1 + $B3"
echo "chr1 150 200 a1 1 + $B1"
echo "chr1 150 200 a1 1 + $B1"
echo "chr1 190 195 a1 1 + $B2"
echo "chr1 300 400 a2 2 - $B4"
echo "chr1 350 400 a2 2 - $B5"
echo "chr1 500 600 a3 3 . $B5"
echo "chr1 590 600 a3 3 . $B7"
} | rows > $EXP
ck $EXP -wb

# -wa -wb: whole A then whole B per hit.
{
echo "$A1 $B3"; echo "$A1 $B1"; echo "$A1 $B1"; echo "$A1 $B2"
echo "$A2 $B4"; echo "$A2 $B5"
echo "$A3 $B5"; echo "$A3 $B7"
} | rows > $EXP
ck $EXP -wa -wb

# -wo: whole A, whole B, overlap in bases; no line for A without a hit.
{
echo "$A1 $B3 100"; echo "$A1 $B1 50"; echo "$A1 $B1 50"; echo "$A1 $B2 5"
echo "$A2 $B4 100"; echo "$A2 $B5 50"
echo "$A3 $B5 100"; echo "$A3 $B7 10"
} | rows > $EXP
ck $EXP -wo

# -wao: as -wo, plus A without a hit with the BED6 NULL record and 0.
{
echo "$A1 $B3 100"; echo "$A1 $B1 50"; echo "$A1 $B1 50"; echo "$A1 $B2 5"
echo "$A2 $B4 100"; echo "$A2 $B5 50"
echo "$A3 $B5 100"; echo "$A3 $B7 10"
echo "$A4 $NUL 0"; echo "$A5 $NUL 0"
} | rows > $EXP
ck $EXP -wao

# -loj: whole A and B per hit; A without a hit gets the NULL record only.
{
echo "$A1 $B3"; echo "$A1 $B1"; echo "$A1 $B1"; echo "$A1 $B2"
echo "$A2 $B4"; echo "$A2 $B5"
echo "$A3 $B5"; echo "$A3 $B7"
echo "$A4 $NUL"; echo "$A5 $NUL"
} | rows > $EXP
ck $EXP -loj

# -u reports each A with at least one hit once.
{ echo "$A1"; echo "$A2"; echo "$A3"; } | rows > $EXP
ck $EXP -u

# -c reports the hit count (duplicate B counted twice), 0 when none.
{ echo "$A1 4"; echo "$A2 2"; echo "$A3 2"; echo "$A4 0"; echo "$A5 0"; } | rows > $EXP
ck $EXP -c

# -v reports only the A records without a hit, including absent chromosome.
{ echo "$A4"; echo "$A5"; } | rows > $EXP
ck $EXP -v

# -s with -wo: only same-strand B. A1(+) keeps the two B1(+); A2(-) has B4(+)
# and B5(.) so nothing; A3(.) never matches.
{ echo "$A1 $B1 50"; echo "$A1 $B1 50"; } | rows > $EXP
ck $EXP -s -wo

# -S with -wo: different strand needs both known: A1(+) with B3(-), B2(-);
# A2(-) with B4(+); A3(.) none; B5(.) never counts.
{ echo "$A1 $B3 100"; echo "$A1 $B2 5"; echo "$A2 $B4 100"; } | rows > $EXP
ck $EXP -S -wo

# -s -c: zero for everything but A1.
{ echo "$A1 2"; echo "$A2 0"; echo "$A3 0"; echo "$A4 0"; echo "$A5 0"; } | rows > $EXP
ck $EXP -s -c

# -S -c
{ echo "$A1 2"; echo "$A2 1"; echo "$A3 0"; echo "$A4 0"; echo "$A5 0"; } | rows > $EXP
ck $EXP -S -c

# -s -v: everything except A1.
{ echo "$A2"; echo "$A3"; echo "$A4"; echo "$A5"; } | rows > $EXP
ck $EXP -s -v

# -S -v: A1 and A2 have a different-strand hit.
{ echo "$A3"; echo "$A4"; echo "$A5"; } | rows > $EXP
ck $EXP -S -v

# -s -wao: unmatched A get the NULL record with 0.
{ echo "$A1 $B1 50"; echo "$A1 $B1 50"; echo "$A2 $NUL 0"; echo "$A3 $NUL 0"; echo "$A4 $NUL 0"; echo "$A5 $NUL 0"; } | rows > $EXP
ck $EXP -s -wao

# -S -loj
{ echo "$A1 $B3"; echo "$A1 $B2"; echo "$A2 $B4"; echo "$A3 $NUL"; echo "$A4 $NUL"; echo "$A5 $NUL"; } | rows > $EXP
ck $EXP -S -loj

# -f 0.5 -wo: fraction of A. A1 keeps B3 (1.0), B1 (0.5, boundary), B1; drops
# B2 (0.05). A2: B4 (1.0), B5 (0.5). A3: B5 (1.0), B7 (0.1) dropped.
{
echo "$A1 $B3 100"; echo "$A1 $B1 50"; echo "$A1 $B1 50"
echo "$A2 $B4 100"; echo "$A2 $B5 50"
echo "$A3 $B5 100"
} | rows > $EXP
ck $EXP -f 0.5 -wo

# -F 1 -wao: B must lie entirely inside A. B3 and B2 for A1, B4 for A2.
{
echo "$A1 $B3 100"; echo "$A1 $B2 5"
echo "$A2 $B4 100"
echo "$A3 $NUL 0"; echo "$A4 $NUL 0"; echo "$A5 $NUL 0"
} | rows > $EXP
ck $EXP -F 1 -wao

# -f 0.5 -F 1: both fractions. Only identical-or-contained pairs: A1-B3, A2-B4.
{ echo "$A1 $B3 100"; echo "$A2 $B4 100"; } | rows > $EXP
ck $EXP -f 0.5 -F 1 -wo

# -f 0.5 -F 1 -e: either fraction. Adds A1-B1 twice and A1-B2 (B fully
# inside), A2-B5 (0.5 of A2) and A3-B5 (all of A3); A3-B7 fails both.
{
echo "$A1 $B3 100"; echo "$A1 $B1 50"; echo "$A1 $B1 50"; echo "$A1 $B2 5"
echo "$A2 $B4 100"; echo "$A2 $B5 50"
echo "$A3 $B5 100"
} | rows > $EXP
ck $EXP -f 0.5 -F 1 -e -wo

# -f 0.5 -r -wo: both at least 0.5. A1-B3 (1,1), A1-B1 twice (0.5,0.5), A2-B4.
{ echo "$A1 $B3 100"; echo "$A1 $B1 50"; echo "$A1 $B1 50"; echo "$A2 $B4 100"; } | rows > $EXP
ck $EXP -f 0.5 -r -wo

# -f 0.5 -v: A without a hit passing the fraction. All five A with a hit
# pass for A1-A3, so only A4 and A5 remain.
{ echo "$A4"; echo "$A5"; } | rows > $EXP
ck $EXP -f 0.5 -v

# -F 1 -v: A3 has no B inside it, so it joins A4 and A5.
{ echo "$A3"; echo "$A4"; echo "$A5"; } | rows > $EXP
ck $EXP -F 1 -v

# -F 1 -u: A1 (B3, B2) and A2 (B4).
{ echo "$A1"; echo "$A2"; } | rows > $EXP
ck $EXP -F 1 -u

# -s -f 0.5 -wo: same strand and A fraction >= 0.5: A1-B1 twice.
{ echo "$A1 $B1 50"; echo "$A1 $B1 50"; } | rows > $EXP
ck $EXP -s -f 0.5 -wo

# -S -F 1 -c: different strand and B inside A: A1 has B3 and B2, A2 has B4.
{ echo "$A1 2"; echo "$A2 1"; echo "$A3 0"; echo "$A4 0"; echo "$A5 0"; } | rows > $EXP
ck $EXP -S -F 1 -c

# -wa with -f: A repeated per surviving hit (A1: B3, B1, B1).
{ echo "$A1"; echo "$A1"; echo "$A1"; echo "$A2"; echo "$A2"; echo "$A3"; } | rows > $EXP
ck $EXP -wa -f 0.5

# -wb with -F 1: A overlap part then B for B fully inside A.
{ echo "chr1 100 200 a1 1 + $B3"; echo "chr1 190 195 a1 1 + $B2"; echo "chr1 300 400 a2 2 - $B4"; } | rows > $EXP
ck $EXP -wb -F 1

# ===========================================================================
## fixture F2: one A, six B with distinct overlap fractions
CLASS="hand F2: single A, fraction boundaries, exact 0.25/0.5/1"
A=$TMPD/f2a.bed
B=$TMPD/f2b.bed
rows > $A <<'EOF'
chr1 1000 1100 A 1 +
EOF
rows > $B <<'EOF'
chr1 900 1010 b3 3 +
chr1 990 1110 b6 6 -
chr1 1000 1050 b1 1 +
chr1 1000 1100 b7 7 .
chr1 1020 1060 b4 4 -
chr1 1075 1200 b2 2 +
chr1 1100 1150 b5 5 +
EOF
# Overlap (bases, fraction of A, fraction of B) of A=[1000,1100) with each B:
#   b3 [900,1010)   10  0.1     0.0909
#   b6 [990,1110)   100 1.0     0.8333
#   b1 [1000,1050)  50  0.5     1.0
#   b7 [1000,1100)  100 1.0     1.0
#   b4 [1020,1060)  40  0.4     1.0
#   b2 [1075,1200)  25  0.25    0.2
#   b5 [1100,1150)  abuts, never a hit
AL="chr1 1000 1100 A 1 +"
W3="chr1 900 1010 b3 3 + 10"
W6="chr1 990 1110 b6 6 - 100"
W1="chr1 1000 1050 b1 1 + 50"
W7="chr1 1000 1100 b7 7 . 100"
W4="chr1 1020 1060 b4 4 - 40"
W2="chr1 1075 1200 b2 2 + 25"

# fcase <options> <hit lines...>: -wo and -c with the listed hits (in B order).
fcase()
{
	local opts=$1 w
	shift
	for w in "$@"; do echo "$AL $w"; done | rows > $EXP
	ck $EXP -wo $opts
	echo "$AL $#" | rows > $EXP
	cku $EXP -c $opts
}

# No fraction: every B with a positive overlap; b5 abuts at 1100 and is out.
fcase "" "$W3" "$W6" "$W1" "$W7" "$W4" "$W2"
# -f 1.0 (the largest allowed): A fully covered, b6 and b7 only.
fcase "-f 1.0" "$W6" "$W7"
# -F 1.0: B fully inside A: b1, b7, b4.
fcase "-F 1.0" "$W1" "$W7" "$W4"
# -f 1 -F 1: only the identical interval.
fcase "-f 1 -F 1" "$W7"
# -f 1 -r: b6 covers all of A but only 0.83 of itself lies in A.
fcase "-f 1 -r" "$W7"
# -f 0.5 boundary: b1 is exactly 0.5 of A and is kept; b4 (0.4) is not.
fcase "-f 0.5" "$W6" "$W1" "$W7"
# -f 0.51 just above the boundary drops b1.
fcase "-f 0.51" "$W6" "$W7"
# -f 0.25 boundary: b2 is exactly 0.25 and is kept.
fcase "-f 0.25" "$W6" "$W1" "$W7" "$W4" "$W2"
# -f 0.26 drops b2.
fcase "-f 0.26" "$W6" "$W1" "$W7" "$W4"
# -f 0.1 boundary: b3 is exactly 10/100, kept.
fcase "-f 0.1" "$W3" "$W6" "$W1" "$W7" "$W4" "$W2"
# -f 0.11 drops b3.
fcase "-f 0.11" "$W6" "$W1" "$W7" "$W4" "$W2"
# -f 0.01 very small fraction keeps all hits.
fcase "-f 0.01" "$W3" "$W6" "$W1" "$W7" "$W4" "$W2"
# -F 0.2 boundary: b2 is exactly 25/125 of B and is kept; b3 (0.09) is not.
fcase "-F 0.2" "$W6" "$W1" "$W7" "$W4" "$W2"
# -F 0.21 drops b2.
fcase "-F 0.21" "$W6" "$W1" "$W7" "$W4"
# -f 0.25 -F 0.5 requires both: drops b2 (0.2 of B) and b3.
fcase "-f 0.25 -F 0.5" "$W6" "$W1" "$W7" "$W4"
# -f 0.25 -F 0.5 -e requires either: b2 passes on A (0.25).
fcase "-f 0.25 -F 0.5 -e" "$W6" "$W1" "$W7" "$W4" "$W2"
# -f 0.9 -F 0.9 (both): only b7.
fcase "-f 0.9 -F 0.9" "$W7"
# -f 0.9 -F 0.9 -e (either): b6 (A 1.0), b1 and b4 (B 1.0), b7.
fcase "-f 0.9 -F 0.9 -e" "$W6" "$W1" "$W7" "$W4"
# -f 0.6 -r: b6 (1.0, 0.83) and b7; b1 fails on A (0.5).
fcase "-f 0.6 -r" "$W6" "$W7"
# -f 0.9 -r: b6 fails on B (0.83), only b7.
fcase "-f 0.9 -r" "$W7"
# -f 0.4 -r boundary: b4 is 0.4 of A and 1.0 of B, kept; b2 (0.25, 0.2) not.
fcase "-f 0.4 -r" "$W6" "$W1" "$W7" "$W4"
# -f 0.5 -F 0.01 -e: either fraction; b3 passes on B (0.09), b2 on B (0.2);
# abutting b5 has 0 bases and passes neither.
fcase "-f 0.5 -F 0.01 -e" "$W3" "$W6" "$W1" "$W7" "$W4" "$W2"
# -s: A is +, same strand B are b3, b1, b2 (b7 is '.', never matches).
fcase "-s" "$W3" "$W1" "$W2"
# -S: opposite strand needs both known: b6 and b4.
fcase "-S" "$W6" "$W4"
# -s -f 0.25: b3 fails (0.1).
fcase "-s -f 0.25" "$W1" "$W2"
# -S -F 1: b4 is inside A on the other strand; b6 is not inside A.
fcase "-S -F 1" "$W4"

# -u and -v on F2: A has hits so -u prints it, -v prints nothing; with
# -f 1 -F 1 -s the only candidate b7 is strand '.', so -v prints A.
echo "$AL" | rows > $EXP
ck $EXP -u -f 0.5
: > $EXP
ck $EXP -v -f 0.5
echo "$AL" | rows > $EXP
ck $EXP -v -f 1 -F 1 -s
: > $EXP
ck $EXP -u -f 1 -F 1 -s

# -loj with -f 1 -F 1 -s: no surviving hit, so A and the NULL record.
echo "$AL $NUL" | rows > $EXP
ck $EXP -loj -f 1 -F 1 -s
# -wao with the same options prints the NULL record and 0.
echo "$AL $NUL 0" | rows > $EXP
ck $EXP -wao -f 1 -F 1 -s

# ===========================================================================
## fixture F3: zero-length records (a point p is the interval [p-1,p+1))
CLASS="hand F3: zero-length records, abutting"
DERIV="hand (source model: zero-length record treated as [p-1,p+1))"
# A zero-length A with nonzero B. a1 at 50: B [40,50) and [50,60) overlap
# the widened [49,51); [51,60) only touches. a2 at 100: B [90,101) and
# [100,110) hit, [101,120) touches. a3 at 200: [180,200) and [200,205) hit,
# [201,210) touches. a4 has no neighbour.
A=$TMPD/f3a.bed
B=$TMPD/f3b.bed
rows > $A <<'EOF'
chr1 50 50 a1
chr1 100 100 a2
chr1 200 200 a3
chr1 300 300 a4
chr2 5 5 a5
EOF
rows > $B <<'EOF'
chr1 40 50 b1
chr1 50 60 b2
chr1 51 60 b3
chr1 90 101 b4
chr1 100 110 b5
chr1 101 120 b6
chr1 180 200 b7
chr1 200 205 b8
chr1 201 210 b9
EOF
# Default: the zero-length A is printed unchanged once per hit.
rows > $EXP <<'EOF'
chr1 50 50 a1
chr1 50 50 a1
chr1 100 100 a2
chr1 100 100 a2
chr1 200 200 a3
chr1 200 200 a3
EOF
ck $EXP
# -wa
ck $EXP -wa
# -wb: A (unchanged) and B record.
rows > $EXP <<'EOF'
chr1 50 50 a1 chr1 40 50 b1
chr1 50 50 a1 chr1 50 60 b2
chr1 100 100 a2 chr1 90 101 b4
chr1 100 100 a2 chr1 100 110 b5
chr1 200 200 a3 chr1 180 200 b7
chr1 200 200 a3 chr1 200 205 b8
EOF
ck $EXP -wb
# -wo: the zero-length A overlaps by 0 bases.
rows > $EXP <<'EOF'
chr1 50 50 a1 chr1 40 50 b1 0
chr1 50 50 a1 chr1 50 60 b2 0
chr1 100 100 a2 chr1 90 101 b4 0
chr1 100 100 a2 chr1 100 110 b5 0
chr1 200 200 a3 chr1 180 200 b7 0
chr1 200 200 a3 chr1 200 205 b8 0
EOF
ck $EXP -wo
# -wao adds the BED4 NULL record for a4 and a5.
rows > $EXP <<'EOF'
chr1 50 50 a1 chr1 40 50 b1 0
chr1 50 50 a1 chr1 50 60 b2 0
chr1 100 100 a2 chr1 90 101 b4 0
chr1 100 100 a2 chr1 100 110 b5 0
chr1 200 200 a3 chr1 180 200 b7 0
chr1 200 200 a3 chr1 200 205 b8 0
chr1 300 300 a4 . -1 -1 . 0
chr2 5 5 a5 . -1 -1 . 0
EOF
ck $EXP -wao
# -loj
rows > $EXP <<'EOF'
chr1 50 50 a1 chr1 40 50 b1
chr1 50 50 a1 chr1 50 60 b2
chr1 100 100 a2 chr1 90 101 b4
chr1 100 100 a2 chr1 100 110 b5
chr1 200 200 a3 chr1 180 200 b7
chr1 200 200 a3 chr1 200 205 b8
chr1 300 300 a4 . -1 -1 .
chr2 5 5 a5 . -1 -1 .
EOF
ck $EXP -loj
# -c
rows > $EXP <<'EOF'
chr1 50 50 a1 2
chr1 100 100 a2 2
chr1 200 200 a3 2
chr1 300 300 a4 0
chr2 5 5 a5 0
EOF
ck $EXP -c
# -u
rows > $EXP <<'EOF'
chr1 50 50 a1
chr1 100 100 a2
chr1 200 200 a3
EOF
ck $EXP -u
# -v
rows > $EXP <<'EOF'
chr1 300 300 a4
chr2 5 5 a5
EOF
ck $EXP -v

# Nonzero A with zero-length B (point q is [q-1,q+1)). B points 99 (touches
# A1 at 100, no hit), 100, 150, 200 (hit A1 [100,200)), 201 (touches, no),
# 300, 399, 400 (hit A2 [300,400); 400 overlaps [399,400)), 401 (touches).
A=$TMPD/f3c.bed
B=$TMPD/f3d.bed
rows > $A <<'EOF'
chr1 100 200 a1 1 +
chr1 300 400 a2 2 -
chr1 500 600 a3 3 .
chr1 700 800 a4 4 +
EOF
rows > $B <<'EOF'
chr1 99 99 p1 1 +
chr1 100 100 p2 2 +
chr1 150 150 p3 3 -
chr1 200 200 p4 4 +
chr1 201 201 p5 5 +
chr1 300 300 p6 6 +
chr1 399 399 p7 7 +
chr1 400 400 p8 8 -
chr1 401 401 p9 9 +
EOF
# -wa
rows > $EXP <<'EOF'
chr1 100 200 a1 1 +
chr1 100 200 a1 1 +
chr1 100 200 a1 1 +
chr1 300 400 a2 2 -
chr1 300 400 a2 2 -
chr1 300 400 a2 2 -
EOF
ck $EXP -wa
# -wa -wb: the B point is printed as written.
rows > $EXP <<'EOF'
chr1 100 200 a1 1 + chr1 100 100 p2 2 +
chr1 100 200 a1 1 + chr1 150 150 p3 3 -
chr1 100 200 a1 1 + chr1 200 200 p4 4 +
chr1 300 400 a2 2 - chr1 300 300 p6 6 +
chr1 300 400 a2 2 - chr1 399 399 p7 7 +
chr1 300 400 a2 2 - chr1 400 400 p8 8 -
EOF
ck $EXP -wa -wb
# -loj
rows > $EXP <<'EOF'
chr1 100 200 a1 1 + chr1 100 100 p2 2 +
chr1 100 200 a1 1 + chr1 150 150 p3 3 -
chr1 100 200 a1 1 + chr1 200 200 p4 4 +
chr1 300 400 a2 2 - chr1 300 300 p6 6 +
chr1 300 400 a2 2 - chr1 399 399 p7 7 +
chr1 300 400 a2 2 - chr1 400 400 p8 8 -
chr1 500 600 a3 3 . . -1 -1 . -1 .
chr1 700 800 a4 4 + . -1 -1 . -1 .
EOF
ck $EXP -loj
# -c
rows > $EXP <<'EOF'
chr1 100 200 a1 1 + 3
chr1 300 400 a2 2 - 3
chr1 500 600 a3 3 . 0
chr1 700 800 a4 4 + 0
EOF
ck $EXP -c
# -u and -v
rows > $EXP <<'EOF'
chr1 100 200 a1 1 +
chr1 300 400 a2 2 -
EOF
ck $EXP -u
rows > $EXP <<'EOF'
chr1 500 600 a3 3 .
chr1 700 800 a4 4 +
EOF
ck $EXP -v
# -s: A1(+) with + points p2, p4; A2(-) with - point p8.
rows > $EXP <<'EOF'
chr1 100 200 a1 1 + 2
chr1 300 400 a2 2 - 1
chr1 500 600 a3 3 . 0
chr1 700 800 a4 4 + 0
EOF
ck $EXP -s -c
# -S: A1(+) with - point p3; A2(-) with + points p6, p7.
rows > $EXP <<'EOF'
chr1 100 200 a1 1 + 1
chr1 300 400 a2 2 - 2
chr1 500 600 a3 3 . 0
chr1 700 800 a4 4 + 0
EOF
ck $EXP -S -c

# Zero-length A and zero-length B: points at distance <= 1 widen to
# overlapping intervals; distance 2 only touches.
A=$TMPD/f3e.bed
B=$TMPD/f3f.bed
rows > $A <<'EOF'
chr1 100 100 a1
chr1 102 102 a2
chr1 500 500 a3
EOF
rows > $B <<'EOF'
chr1 99 99 q1
chr1 100 100 q2
chr1 101 101 q3
chr1 102 102 q4
chr1 103 103 q5
chr1 104 104 q6
EOF
rows > $EXP <<'EOF'
chr1 100 100 a1 3
chr1 102 102 a2 3
chr1 500 500 a3 0
EOF
ck $EXP -c

# ===========================================================================
## fixture F4: hit order across UCSC bin levels (16 kb, 128 kb, 1 Mb, 8 Mb)
CLASS="hand F4: bin levels, coordinates across 16384, 131072, 1048576"
DERIV="hand (bin order)"
A=$TMPD/f4a.bed
B=$TMPD/f4b.bed
rows > $A <<'EOF'
chr1 0 2000000 big
chr1 16383 16385 x16k
chr1 131072 131073 x128k
chr1 1048576 1048577 x1m
EOF
# B bins: m1 level0 bin0; m2 [16383,16385) crosses 16384 -> level1 bin0;
# m3 level0 bin1; m4 level0 bin2; m5 [131071,131073) crosses 131072 ->
# level2 bin0; m6 level0 bin12; m7 [1048575,1048577) crosses 1048576 ->
# level3 bin0; m8 level0 bin64; m9 level0 bin115.
rows > $B <<'EOF'
chr1 100 200 m1
chr1 16383 16385 m2
chr1 16384 16400 m3
chr1 40000 40010 m4
chr1 131071 131073 m5
chr1 200000 200010 m6
chr1 1048575 1048577 m7
chr1 1048576 1048600 m8
chr1 1900000 1900010 m9
EOF
# This test pins the current output: the order of hits inside one A record
# in unsorted mode (smallest bin level first, ascending bin number, then B
# file order), which the help text does not specify. The order is level 0
# bins ascending (m1 m3 m4 m6 m8 m9), then level 1 m2, level 2 m5, level 3 m7.
rows > $EXP <<'EOF'
chr1 0 2000000 big chr1 100 200 m1 100
chr1 0 2000000 big chr1 16384 16400 m3 16
chr1 0 2000000 big chr1 40000 40010 m4 10
chr1 0 2000000 big chr1 200000 200010 m6 10
chr1 0 2000000 big chr1 1048576 1048600 m8 24
chr1 0 2000000 big chr1 1900000 1900010 m9 10
chr1 0 2000000 big chr1 16383 16385 m2 2
chr1 0 2000000 big chr1 131071 131073 m5 2
chr1 0 2000000 big chr1 1048575 1048577 m7 2
chr1 16383 16385 x16k chr1 16384 16400 m3 1
chr1 16383 16385 x16k chr1 16383 16385 m2 2
chr1 131072 131073 x128k chr1 131071 131073 m5 1
chr1 1048576 1048577 x1m chr1 1048576 1048600 m8 1
chr1 1048576 1048577 x1m chr1 1048575 1048577 m7 1
EOF
DERIV=pinned
cku $EXP -wo
DERIV=hand
# -sorted visits B in file order.
rows > $EXP <<'EOF'
chr1 0 2000000 big chr1 100 200 m1 100
chr1 0 2000000 big chr1 16383 16385 m2 2
chr1 0 2000000 big chr1 16384 16400 m3 16
chr1 0 2000000 big chr1 40000 40010 m4 10
chr1 0 2000000 big chr1 131071 131073 m5 2
chr1 0 2000000 big chr1 200000 200010 m6 10
chr1 0 2000000 big chr1 1048575 1048577 m7 2
chr1 0 2000000 big chr1 1048576 1048600 m8 24
chr1 0 2000000 big chr1 1900000 1900010 m9 10
chr1 16383 16385 x16k chr1 16383 16385 m2 2
chr1 16383 16385 x16k chr1 16384 16400 m3 1
chr1 131072 131073 x128k chr1 131071 131073 m5 1
chr1 1048576 1048577 x1m chr1 1048575 1048577 m7 1
chr1 1048576 1048577 x1m chr1 1048576 1048600 m8 1
EOF
cks $EXP -wo
# This test pins the current output: unsorted hit order for short A records
# whose hits sit in different bin levels (the smaller bin level comes first,
# so m3 before m2, and m8 before m7), which the help text does not specify.
A=$TMPD/f4e.bed
rows > $A <<'EOF2'
chr1 16383 16385 x16k
chr1 1048576 1048577 x1m
EOF2
rows > $EXP <<'EOF2'
chr1 16383 16385 x16k chr1 16384 16400 m3
chr1 16383 16385 x16k chr1 16383 16385 m2
chr1 1048576 1048577 x1m chr1 1048576 1048600 m8
chr1 1048576 1048577 x1m chr1 1048575 1048577 m7
EOF2
DERIV=pinned
cku $EXP -wa -wb
DERIV=hand
# This test pins the current output: same A records with -loj, unsorted
# order as above (smaller bin level first).
rows > $EXP <<'EOF2'
chr1 16383 16385 x16k chr1 16384 16400 m3
chr1 16383 16385 x16k chr1 16383 16385 m2
chr1 1048576 1048577 x1m chr1 1048576 1048600 m8
chr1 1048576 1048577 x1m chr1 1048575 1048577 m7
EOF2
DERIV=pinned
cku $EXP -loj
DERIV=hand
A=$TMPD/f4a.bed
# Same hits in either mode after sorting both outputs: -wa -wb.
rows > $EXP <<'EOF'
chr1 0 2000000 big chr1 100 200 m1
chr1 0 2000000 big chr1 16383 16385 m2
chr1 0 2000000 big chr1 16384 16400 m3
chr1 0 2000000 big chr1 40000 40010 m4
chr1 0 2000000 big chr1 131071 131073 m5
chr1 0 2000000 big chr1 200000 200010 m6
chr1 0 2000000 big chr1 1048575 1048577 m7
chr1 0 2000000 big chr1 1048576 1048600 m8
chr1 0 2000000 big chr1 1900000 1900010 m9
chr1 1048576 1048577 x1m chr1 1048575 1048577 m7
chr1 1048576 1048577 x1m chr1 1048576 1048600 m8
chr1 131072 131073 x128k chr1 131071 131073 m5
chr1 16383 16385 x16k chr1 16383 16385 m2
chr1 16383 16385 x16k chr1 16384 16400 m3
EOF
sort $EXP > $TMPD/f4sorted
for m in unsorted -sorted; do
	tbegin $m -wa -wb "(output sorted)"
	if [ $m = unsorted ]; then
		$BT intersect -wa -wb -a $A -b $B | sort > $OBS
	else
		$BT intersect -wa -wb -sorted -a $A -b $B | sort > $OBS
	fi
	check $OBS $TMPD/f4sorted
done
# -c and -v are order independent.
rows > $EXP <<'EOF'
chr1 0 2000000 big 9
chr1 16383 16385 x16k 2
chr1 131072 131073 x128k 1
chr1 1048576 1048577 x1m 2
EOF
ck $EXP -c
: > $EXP
ck $EXP -v
# Long B over many short A: each x record is a hit of the one long B.
A=$TMPD/f4c.bed
B=$TMPD/f4d.bed
rows > $A <<'EOF'
chr1 10 20 s1
chr1 16380 16390 s2
chr1 131070 131080 s3
chr1 1048570 1048580 s4
chr1 3000000 3000010 s5
EOF
rows > $B <<'EOF'
chr1 0 2000000 long
EOF
rows > $EXP <<'EOF'
chr1 10 20 s1 chr1 0 2000000 long 10
chr1 16380 16390 s2 chr1 0 2000000 long 10
chr1 131070 131080 s3 chr1 0 2000000 long 10
chr1 1048570 1048580 s4 chr1 0 2000000 long 10
chr1 3000000 3000010 s5 . -1 -1 . 0
EOF
ck $EXP -wao
# -F 0.00001 on the long B (hits cover 10/2000000 = 5e-6): none pass.
rows > $EXP <<'EOF'
chr1 10 20 s1 0
chr1 16380 16390 s2 0
chr1 131070 131080 s3 0
chr1 1048570 1048580 s4 0
chr1 3000000 3000010 s5 0
EOF
ck $EXP -c -F 0.00001

# ===========================================================================
## fixture F5: several chromosomes, lexicographic chromosome order
CLASS="hand F5: multiple chromosomes, chr1/chr10/chr2 ordering"
DERIV=hand
A=$TMPD/f5a.bed
B=$TMPD/f5b.bed
rows > $A <<'EOF'
chr1 100 200 a1
chr10 100 200 a2
chr2 100 200 a3
chr3 100 200 a4
EOF
rows > $B <<'EOF'
chr1 150 160 b1
chr10 190 300 b2
chr2 50 120 b3
chr2 150 160 b4
chrX 1 1000 b5
EOF
rows > $EXP <<'EOF'
chr1 150 160 a1
chr10 190 200 a2
chr2 100 120 a3
chr2 150 160 a3
EOF
ck $EXP
rows > $EXP <<'EOF'
chr1 100 200 a1 chr1 150 160 b1 10
chr10 100 200 a2 chr10 190 300 b2 10
chr2 100 200 a3 chr2 50 120 b3 20
chr2 100 200 a3 chr2 150 160 b4 10
chr3 100 200 a4 . -1 -1 . 0
EOF
ck $EXP -wao
rows > $EXP <<'EOF'
chr3 100 200 a4
EOF
ck $EXP -v
rows > $EXP <<'EOF'
chr1 100 200 a1 1
chr10 100 200 a2 1
chr2 100 200 a3 2
chr3 100 200 a4 0
EOF
ck $EXP -c
# -f 0.5 -u: every hit covers at most 20 of the 100 bases of its A record
# (a1 10, a2 10, a3 20 and 10), so nothing survives.
: > $EXP
ck $EXP -u -f 0.5
# -F 0.5 -u: b1 (10 of 10), b2 (10 of 110 no), b3 (20 of 70 no), b4 (10 of 10).
rows > $EXP <<'EOF'
chr1 100 200 a1
chr2 100 200 a3
EOF
ck $EXP -u -F 0.5

# ===========================================================================
## generated datasets compared with the awk reference
CLASS="reference awk: generated BED6, 4 A and 3 B chromosomes, 3 Mb, lengths up to 300 kb (bin levels 0-3)"
DERIV="reference awk"
awk -v seed=12345 -v n=240 -v nchr=4 -v maxc=3000000 -v maxlen=300000 -f $TMPD/gen.awk > $TMPD/wa.bed
awk -v seed=777 -v n=300 -v nchr=3 -v maxc=3000000 -v maxlen=300000 -f $TMPD/gen.awk > $TMPD/wb.bed
sort -k1,1 -k2,2n $TMPD/wa.bed > $TMPD/wa.sorted.bed
sort -k1,1 -k2,2n $TMPD/wb.bed > $TMPD/wb.sorted.bed
# Narrow set: every record inside one 16 kb bin (B file order in both modes).
awk -v seed=4242 -v n=250 -v nchr=3 -v maxc=12000 -v maxlen=400 -f $TMPD/gen.awk > $TMPD/na.bed
awk -v seed=99 -v n=300 -v nchr=2 -v maxc=12000 -v maxlen=400 -f $TMPD/gen.awk > $TMPD/nb.bed
sort -k1,1 -k2,2n $TMPD/na.bed > $TMPD/na.sorted.bed
sort -k1,1 -k2,2n $TMPD/nb.bed > $TMPD/nb.sorted.bed

# Sanity: the generated sets have the intended shape (hits exist, several
# levels occur); a degenerate generator would make every comparison vacuous.
tbegin unsorted generator-shape
nh=$(awk -v opts="-wo" -v order=file -f $TMPD/ref.awk $TMPD/wb.bed $TMPD/wa.bed | wc -l | tr -d ' ')
nl=$(awk -v shape=1 -f $TMPD/ref.awk $TMPD/wb.bed /dev/null)
if [ "$nh" -gt 300 ] && [ "$nl" -ge 3 ]; then echo ok; else FAILURES=$((FAILURES + 1)); echo "fail ($nh hits, $nl bin levels in B)"; fi

MODES="|-wa|-wb|-wa -wb|-wo|-wao|-loj|-u|-c|-v"
OLDIFS=$IFS
IFS='|'
for m in $MODES; do
	IFS=$OLDIFS
	# Wide set, unsorted input, unsorted mode: hit order by bin level.
	CLASS="reference awk: wide set (bin levels 0-3), unsorted input"
	refck $TMPD/wa.bed $TMPD/wb.bed unsorted $m
	# Wide set, sorted input, -sorted mode.
	CLASS="reference awk: wide set, sorted input"
	refck $TMPD/wa.sorted.bed $TMPD/wb.sorted.bed sorted $m
	# Narrow set (single bin), unsorted input, unsorted mode.
	CLASS="reference awk: narrow set (single 16 kb bin), unsorted input"
	refck $TMPD/na.bed $TMPD/nb.bed unsorted $m
	# Narrow set, sorted input, -sorted mode.
	CLASS="reference awk: narrow set, sorted input"
	refck $TMPD/na.sorted.bed $TMPD/nb.sorted.bed sorted $m
	IFS='|'
done
IFS=$OLDIFS

# Fractions, strands and their combinations on both sets and both modes.
# Fractions are far from rounding ties except where ov/len is an exact ratio
# of the same decimal (compared as real numbers in both models).
OPTS="-f 0.5|-F 0.5|-f 0.3 -F 0.6|-f 0.3 -F 0.6 -e|-f 0.5 -r|-f 0.2 -F 0.01 -e|-f 1.0|-F 1.0|-f 1 -F 0.2|-f 0.01|-F 0.001|-f 0.6 -r|-s|-S|-s -f 0.5|-S -F 0.25|-s -wo|-S -wao|-s -c|-S -c|-s -v|-S -v|-s -r -f 0.4|-f 0.5 -v|-f 0.5 -r -v|-f 0.5 -c|-F 0.5 -u|-f 0.5 -wao|-F 0.3 -loj|-f 0.5 -r -wo|-s -f 0.5 -wo|-S -f 0.5 -F 0.5 -e -wao|-s -F 0.25 -loj|-S -f 0.1 -u|-s -wa -wb|-f 0.25 -wb|-F 0.75 -wa"
IFS='|'
for o in $OPTS; do
	IFS=$OLDIFS
	CLASS="reference awk: wide set, unsorted input"
	refck $TMPD/wa.bed $TMPD/wb.bed unsorted $o
	CLASS="reference awk: wide set, sorted input"
	refck $TMPD/wa.sorted.bed $TMPD/wb.sorted.bed sorted $o
	IFS='|'
done
IFS=$OLDIFS

NOPTS="-f 0.5|-F 0.5 -wo|-f 0.3 -F 0.6 -e -wa|-f 0.75 -r -c|-f 1 -wao|-F 1 -v|-s -wo|-S -wao|-s -f 0.4 -c|-S -F 0.4 -u|-s -v|-S -v|-f 0.25 -F 0.25 -loj|-f 0.5 -r -wb"
IFS='|'
for o in $NOPTS; do
	IFS=$OLDIFS
	CLASS="reference awk: narrow set (single bin), unsorted input"
	refck $TMPD/na.bed $TMPD/nb.bed unsorted $o
	CLASS="reference awk: narrow set, sorted input"
	refck $TMPD/na.sorted.bed $TMPD/nb.sorted.bed sorted $o
	IFS='|'
done
IFS=$OLDIFS

# Exact-fraction grid: A is 100 long; B of length 100 shifted by multiples of
# 25 gives fractions exactly 0, .25, .5, .75, 1 of both A and B.
CLASS="reference awk: grid, exact fractions 0.25/0.5/0.75/1"
awk 'BEGIN { for (i = 0; i < 20; i++) printf "chr1\t%d\t%d\tA%d\t%d\t%s\n", (i + 1) * 1000, (i + 1) * 1000 + 100, i, i, (i % 2) ? "-" : "+" }' > $TMPD/ga.bed
awk 'BEGIN { for (i = 0; i < 20; i++) for (d = -100; d <= 100; d += 25) printf "chr1\t%d\t%d\tB%d_%d\t0\t%s\n", (i + 1) * 1000 + d, (i + 1) * 1000 + d + 100, i, d, ((i + d / 25 + 8) % 3 == 0) ? "+" : (((i + d / 25 + 8) % 3 == 1) ? "-" : ".") }' > $TMPD/gb.bed
# B starts d bases from the A start with the same length, so it overlaps A by
# 100-|d| bases (|d| = 100 only abuts).
# Sort B (starts are not monotone in the file) for -sorted mode.
sort -k1,1 -k2,2n $TMPD/gb.bed > $TMPD/gb.sorted.bed
GOPTS="-wo|-f 0.25 -wo|-f 0.5 -wo|-f 0.75 -wo|-f 1.0 -wo|-F 0.25 -wo|-F 0.5 -wo|-F 1 -wo|-f 0.25 -r -wo|-f 0.75 -r -wo|-f 0.5 -F 0.75 -wo|-f 0.5 -F 0.75 -e -wo|-f 0.75 -F 0.75 -e -c|-s -f 0.5 -wo|-S -f 0.5 -wo|-f 0.5 -v|-f 0.5 -u"
IFS='|'
for o in $GOPTS; do
	IFS=$OLDIFS
	refck $TMPD/ga.bed $TMPD/gb.sorted.bed unsorted $o
	refck $TMPD/ga.bed $TMPD/gb.sorted.bed sorted $o
	IFS='|'
done
IFS=$OLDIFS

# ===========================================================================
## option errors tied to the overlap rules
CLASS="hand: option validation"
DERIV=hand
A=$TMPD/f1a.bed
B=$TMPD/f1b.bed

# -e without -f is rejected.
tbegin unsorted -e -F 0.5
expect_error "-e must be used with -f" $BT intersect -e -F 0.5 -a $A -b $B
# -r with -F is rejected.
tbegin unsorted -f 0.5 -F 0.5 -r
expect_error "-r must be used solely with -f" $BT intersect -f 0.5 -F 0.5 -r -a $A -b $B
# -f above 1 is rejected.
tbegin unsorted -f 1.5
expect_error "-f must be in the range" $BT intersect -f 1.5 -a $A -b $B
# -F of 0 is rejected.
tbegin unsorted -F 0
expect_error "-F must be in the range" $BT intersect -F 0 -a $A -b $B
# -u with -v is rejected.
tbegin unsorted -u -v
expect_error "-u for anyHit OR -v for noHit" $BT intersect -u -v -a $A -b $B
# -wo with -wa is rejected.
tbegin unsorted -wo -wa
expect_error "-wa for writeA OR -wo" $BT intersect -wo -wa -a $A -b $B
# -s with -S is rejected.
tbegin unsorted -s -S
expect_error "-s for sameStrand, or -S for diffStrand" $BT intersect -s -S -a $A -b $B

finish
