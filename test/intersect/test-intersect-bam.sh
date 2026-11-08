. ./helpers-intersect.sh bam
set -e

# BAM/CRAM input and output and -split (blocked intervals) for
# `bedtools intersect`. Every BAM is built from SAM text written below, so
# every expected value is worked out from the coordinates in that text.
# BAM -> 0-based half-open: start = POS-1, end = start + reference span of the
# CIGAR (M, =, X, D, N; not I, S, H, P).

# Never let htslib look for references on the network or in a host cache.
unset CRAM_REFERENCE
export REF_PATH="$TMPD/no-ref-path" REF_CACHE="$TMPD/no-ref-cache"

# ---------------------------------------------------------------------------
# Alignments (SAM fields; SEQ and QUAL are *). Block structure per read:
#   fwd  [100,150)                  one block, +
#   rev  [200,240)                  one block, - (flag 0x10)
#   pe1  [300,330)  flag 99 = first in pair, name prints as pe1/1 with -bed
#   pe2  [400,430)  flag 147 = second in pair, - strand, pe2/2 with -bed
#   spl  10M100N10M [500,620)       blocks [500,510) [610,620)
#   del  10M20D10M  [700,740)       one block [700,740) (D stays in a block)
#   cli  2S10M3I10M2H [800,820)     one block [800,820), S/I/H not counted
#   c2   chr2 [10,30)
#   dn   5M10D5M20N5M chr2 [100,145) blocks [100,120) [140,145)
#   xe   5=2X5=  chr2 [300,312)     one block
#   plc  unmapped (flag 4) but placed at chr1:111
#   unm  unmapped, no chromosome
# ---------------------------------------------------------------------------
FWD="fwd 0 chr1 101 30 50M * 0 0 * *"
REV="rev 16 chr1 201 20 40M * 0 0 * *"
PE1="pe1 99 chr1 301 40 30M = 401 130 * *"
PE2="pe2 147 chr1 401 41 30M = 301 -130 * *"
SPL="spl 0 chr1 501 50 10M100N10M * 0 0 * *"
DEL="del 0 chr1 701 25 10M20D10M * 0 0 * *"
CLI="cli 0 chr1 801 15 2S10M3I10M2H * 0 0 * *"
C2="c2 0 chr2 11 5 20M * 0 0 * *"
DN="dn 0 chr2 101 10 5M10D5M20N5M * 0 0 * *"
XE="xe 0 chr2 301 10 5=2X5= * 0 0 * *"
PLC="plc 4 chr1 111 0 * * 0 0 * *"
UNM="unm 4 * 0 0 * * 0 0 * *"

# BED12 text of the same alignments as printed by -bed.
E_FWD="chr1 100 150 fwd 30 + 100 150 0,0,0 1 50, 0,"
E_REV="chr1 200 240 rev 20 - 200 240 0,0,0 1 40, 0,"
E_PE1="chr1 300 330 pe1/1 40 + 300 330 0,0,0 1 30, 0,"
E_PE2="chr1 400 430 pe2/2 41 - 400 430 0,0,0 1 30, 0,"
E_SPL="chr1 500 620 spl 50 + 500 620 0,0,0 2 10,10, 0,110,"
E_DEL="chr1 700 740 del 25 + 700 740 0,0,0 1 40, 0,"
E_CLI="chr1 800 820 cli 15 + 800 820 0,0,0 1 20, 0,"
E_C2="chr2 10 30 c2 5 + 10 30 0,0,0 1 20, 0,"
E_DN="chr2 100 145 dn 10 + 100 145 0,0,0 2 20,5, 0,40,"
E_XE="chr2 300 312 xe 10 + 300 312 0,0,0 1 12, 0,"

HDR1="@HD VN:1.6 SO:coordinate"
HDR2="@SQ SN:chr1 LN:1000"
HDR3="@SQ SN:chr2 LN:1000"

# mkbam <name> <records...>: write $TMPD/<name>.sam, make BAM.
mkbam()
{
	local name=$1
	shift
	{
		printf '%s\n' "$HDR1" "$HDR2" "$HDR3"
		printf '%s\n' "$@"
	} | rows > "$TMPD/$name.sam"
	$htsutil samtobam "$TMPD/$name.sam" "$TMPD/$name.bam"
}

# all.bam: every read incl. unmapped; m.bam: mapped reads only;
# m1.bam: mapped reads of chr1 only.
mkbam all "$FWD" "$PLC" "$REV" "$PE1" "$PE2" "$SPL" "$DEL" "$CLI" "$C2" "$DN" "$XE" "$UNM"
mkbam m "$FWD" "$REV" "$PE1" "$PE2" "$SPL" "$DEL" "$CLI" "$C2" "$DN" "$XE"
mkbam m1 "$FWD" "$REV" "$PE1" "$PE2" "$SPL" "$DEL" "$CLI"
A=$TMPD/all.bam
M=$TMPD/m.bam
M1=$TMPD/m1.bam

# B intervals. Overlap with the alignment spans (no -split) / with blocks (-split):
#   chr1 120-130  fwd 10/10        chr1 230-250  rev 10/10
#   chr1 405-415  pe2 10/10        chr1 550-560  spl 10 / none (N gap)
#   chr1 605-615  spl 10 / 5       chr1 715-725  del 10/10 (inside the D)
#   chr2 108-112  dn 4/4 (D)       chr2 125-135  dn 10 / none (N gap)
#   chr2 305-308  xe 3/3
printf 'chr1\t120\t130\nchr1\t230\t250\nchr1\t405\t415\nchr1\t550\t560\nchr1\t605\t615\nchr1\t715\t725\nchr2\t108\t112\nchr2\t125\t135\nchr2\t305\t308\n' > "$TMPD/b.bed"
B=$TMPD/b.bed
B1=$TMPD/b1.bed
printf 'chr1\t120\t130\nchr1\t230\t250\nchr1\t405\t415\nchr1\t550\t560\nchr1\t605\t615\nchr1\t715\t725\n' > "$TMPD/b1.bed"
printf 'chr1\t0\t1000\nchr2\t0\t1000\n' > "$TMPD/big.bed"
# Stranded B: x + and y - over chr1 100-300, z - over chr1 300-450, w + over chr2 0-100.
printf 'chr1\t100\t300\tx\t0\t+\nchr1\t100\t300\ty\t0\t-\nchr1\t300\t450\tz\t0\t-\nchr2\t0\t100\tw\t0\t+\n' > "$TMPD/s.bed"

# check_warning <stderr file> <text>: obs/exp comparison that also requires the
# warning text on stderr; a missing warning counts as a failure.
check_warning()
{
	if grep -F -e "$2" "$1" > /dev/null; then
		check "$TMPD/obs" "$TMPD/exp"
	else
		FAILURES=$((FAILURES + 1))
		echo fail
		echo "      missing warning: $2"
	fi
}

# run_bam <intersect args>: BAM output as SAM records (no header) into obs.
run_bam()
{
	$BT intersect "$@" 2> /dev/null | $htsutil viewbamrecords > "$TMPD/obs"
}
# run_txt <intersect args>: text output into obs.
run_txt()
{
	$BT intersect "$@" 2> /dev/null > "$TMPD/obs"
}
# lines <text...>: expected file of tab-joined lines (spaces become tabs).
lines()
{
	printf '%s\n' "$@" | rows > "$TMPD/exp"
}
# none: expected file is empty.
none()
{
	: > "$TMPD/exp"
}
# checkgrouped: compare obs and exp in output order, except that lines of the
# same A record (same field 4) are sorted among themselves, because the order
# of several hits of one alignment is not specified. Expected lines are written
# in A order.
checkgrouped()
{
	for f in obs exp; do
		awk -F"\t" 'function flush(  i, j, t) { for (i = 2; i <= n; i++) { t = g[i]; for (j = i - 1; j >= 1 && g[j] > t; j--) g[j+1] = g[j]; g[j+1] = t } for (i = 1; i <= n; i++) print g[i]; n = 0 } { if (n > 0 && $4 != k) flush(); k = $4; g[++n] = $0 } END { flush() }' "$TMPD/$f" > "$TMPD/$f.g"
	done
	check "$TMPD/obs.g" "$TMPD/exp.g"
}

###########################################################################
# BAM as A, BAM output
###########################################################################

# BAM -a, BED -b: output is BAM holding each overlapping alignment once, in A order; unmapped reads are never hits.
# targets: -a bam -b bed | all.bam b.bed | hand
begin 001
run_bam -a $A -b $B
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# Obsolete -abam spelling gives the same BAM output as -a.
# targets: -abam | all.bam b.bed | hand
begin 002
run_bam -abam $A -b $B
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -abam with no file following it is an error.
# targets: -abam | missing argument | hand
begin 003
expect_error "-abam option given, but no query BAM file specified" $BT intersect -abam

# -wa with BAM A: the record is written once even when it has two hits (spl).
# targets: -wa | all.bam b.bed | hand
begin 004
run_bam -a $A -b $B -wa
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -u with BAM A: each overlapping alignment once.
# targets: -u | all.bam b.bed | hand
begin 005
run_bam -a $A -b $B -u
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -v with BAM A: reads without any hit, including the placed and the unplaced unmapped reads, in A order.
# targets: -v | all.bam b.bed | hand
begin 006
run_bam -a $A -b $B -v
lines "$PLC" "$PE1" "$CLI" "$C2" "$UNM"
check "$TMPD/obs" "$TMPD/exp"

# -v with a B that covers every mapped read leaves exactly the unmapped reads.
# targets: -v unmapped | all.bam big.bed | hand
begin 007
run_bam -a $A -b $TMPD/big.bed -v
lines "$PLC" "$UNM"
check "$TMPD/obs" "$TMPD/exp"

# Without -v the unmapped reads (placed or not) are never reported as hits.
# targets: unmapped in A | all.bam big.bed | hand
begin 008
run_bam -a $A -b $TMPD/big.bed
lines "$FWD" "$REV" "$PE1" "$PE2" "$SPL" "$DEL" "$CLI" "$C2" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -c with BAM output is an error unless -bed is given.
# targets: -c | bam A without -bed | hand
begin 009
expect_error "writeCount option is not valid with BAM query input" $BT intersect -a $A -b $B -c

# -wo with BAM output is an error unless -bed is given.
# targets: -wo | bam A without -bed | hand
begin 010
expect_error "writeAllOverlap option is not valid with BAM query input" $BT intersect -a $A -b $B -wo

# -wao with BAM output is an error unless -bed is given.
# targets: -wao | bam A without -bed | hand
begin 011
expect_error "writeAllOverlap option is not valid with BAM query input" $BT intersect -a $A -b $B -wao

# -wb with BAM output warns that it is ignored and writes the same records as without it.
# targets: -wb | all.bam b.bed | hand
begin 012
expect_status 0 $BT intersect -a $A -b $B -wb
$htsutil viewbamrecords "$TMPD/st.stdout" > "$TMPD/obs"
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check_warning "$TMPD/st.stderr" "-wb and -loj are ignored with bam input"

# -header with BAM input warns that it is ignored; the BAM header is kept unchanged and the records are as without it.
# targets: -header | all.bam b.bed | hand
begin 013
$BT intersect -a $A -b $B -header 2> "$TMPD/err" | $htsutil viewbam > "$TMPD/obs"
lines "$HDR1" "$HDR2" "$HDR3" "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check_warning "$TMPD/err" "-header option is not valid for BAM input"

# -header -bed with BAM input warns and writes BED12 with no header lines.
# targets: -header -bed | m.bam b.bed | hand
begin 014
$BT intersect -a $M -b $B1 -header -bed -u 2> "$TMPD/err" > "$TMPD/obs"
lines "$E_FWD" "$E_REV" "$E_PE2" "$E_SPL" "$E_DEL"
check_warning "$TMPD/err" "-header option is not valid for BAM input"

# The BAM header text and reference list are copied to the output.
# targets: header copy | all.bam b.bed | hand
begin 015
$BT intersect -a $A -b $B 2> /dev/null | $htsutil viewbam | grep '^@' > "$TMPD/obs"
lines "$HDR1" "$HDR2" "$HDR3"
check "$TMPD/obs" "$TMPD/exp"

# Output records are the input records unchanged: flags, mate fields and TLEN are kept and no /1 /2 suffix is added to names in BAM output.
# targets: BAM output fidelity | pe1 pe2 | hand
begin 016
printf 'chr1\t290\t440\n' > "$TMPD/pe.bed"
run_bam -a $M -b "$TMPD/pe.bed"
lines "$PE1" "$PE2"
check "$TMPD/obs" "$TMPD/exp"

# -f 0.25 with BAM output (no -split): overlap / alignment span >= 0.25. rev 10/40, pe2 10/30, del 10/40, xe 3/12 pass; fwd 10/50, spl 10/120, dn 4/45 and 10/45 fail.
# targets: -f | m.bam b.bed | hand
begin 017
run_bam -a $M -b $B -f 0.25
lines "$REV" "$PE2" "$DEL" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -F 0.6 with BAM output: overlap / B length >= 0.6. Only chr1 230-250 (10/20 with rev) fails.
# targets: -F | m.bam b.bed | hand
begin 018
run_bam -a $M -b $B -F 0.6
lines "$FWD" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -f 0.5 -r with BAM output: B must cover half of the alignment and the alignment half of B; no read is half covered by any 3-20 bp B interval, so the output is empty.
# targets: -f -r | m.bam b.bed | hand
begin 019
run_bam -a $M -b $B -f 0.5 -r
none
check "$TMPD/obs" "$TMPD/exp"

# -f 0.5 -F 0.5 -e: either fraction suffices. B chr1 230-250 vs rev: 10/20 of B passes; every other B interval is fully inside its read, so -F is 1.0 and all hit reads pass.
# targets: -f -F -e | m.bam b.bed | hand
begin 020
run_bam -a $M -b $B -f 0.5 -F 0.5 -e
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -s with BAM A: strand comes from flag 0x10. fwd(+) x(+); rev(-) y(-); pe2(-) z(-); c2(+) w(+). pe1 is + but z is -.
# targets: -s | all.bam s.bed | hand
begin 021
run_bam -a $A -b $TMPD/s.bed -s
lines "$FWD" "$REV" "$PE2" "$C2"
check "$TMPD/obs" "$TMPD/exp"

# -S with BAM A: opposite strands: fwd(+) y(-); rev(-) x(+); pe1(+) z(-).
# targets: -S | all.bam s.bed | hand
begin 022
run_bam -a $A -b $TMPD/s.bed -S
lines "$FWD" "$REV" "$PE1"
check "$TMPD/obs" "$TMPD/exp"

# Strand from the flag for paired and other flag combinations: 0x10 (flags 16, 83) is -, others (0, 1, 33, 65, 129, 193) are +.
# targets: -s -S | pairs.bam | hand
begin 023
P0="p0 1 chr1 11 7 10M * 0 0 * *"
P1="p1 65 chr1 11 8 10M * 0 0 * *"
P2="p2 129 chr1 11 9 10M * 0 0 * *"
P3="p3 193 chr1 11 10 10M * 0 0 * *"
P4="p4 0 chr1 11 255 10M * 0 0 * *"
P5="p5 16 chr1 11 0 10M * 0 0 * *"
P6="p6 83 chr1 11 11 10M * 0 0 * *"
P7="p7 33 chr1 11 12 10M * 0 0 * *"
mkbam pairs "$P0" "$P1" "$P2" "$P3" "$P4" "$P5" "$P6" "$P7"
printf 'chr1\t0\t100\tplus\t0\t+\n' > "$TMPD/plus.bed"
printf 'chr1\t0\t100\tminus\t0\t-\n' > "$TMPD/minus.bed"
run_bam -a "$TMPD/pairs.bam" -b "$TMPD/plus.bed" -s
lines "$P0" "$P1" "$P2" "$P3" "$P4" "$P7"
check "$TMPD/obs" "$TMPD/exp"

# -S against a + B keeps only the reverse-strand reads (flags 16 and 83).
# targets: -S | pairs.bam | hand
begin 024
run_bam -a "$TMPD/pairs.bam" -b "$TMPD/plus.bed" -S
lines "$P5" "$P6"
check "$TMPD/obs" "$TMPD/exp"

# -s against a - B keeps only the reverse-strand reads.
# targets: -s | pairs.bam | hand
begin 025
run_bam -a "$TMPD/pairs.bam" -b "$TMPD/minus.bed" -s
lines "$P5" "$P6"
check "$TMPD/obs" "$TMPD/exp"

# -ubam: the output is a valid BAM whose records equal the default output (the option is documented as "uncompressed"; compression is not checked here).
# targets: -ubam | all.bam b.bed | hand
begin 026
run_bam -a $A -b $B -ubam
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -ubam with -v keeps the same records as -v.
# targets: -ubam -v | all.bam b.bed | hand
begin 027
run_bam -a $A -b $B -ubam -v
lines "$PLC" "$PE1" "$CLI" "$C2" "$UNM"
check "$TMPD/obs" "$TMPD/exp"

# -fbam is accepted and does not change which records are written.
# targets: -fbam | all.bam b.bed | hand
begin 028
run_bam -a $A -b $B -fbam
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -fbam with -bed gives the same BED12 as without it.
# targets: -fbam -bed | m.bam b.bed | hand
begin 029
run_txt -a $M -b $B -fbam -bed -u
lines "$E_FWD" "$E_REV" "$E_PE2" "$E_SPL" "$E_DEL" "$E_DN" "$E_XE"
check "$TMPD/obs" "$TMPD/exp"

# -sorted with a coordinate-sorted BAM A gives the same records as the default algorithm.
# targets: -sorted | all.bam b.bed | hand
begin 030
run_bam -a $A -b $B -sorted
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -sorted -v with BAM A, including unmapped reads at the end of the file.
# targets: -sorted -v | all.bam b.bed | hand
begin 031
run_bam -a $A -b $B -sorted -v
lines "$PLC" "$PE1" "$CLI" "$C2" "$UNM"
check "$TMPD/obs" "$TMPD/exp"

# -sorted -s with BAM A.
# targets: -sorted -s | all.bam s.bed | hand
begin 032
run_bam -a $A -b $TMPD/s.bed -sorted -s
lines "$FWD" "$REV" "$PE2" "$C2"
check "$TMPD/obs" "$TMPD/exp"

# BAM A read from standard input.
# targets: -a stdin | all.bam b.bed | hand
begin 033
cat $A | $BT intersect -a stdin -b $B 2> /dev/null | $htsutil viewbamrecords > "$TMPD/obs"
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

###########################################################################
# BAM as A, -bed output
###########################################################################

# -bed alone: each hit prints A trimmed to the overlap (BED12 line keeps the thick and block fields of the whole alignment); name gets /2 for pe2; score is MAPQ; strand from flag. Two hits (spl, dn) give two lines.
# targets: -bed | m.bam b.bed | hand
begin 034
run_txt -a $M -b $B -bed
lines "chr1 120 130 fwd 30 + 100 150 0,0,0 1 50, 0," \
	"chr1 230 240 rev 20 - 200 240 0,0,0 1 40, 0," \
	"chr1 405 415 pe2/2 41 - 400 430 0,0,0 1 30, 0," \
	"chr1 550 560 spl 50 + 500 620 0,0,0 2 10,10, 0,110," \
	"chr1 605 615 spl 50 + 500 620 0,0,0 2 10,10, 0,110," \
	"chr1 715 725 del 25 + 700 740 0,0,0 1 40, 0," \
	"chr2 108 112 dn 10 + 100 145 0,0,0 2 20,5, 0,40," \
	"chr2 125 135 dn 10 + 100 145 0,0,0 2 20,5, 0,40," \
	"chr2 305 308 xe 10 + 300 312 0,0,0 1 12, 0,"
checkgrouped

# -bed -wa: the whole alignment per hit, so spl and dn print twice.
# targets: -bed -wa | m.bam b.bed | hand
begin 035
run_txt -a $M -b $B -bed -wa
lines "$E_FWD" "$E_REV" "$E_PE2" "$E_SPL" "$E_SPL" "$E_DEL" "$E_DN" "$E_DN" "$E_XE"
checkgrouped

# -bed -wa -wb: BED12 of the alignment followed by the B record.
# targets: -bed -wa -wb | m.bam b.bed | hand
begin 036
run_txt -a $M -b $B -bed -wa -wb
lines "$E_FWD chr1 120 130" "$E_REV chr1 230 250" "$E_PE2 chr1 405 415" \
	"$E_SPL chr1 550 560" "$E_SPL chr1 605 615" "$E_DEL chr1 715 725" \
	"$E_DN chr2 108 112" "$E_DN chr2 125 135" "$E_XE chr2 305 308"
checkgrouped

# -bed -wo: BED12, B record and the overlap in bases.
# targets: -bed -wo | m.bam b.bed | hand
begin 037
run_txt -a $M -b $B -bed -wo
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 550 560 10" "$E_SPL chr1 605 615 10" "$E_DEL chr1 715 725 10" \
	"$E_DN chr2 108 112 4" "$E_DN chr2 125 135 10" "$E_XE chr2 305 308 3"
checkgrouped

# -bed -wao: as -wo plus unhit alignments (pe1, cli, c2) with a null B and 0 overlap.
# targets: -bed -wao | m.bam b.bed | hand
begin 038
run_txt -a $M -b $B -bed -wao
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE1 . -1 -1 0" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 550 560 10" "$E_SPL chr1 605 615 10" "$E_DEL chr1 715 725 10" \
	"$E_CLI . -1 -1 0" "$E_C2 . -1 -1 0" \
	"$E_DN chr2 108 112 4" "$E_DN chr2 125 135 10" "$E_XE chr2 305 308 3"
checkgrouped

# -bed -loj: every alignment with each hit, or a null B (. -1 -1) when there is none.
# targets: -bed -loj | m.bam b.bed | hand
begin 039
run_txt -a $M -b $B -bed -loj
lines "$E_FWD chr1 120 130" "$E_REV chr1 230 250" "$E_PE1 . -1 -1" "$E_PE2 chr1 405 415" \
	"$E_SPL chr1 550 560" "$E_SPL chr1 605 615" "$E_DEL chr1 715 725" \
	"$E_CLI . -1 -1" "$E_C2 . -1 -1" \
	"$E_DN chr2 108 112" "$E_DN chr2 125 135" "$E_XE chr2 305 308"
checkgrouped

# -bed -c: number of B hits per alignment, 0 for none. spl and dn have 2.
# targets: -bed -c | m.bam b.bed | hand
begin 040
run_txt -a $M -b $B -bed -c
lines "$E_FWD 1" "$E_REV 1" "$E_PE1 0" "$E_PE2 1" "$E_SPL 2" "$E_DEL 1" "$E_CLI 0" "$E_C2 0" "$E_DN 2" "$E_XE 1"
check "$TMPD/obs" "$TMPD/exp"

# -bed -u: each hit alignment once as BED12.
# targets: -bed -u | m.bam b.bed | hand
begin 041
run_txt -a $M -b $B -bed -u
lines "$E_FWD" "$E_REV" "$E_PE2" "$E_SPL" "$E_DEL" "$E_DN" "$E_XE"
check "$TMPD/obs" "$TMPD/exp"

# -bed -v: alignments without any hit.
# targets: -bed -v | m.bam b.bed | hand
begin 042
run_txt -a $M -b $B -bed -v
lines "$E_PE1" "$E_CLI" "$E_C2"
check "$TMPD/obs" "$TMPD/exp"

# -bed -wo with -f 0.25: only hits whose overlap is >= 25% of the alignment span (rev 10/40, pe2 10/30, del 10/40, xe 3/12).
# targets: -bed -wo -f | m.bam b.bed | hand
begin 043
run_txt -a $M -b $B -bed -wo -f 0.25
lines "$E_REV chr1 230 250 10" "$E_PE2 chr1 405 415 10" "$E_DEL chr1 715 725 10" "$E_XE chr2 305 308 3"
checkgrouped

# -bed -wo -s: only the B on the same strand: fwd(+)/x, rev(-)/y, pe2(-)/z, c2(+)/w.
# targets: -bed -wo -s | m.bam s.bed | hand
begin 044
run_txt -a $M -b $TMPD/s.bed -bed -wo -s
lines "$E_FWD chr1 100 300 x 0 + 50" "$E_REV chr1 100 300 y 0 - 40" "$E_PE2 chr1 300 450 z 0 - 30" "$E_C2 chr2 0 100 w 0 + 20"
checkgrouped

# -bed -wo -S: opposite strands: fwd(+)/y, rev(-)/x, pe1(+)/z.
# targets: -bed -wo -S | m.bam s.bed | hand
begin 045
run_txt -a $M -b $TMPD/s.bed -bed -wo -S
lines "$E_FWD chr1 100 300 y 0 - 50" "$E_REV chr1 100 300 x 0 + 40" "$E_PE1 chr1 300 450 z 0 - 30"
checkgrouped

# -bed names: flag 0x40 gives /1, 0x80 gives /2, both give /1 (first checked first), neither gives no suffix; score is MAPQ; strand follows 0x10.
# targets: -bed names/score/strand | pairs.bam big.bed | hand
begin 046
run_txt -a "$TMPD/pairs.bam" -b "$TMPD/big.bed" -bed -wa
lines "chr1 10 20 p0 7 + 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p1/1 8 + 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p2/2 9 + 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p3/1 10 + 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p4 255 + 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p5 0 - 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p6/1 11 - 10 20 0,0,0 1 10, 0," \
	"chr1 10 20 p7 12 + 10 20 0,0,0 1 10, 0,"
check "$TMPD/obs" "$TMPD/exp"

# -bed with -wb and no -wa prints A trimmed to the overlap, then B (spl, N gap hit included without -split).
# targets: -bed -wb | m.bam b.bed | hand
begin 047
run_txt -a $M -b $B1 -bed -wb
lines "chr1 120 130 fwd 30 + 100 150 0,0,0 1 50, 0, chr1 120 130" \
	"chr1 230 240 rev 20 - 200 240 0,0,0 1 40, 0, chr1 230 250" \
	"chr1 405 415 pe2/2 41 - 400 430 0,0,0 1 30, 0, chr1 405 415" \
	"chr1 550 560 spl 50 + 500 620 0,0,0 2 10,10, 0,110, chr1 550 560" \
	"chr1 605 615 spl 50 + 500 620 0,0,0 2 10,10, 0,110, chr1 605 615" \
	"chr1 715 725 del 25 + 700 740 0,0,0 1 40, 0, chr1 715 725"
checkgrouped

# -bed -sorted -wo gives the same lines as the default algorithm.
# targets: -bed -sorted -wo | m.bam b.bed | hand
begin 048
run_txt -a $M -b $B -bed -sorted -wo
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 550 560 10" "$E_SPL chr1 605 615 10" "$E_DEL chr1 715 725 10" \
	"$E_DN chr2 108 112 4" "$E_DN chr2 125 135 10" "$E_XE chr2 305 308 3"
checkgrouped

# CIGAR spans: P is ignored and the pad read forms one block; S, I and H do not extend the span (cli [800,820)); N splits into three blocks for tri.
# targets: -bed cigar P N | cigar.bam big.bed | hand
begin 049
TP="tp 0 chr1 901 9 5M2P5M * 0 0 * *"
TRI="tri 0 chr1 951 9 5M10N5M10N5M * 0 0 * *"
mkbam cigar "$CLI" "$TP" "$TRI"
run_txt -a "$TMPD/cigar.bam" -b "$TMPD/big.bed" -bed -wa
lines "$E_CLI" "chr1 900 910 tp 9 + 900 910 0,0,0 1 10, 0," "chr1 950 985 tri 9 + 950 985 0,0,0 3 5,5,5, 0,15,30,"
check "$TMPD/obs" "$TMPD/exp"

###########################################################################
# BAM as B
###########################################################################

# BED A with BAM B, -wa -wb: the BAM record is printed as BED6 (name with /2, MAPQ, strand).
# targets: -wa -wb bam B | b.bed m.bam | hand
begin 050
run_txt -a $B -b $M -wa -wb
lines "chr1 120 130 chr1 100 150 fwd 30 +" "chr1 230 250 chr1 200 240 rev 20 -" \
	"chr1 405 415 chr1 400 430 pe2/2 41 -" "chr1 550 560 chr1 500 620 spl 50 +" \
	"chr1 605 615 chr1 500 620 spl 50 +" "chr1 715 725 chr1 700 740 del 25 +" \
	"chr2 108 112 chr2 100 145 dn 10 +" "chr2 125 135 chr2 100 145 dn 10 +" \
	"chr2 305 308 chr2 300 312 xe 10 +"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, -wo: BED6 of the alignment span and the overlap with the whole span (N gap counts without -split).
# targets: -wo bam B | b.bed m.bam | hand
begin 051
run_txt -a $B -b $M -wo
lines "chr1 120 130 chr1 100 150 fwd 30 + 10" "chr1 230 250 chr1 200 240 rev 20 - 10" \
	"chr1 405 415 chr1 400 430 pe2/2 41 - 10" "chr1 550 560 chr1 500 620 spl 50 + 10" \
	"chr1 605 615 chr1 500 620 spl 50 + 10" "chr1 715 725 chr1 700 740 del 25 + 10" \
	"chr2 108 112 chr2 100 145 dn 10 + 4" "chr2 125 135 chr2 100 145 dn 10 + 10" \
	"chr2 305 308 chr2 300 312 xe 10 + 3"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, default output: A trimmed to the overlap with the alignment span.
# targets: default bam B | b.bed m.bam | hand
begin 052
run_txt -a "$TMPD/b1.bed" -b $M
lines "chr1 120 130" "chr1 230 240" "chr1 405 415" "chr1 550 560" "chr1 605 615" "chr1 715 725"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, -c: number of BAM records overlapping each A interval. A spanning fwd and rev counts 2.
# targets: -c bam B | bed m.bam | hand
begin 053
printf 'chr1\t110\t210\tq1\nchr1\t1\t20\tq2\nchr1\t0\t1000\tq3\n' > "$TMPD/q.bed"
run_txt -a "$TMPD/q.bed" -b $M -c
lines "chr1 110 210 q1 2" "chr1 1 20 q2 0" "chr1 0 1000 q3 7"
check "$TMPD/obs" "$TMPD/exp"

# Unmapped reads in B are ignored: a chromosome-wide A interval counts the 7 mapped chr1 reads (not plc) and 3 mapped chr2 reads (not unm).
# targets: -c unmapped in B | big.bed all.bam | hand
begin 054
run_txt -a "$TMPD/big.bed" -b $A -c
lines "chr1 0 1000 7" "chr2 0 1000 3"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, -v: A intervals not overlapping any mapped read.
# targets: -v bam B | bed all.bam | hand
begin 055
printf 'chr1\t0\t50\nchr1\t110\t112\nchr1\t160\t190\nchr1\t900\t950\n' > "$TMPD/v.bed"
run_txt -a "$TMPD/v.bed" -b $A -v
lines "chr1 0 50" "chr1 160 190" "chr1 900 950"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, -loj where every A interval has a hit: A then the BED6 of the alignment.
# targets: -loj bam B | bed m.bam | hand
begin 056
printf 'chr1\t110\t120\tq1\nchr1\t410\t420\tq2\n' > "$TMPD/lj.bed"
run_txt -a "$TMPD/lj.bed" -b $M -loj
lines "chr1 110 120 q1 chr1 100 150 fwd 30 +" "chr1 410 420 q2 chr1 400 430 pe2/2 41 -"
check "$TMPD/obs" "$TMPD/exp"

# BED A with BAM B, -wao where every A interval has a hit: adds the overlap.
# targets: -wao bam B | bed m.bam | hand
begin 057
run_txt -a "$TMPD/lj.bed" -b $M -wao
lines "chr1 110 120 q1 chr1 100 150 fwd 30 + 10" "chr1 410 420 q2 chr1 400 430 pe2/2 41 - 10"
check "$TMPD/obs" "$TMPD/exp"

# Stranded BED A with BAM B, -s -wa -wb: BAM strand from flag 0x10.
# targets: -s bam B | s.bed m.bam | hand
begin 058
run_txt -a "$TMPD/s.bed" -b $M -wa -wb -s
lines "chr1 100 300 x 0 + chr1 100 150 fwd 30 +" "chr1 100 300 y 0 - chr1 200 240 rev 20 -" \
	"chr1 300 450 z 0 - chr1 400 430 pe2/2 41 -" "chr2 0 100 w 0 + chr2 10 30 c2 5 +"
check "$TMPD/obs" "$TMPD/exp"

# Same with -S.
# targets: -S bam B | s.bed m.bam | hand
begin 059
run_txt -a "$TMPD/s.bed" -b $M -wa -wb -S
lines "chr1 100 300 x 0 + chr1 200 240 rev 20 -" "chr1 100 300 y 0 - chr1 100 150 fwd 30 +" \
	"chr1 300 450 z 0 - chr1 300 330 pe1/1 40 +"
checkgrouped

# BED A with BAM B and -sorted, -wa -wb (B restricted to chr1 so both files have the same chromosomes).
# targets: -sorted bam B | b1.bed m1.bam | hand
begin 060
run_txt -a "$TMPD/b1.bed" -b $M1 -wa -wb -sorted
lines "chr1 120 130 chr1 100 150 fwd 30 +" "chr1 230 250 chr1 200 240 rev 20 -" \
	"chr1 405 415 chr1 400 430 pe2/2 41 -" "chr1 550 560 chr1 500 620 spl 50 +" \
	"chr1 605 615 chr1 500 620 spl 50 +" "chr1 715 725 chr1 700 740 del 25 +"
check "$TMPD/obs" "$TMPD/exp"

# Several -b files mixing BAM and BED with -names, -wa -wb -wo: output carries the B name and BED6 for the BAM; q1 hits fwd (20 bp) and c1 (15 bp), q2 hits pe1 (15 bp) only.
# targets: multiple -b -names -wo | bed + bam + bed | hand
begin 061
printf 'chr1\t110\t130\tq1\nchr1\t305\t320\tq2\n' > "$TMPD/mq.bed"
printf 'chr1\t100\t125\tc1\n' > "$TMPD/mc.bed"
run_txt -a "$TMPD/mq.bed" -b $M "$TMPD/mc.bed" -names bam bed -wo
lines "chr1 110 130 q1 bam chr1 100 150 fwd 30 + 20" "chr1 110 130 q1 bed chr1 100 125 c1 15" \
	"chr1 305 320 q2 bam chr1 300 330 pe1/1 40 + 15"
checkgrouped

# BAM A with several -b BED files, -bed -wo -names.
# targets: multiple -b -names bam A | m.bam two beds | hand
begin 062
printf 'chr1\t110\t130\tq1\n' > "$TMPD/mq1.bed"
printf 'chr1\t120\t125\tq2\n' > "$TMPD/mq2.bed"
run_txt -a $M -b "$TMPD/mq1.bed" "$TMPD/mq2.bed" -names one two -bed -wo
lines "$E_FWD one chr1 110 130 q1 20" "$E_FWD two chr1 120 125 q2 5"
checkgrouped

###########################################################################
# -split with BAM A (blocks from CIGAR: M/=/X and D extend a block, N ends a block, I/S/H/P are ignored)
###########################################################################

# -split: a B inside the N gap of the spliced read (and inside the N gap of dn) is not a hit; without -split both are hits.
# targets: -split N gap | all.bam gap.bed | hand
begin 063
printf 'chr1\t550\t560\nchr2\t125\t135\n' > "$TMPD/gap.bed"
run_bam -a $A -b "$TMPD/gap.bed" -split
none
check "$TMPD/obs" "$TMPD/exp"

# Same B without -split overlaps both spans.
# targets: no -split control | all.bam gap.bed | hand
begin 064
run_bam -a $A -b "$TMPD/gap.bed"
lines "$SPL" "$DN"
check "$TMPD/obs" "$TMPD/exp"

# -split -v with the gap B reports every read including the spliced ones, since none has a block hit.
# targets: -split -v | all.bam gap.bed | hand
begin 065
run_bam -a $A -b "$TMPD/gap.bed" -split -v
lines "$FWD" "$PLC" "$REV" "$PE1" "$PE2" "$SPL" "$DEL" "$CLI" "$C2" "$DN" "$XE" "$UNM"
check "$TMPD/obs" "$TMPD/exp"

# -split with a B on the second block of spl only (chr1 605-615 overlaps [610,620)).
# targets: -split | all.bam tail.bed | hand
begin 066
printf 'chr1\t605\t615\n' > "$TMPD/tail.bed"
run_bam -a $A -b "$TMPD/tail.bed" -split
lines "$SPL"
check "$TMPD/obs" "$TMPD/exp"

# -split: a B lying entirely inside the 20 bp deletion of del still overlaps its block ([700,740) includes the D); same for dn's D.
# targets: -split D | all.bam delgap.bed | hand
begin 067
printf 'chr1\t715\t725\nchr2\t108\t112\n' > "$TMPD/delgap.bed"
run_bam -a $A -b "$TMPD/delgap.bed" -split
lines "$DEL" "$DN"
check "$TMPD/obs" "$TMPD/exp"

# -split: soft clip, insertion and hard clip are not part of a block: B at 820-830 (right after cli's end) does not overlap even though cli has 2+3+2 non-reference bases.
# targets: -split S I H | all.bam clip.bed | hand
begin 068
printf 'chr1\t820\t830\nchr1\t795\t800\n' > "$TMPD/clip.bed"
run_bam -a $A -b "$TMPD/clip.bed" -split
none
check "$TMPD/obs" "$TMPD/exp"

# -split with the default b.bed: overlapping reads are those with a block overlap (spl via 605-615 only, dn via 108-112 only).
# targets: -split | all.bam b.bed | hand
begin 069
run_bam -a $A -b $B -split
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -split -bed -wo: overlap bases are counted per block: spl overlaps only [610,620) so 5 (not 10 per hit), the gap hit disappears, dn keeps the D-side hit (4) and drops the N-gap hit.
# targets: -split -bed -wo | m.bam b.bed | hand
begin 070
run_txt -a $M -b $B -split -bed -wo
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 605 615 5" "$E_DEL chr1 715 725 10" "$E_DN chr2 108 112 4" "$E_XE chr2 305 308 3"
checkgrouped

# -split -bed -c counts only block hits: spl 1 (not 2), dn 1 (not 2).
# targets: -split -bed -c | m.bam b.bed | hand
begin 071
run_txt -a $M -b $B -split -bed -c
lines "$E_FWD 1" "$E_REV 1" "$E_PE1 0" "$E_PE2 1" "$E_SPL 1" "$E_DEL 1" "$E_CLI 0" "$E_C2 0" "$E_DN 1" "$E_XE 1"
check "$TMPD/obs" "$TMPD/exp"

# -split -bed -wao: alignments with no block hit (pe1, cli, c2) print a null B with 0 overlap; spl and dn print only their block hits.
# targets: -split -bed -wao | m.bam b.bed | hand
begin 072
run_txt -a $M -b $B -split -bed -wao
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE1 . -1 -1 0" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 605 615 5" "$E_DEL chr1 715 725 10" "$E_CLI . -1 -1 0" "$E_C2 . -1 -1 0" \
	"$E_DN chr2 108 112 4" "$E_XE chr2 305 308 3"
checkgrouped

# -split -bed -v: reads without a block hit: pe1, cli, c2 as before; spl and dn still have block hits.
# targets: -split -bed -v | m.bam b.bed | hand
begin 073
run_txt -a $M -b $B -split -bed -v
lines "$E_PE1" "$E_CLI" "$E_C2"
check "$TMPD/obs" "$TMPD/exp"

# -split -bed -wo sums the overlap over blocks: B chr1 505-615 overlaps spl [500,510) by 5 and [610,620) by 5 = 10; without -split it is 110.
# targets: -split -wo sum over blocks | m.bam sum.bed | hand
begin 074
printf 'chr1\t505\t615\n' > "$TMPD/sum.bed"
run_txt -a $M -b "$TMPD/sum.bed" -split -bed -wo
lines "$E_SPL chr1 505 615 10"
check "$TMPD/obs" "$TMPD/exp"

# Same B without -split: overlap is with the whole span, 110 bp.
# targets: -wo no split control | m.bam sum.bed | hand
begin 075
run_txt -a $M -b "$TMPD/sum.bed" -bed -wo
lines "$E_SPL chr1 505 615 110"
check "$TMPD/obs" "$TMPD/exp"

# Three blocks (tri [950,955) [965,970) [980,985)) against B 952-982: 3 + 5 + 2 = 10 bases.
# targets: -split -wo three blocks | cigar.bam | hand
begin 076
printf 'chr1\t952\t982\n' > "$TMPD/tri.bed"
run_txt -a "$TMPD/cigar.bam" -b "$TMPD/tri.bed" -split -bed -wo
lines "chr1 950 985 tri 9 + 950 985 0,0,0 3 5,5,5, 0,15,30, chr1 952 982 10"
check "$TMPD/obs" "$TMPD/exp"

# -split -f: fraction of A is over the summed block length, here 20 (spl) : 10/20 = 0.5 passes -f 0.5.
# targets: -split -f | m.bam sum.bed | hand
begin 077
run_txt -a $M -b "$TMPD/sum.bed" -split -f 0.5 -bed -wo
lines "$E_SPL chr1 505 615 10"
check "$TMPD/obs" "$TMPD/exp"

# -split -f 0.55 rejects it (10/20 = 0.5).
# targets: -split -f | m.bam sum.bed | hand
begin 078
run_txt -a $M -b "$TMPD/sum.bed" -split -f 0.55 -bed -wo
none
check "$TMPD/obs" "$TMPD/exp"

# Without -split the same -f is over the 120 bp span: 110/120 passes -f 0.9 and fails -f 0.95.
# targets: -f no split | m.bam sum.bed | hand
begin 079
run_txt -a $M -b "$TMPD/sum.bed" -f 0.9 -bed -wo
lines "$E_SPL chr1 505 615 110"
check "$TMPD/obs" "$TMPD/exp"

# Control for the previous test: -f 0.95 without -split rejects it.
# targets: -f no split | m.bam sum.bed | hand
begin 080
run_txt -a $M -b "$TMPD/sum.bed" -f 0.95 -bed -wo
none
check "$TMPD/obs" "$TMPD/exp"

# -split -f: a deletion counts toward the block length. del [700,740) has block length 40, B 715-725 overlaps 10: 10/40 = 0.25 passes -f 0.25 and fails -f 0.3 (if the D did not belong to the block, B 715-725 would fall in a gap and there would be no hit).
# targets: -split -f D | m.bam | hand
begin 081
printf 'chr1\t715\t725\n' > "$TMPD/d.bed"
run_txt -a $M -b "$TMPD/d.bed" -split -f 0.25 -bed -wo
lines "$E_DEL chr1 715 725 10"
check "$TMPD/obs" "$TMPD/exp"

# Companion: -f 0.3 rejects it.
# targets: -split -f D | m.bam | hand
begin 082
run_txt -a $M -b "$TMPD/d.bed" -split -f 0.3 -bed -wo
none
check "$TMPD/obs" "$TMPD/exp"

# -split -f: S, I and H are not in the block length. cli has block length 20; B 805-815 overlaps 10: 0.5 passes -f 0.5 (it would be 10/27 if the clips counted).
# targets: -split -f S I H | m.bam | hand
begin 083
printf 'chr1\t805\t815\n' > "$TMPD/cl.bed"
run_txt -a $M -b "$TMPD/cl.bed" -split -f 0.5 -bed -wo
lines "$E_CLI chr1 805 815 10"
check "$TMPD/obs" "$TMPD/exp"

# -split -f over three blocks: tri has 15 block bases, B 952-982 covers 10: 10/15 = 0.667 passes -f 0.6 and fails -f 0.7.
# targets: -split -f | cigar.bam tri.bed | hand
begin 084
run_txt -a "$TMPD/cigar.bam" -b "$TMPD/tri.bed" -split -f 0.6 -bed -wo
lines "chr1 950 985 tri 9 + 950 985 0,0,0 3 5,5,5, 0,15,30, chr1 952 982 10"
check "$TMPD/obs" "$TMPD/exp"

# Companion: -f 0.7 rejects it.
# targets: -split -f | cigar.bam tri.bed | hand
begin 085
run_txt -a "$TMPD/cigar.bam" -b "$TMPD/tri.bed" -split -f 0.7 -bed -wo
none
check "$TMPD/obs" "$TMPD/exp"

# -split -s / -S with BAM A (spl is +): B t1 (+) and t2 (-) over the second block.
# targets: -split -s -S | m.bam strand.bed | hand
begin 086
printf 'chr1\t605\t615\tt1\t0\t+\nchr1\t605\t615\tt2\t0\t-\n' > "$TMPD/strand.bed"
run_txt -a $M -b "$TMPD/strand.bed" -split -s -bed -wo
lines "$E_SPL chr1 605 615 t1 0 + 5"
check "$TMPD/obs" "$TMPD/exp"

# Companion: -S selects t2.
# targets: -split -S | m.bam strand.bed | hand
begin 087
run_txt -a $M -b "$TMPD/strand.bed" -split -S -bed -wo
lines "$E_SPL chr1 605 615 t2 0 - 5"
check "$TMPD/obs" "$TMPD/exp"

# -split with -sorted equals the unsorted result (BAM output).
# targets: -split -sorted | all.bam b.bed | hand
begin 088
run_bam -a $A -b $B -split -sorted
lines "$FWD" "$REV" "$PE2" "$SPL" "$DEL" "$DN" "$XE"
check "$TMPD/obs" "$TMPD/exp"

# -split -sorted -bed -wo equals the unsorted -split -wo lines.
# targets: -split -sorted -wo | m.bam b.bed | hand
begin 089
run_txt -a $M -b $B -split -sorted -bed -wo
lines "$E_FWD chr1 120 130 10" "$E_REV chr1 230 250 10" "$E_PE2 chr1 405 415 10" \
	"$E_SPL chr1 605 615 5" "$E_DEL chr1 715 725 10" "$E_DN chr2 108 112 4" "$E_XE chr2 305 308 3"
checkgrouped

###########################################################################
# -split with BED12 on either side
###########################################################################

# BED12 B g: blocks [100,120) [300,310) [610,620). With -split a BAM A hits it only where blocks meet: fwd 20, pe1 10, spl (block [610,620)) 10; rev and pe2 lie between the blocks.
# targets: -split bed12 B -wo | m.bam g.bed | hand
begin 090
printf 'chr1\t100\t620\tg\t0\t+\t100\t620\t0\t3\t20,10,10,\t0,200,510,\n' > "$TMPD/g.bed"
G="chr1 100 620 g 0 + 100 620 0 3 20,10,10, 0,200,510,"
run_txt -a $M -b "$TMPD/g.bed" -split -bed -wo
lines "$E_FWD $G 20" "$E_PE1 $G 10" "$E_SPL $G 10"
checkgrouped

# Same without -split: overlap is with the whole 520 bp span, so rev and pe2 are hits and spl overlaps 120.
# targets: bed12 B -wo no split | m.bam g.bed | hand
begin 091
run_txt -a $M -b "$TMPD/g.bed" -bed -wo
lines "$E_FWD $G 50" "$E_REV $G 40" "$E_PE1 $G 30" "$E_PE2 $G 30" "$E_SPL $G 120"
checkgrouped

# BAM output with BED12 B and -split: fwd, pe1, spl.
# targets: -split bed12 B bam output | all.bam g.bed | hand
begin 092
run_bam -a $A -b "$TMPD/g.bed" -split
lines "$FWD" "$PE1" "$SPL"
check "$TMPD/obs" "$TMPD/exp"

# -split -F with BED12 B: spl overlaps 10 of g's 40 block bases (0.25), so -F 0.25 passes (spl alone as A).
# targets: -split -F | spl only | hand
begin 093
mkbam spl1 "$SPL"
run_txt -a "$TMPD/spl1.bam" -b "$TMPD/g.bed" -split -F 0.25 -bed -wo
lines "$E_SPL $G 10"
check "$TMPD/obs" "$TMPD/exp"

# Companion: -F 0.3 rejects it.
# targets: -split -F | spl only | hand
begin 094
run_txt -a "$TMPD/spl1.bam" -b "$TMPD/g.bed" -split -F 0.3 -bed -wo
none
check "$TMPD/obs" "$TMPD/exp"

# -split -f 0.5 with BED12 B over all reads: only spl keeps 10/20 of its block bases; fwd 20/50 and pe1 10/30 fall short.
# targets: -split -f | m.bam g.bed | hand
begin 095
run_txt -a $M -b "$TMPD/g.bed" -split -f 0.5 -bed -wo
lines "$E_SPL $G 10"
check "$TMPD/obs" "$TMPD/exp"

# BED12 A with BAM B, -split -wo: A blocks [100,120) [300,310) [610,620) against alignment blocks: fwd 20, pe1 10, spl 10.
# targets: -split bed12 A, bam B | g.bed m.bam | hand
begin 096
run_txt -a "$TMPD/g.bed" -b $M -split -wo
lines "$G chr1 100 150 fwd 30 + 20" "$G chr1 300 330 pe1/1 40 + 10" "$G chr1 500 620 spl 50 + 10"
checkgrouped

# BED12 A with BAM B without -split: overlaps are with the 520 bp span.
# targets: bed12 A, bam B | g.bed m.bam | hand
begin 097
run_txt -a "$TMPD/g.bed" -b $M -wo
lines "$G chr1 100 150 fwd 30 + 50" "$G chr1 200 240 rev 20 - 40" "$G chr1 300 330 pe1/1 40 + 30" \
	"$G chr1 400 430 pe2/2 41 - 30" "$G chr1 500 620 spl 50 + 120"
checkgrouped

# BED12 A with BAM B, -split -c: counts block hits (3), not span hits (5).
# targets: -split -c bed12 A | g.bed m.bam | hand
begin 098
run_txt -a "$TMPD/g.bed" -b $M -split -c
lines "$G 3"
check "$TMPD/obs" "$TMPD/exp"

# BED12 A with BAM B and -sorted -split: same lines as unsorted (B limited to chr1).
# targets: -split -sorted | g.bed m1.bam | hand
begin 099
run_txt -a "$TMPD/g.bed" -b $M1 -split -sorted -wo
lines "$G chr1 100 150 fwd 30 + 20" "$G chr1 300 330 pe1/1 40 + 10" "$G chr1 500 620 spl 50 + 10"
checkgrouped

# BED12 A with BED3 B, -split -wo: A blocks [100,120) [280,300); B 110-130 overlaps block 1 by 10, B 150-270 lies in the gap (no hit), B 290-400 overlaps block 2 by 10.
# targets: -split bed12 A bed B -wo | bed12 + bed | hand
begin 100
printf 'chr1\t100\t300\tqa\t0\t+\t100\t300\t0\t2\t20,20,\t0,180,\n' > "$TMPD/qa.bed"
QA="chr1 100 300 qa 0 + 100 300 0 2 20,20, 0,180,"
printf 'chr1\t110\t130\nchr1\t150\t270\nchr1\t290\t400\n' > "$TMPD/qb.bed"
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -split -wo
lines "$QA chr1 110 130 10" "$QA chr1 290 400 10"
checkgrouped

# Same without -split: overlaps with the whole [100,300): 20, 120, 10.
# targets: bed12 A no split | bed12 + bed | hand
begin 101
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -wo
lines "$QA chr1 110 130 20" "$QA chr1 150 270 120" "$QA chr1 290 400 10"
checkgrouped

# BED12 A with -split -u: written once although two blocks have hits.
# targets: -split -u bed12 A | qa qb | hand
begin 102
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -split -u
lines "$QA"
check "$TMPD/obs" "$TMPD/exp"

# -split -c with BED12 A.
# targets: -split -c bed12 A | qa qb | hand
begin 103
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -split -c
lines "$QA 2"
check "$TMPD/obs" "$TMPD/exp"

# -split -v with BED12 A: it has block hits, so nothing is reported.
# targets: -split -v bed12 A | qa qb | hand
begin 104
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -split -v
none
check "$TMPD/obs" "$TMPD/exp"

# -split -wa with BED12 A: the A record is written once per block hit.
# targets: -split -wa bed12 A | qa qb | hand
begin 105
run_txt -a "$TMPD/qa.bed" -b "$TMPD/qb.bed" -split -wa
lines "$QA" "$QA"
check "$TMPD/obs" "$TMPD/exp"

# -split -f with BED12 A: 10 of 40 block bases is 0.25. B 110-130 alone (no second block hit) passes -f 0.25 and fails -f 0.3.
# targets: -split -f bed12 A | qa one B | hand
begin 106
printf 'chr1\t110\t130\n' > "$TMPD/q1.bed"
run_txt -a "$TMPD/qa.bed" -b "$TMPD/q1.bed" -split -f 0.25 -wo
lines "$QA chr1 110 130 10"
check "$TMPD/obs" "$TMPD/exp"

# Companion: -f 0.3 rejects it.
# targets: -split -f bed12 A | qa one B | hand
begin 107
run_txt -a "$TMPD/qa.bed" -b "$TMPD/q1.bed" -split -f 0.3 -wo
none
check "$TMPD/obs" "$TMPD/exp"

# BED12 A against BED12 B, -split -wo: A blocks [100,120) [280,300), B (g) blocks [100,120) [300,310) [610,620): only the first blocks coincide, 20 bases.
# targets: -split bed12 A and B -wo | qa g | hand
begin 108
run_txt -a "$TMPD/qa.bed" -b "$TMPD/g.bed" -split -wo
lines "$QA $G 20"
check "$TMPD/obs" "$TMPD/exp"

# BED12 A against BED12 B without -split: whole spans [100,300) and [100,620): 200.
# targets: bed12 A and B -wo no split | qa g | hand
begin 109
run_txt -a "$TMPD/qa.bed" -b "$TMPD/g.bed" -wo
lines "$QA $G 200"
check "$TMPD/obs" "$TMPD/exp"

###########################################################################
# CRAM (reference given locally; REF_PATH and REF_CACHE point at empty places so htslib cannot reach a network)
###########################################################################

# CRAM -a with a BED B: output is CRAM when CRAM_REFERENCE is set; the read overlapping 10120-10130 is the second mate (flag 147, 3S97M, [10046,10143)); the first mate [10003,10103) does not reach it.
# targets: cram A | a.cram cb.bed | hand
begin 110
printf 'chr1\t10120\t10130\n' > "$TMPD/cb.bed"
CRAM_REFERENCE=test_ref.fa $BT intersect -a a.cram -b "$TMPD/cb.bed" 2> /dev/null | $htsutil viewcramrecords - test_ref.fa > "$TMPD/obs"
printf '%s\n' "FCC1MK2ACXX:2:2110:4301:28831#	147	chr1	10047	0	3S97M	=	10004	-140	AACCCTACCCCTACCCCTAACCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCCTACCCTAACCCTAACCCTAACCCTAACCCTAACCC	###?<55-@?250&A?882(@?795&B?8;;/?8('9'A?;8?8C;(<A@FB-/7'F?((@0D:)*9JIIHFCJJJIHFJJJIHEJJHHFFHFFFDDCCB	PG:Z:novoalign	AM:i:2	SM:i:70	MQ:i:0	PQ:i:219	UQ:i:194	AS:i:194	MD:Z:4A5A11A5A5A5A5A5A5A3A34	NM:i:10	RG:Z:NCH411GBM_CD133low" > "$TMPD/exp"
check "$TMPD/obs" "$TMPD/exp"

# CRAM -a with -bed (no reference needed): the second mate as BED12 trimmed to the overlap; name gets /2, MAPQ 0, strand -.
# targets: cram A -bed | a.cram cb.bed | hand
begin 111
$BT intersect -a a.cram -b "$TMPD/cb.bed" -bed 2> /dev/null > "$TMPD/obs"
lines "chr1 10120 10130 FCC1MK2ACXX:2:2110:4301:28831#/2 0 - 10046 10143 0,0,0 1 97, 0,"
check "$TMPD/obs" "$TMPD/exp"

# CRAM -a with -bed -wa -wo: the full read as BED12, the B record and 10 bases of overlap.
# targets: cram A -bed -wo | a.cram cb.bed | hand
begin 112
$BT intersect -a a.cram -b "$TMPD/cb.bed" -bed -wo 2> /dev/null > "$TMPD/obs"
lines "chr1 10046 10143 FCC1MK2ACXX:2:2110:4301:28831#/2 0 - 10046 10143 0,0,0 1 97, 0, chr1 10120 10130 10"
check "$TMPD/obs" "$TMPD/exp"

# CRAM as B: BED A with -wa -wb prints the CRAM reads as BED6 (flag 147: /2, strand -).
# targets: cram B -wa -wb | cb.bed a.cram | hand
begin 113
$BT intersect -a "$TMPD/cb.bed" -b a.cram -wa -wb 2> /dev/null > "$TMPD/obs"
lines "chr1 10120 10130 chr1 10046 10143 FCC1MK2ACXX:2:2110:4301:28831#/2 0 -"
check "$TMPD/obs" "$TMPD/exp"

# CRAM -a with -v: the read not overlapping B is the first mate (flag 99, [10003,10103)), reported as a CRAM/BAM record.
# targets: cram A -v | a.cram cb.bed | hand
begin 114
$BT intersect -a a.cram -b "$TMPD/cb.bed" -v 2> /dev/null | $htsutil viewbamrecords - 2> /dev/null | cut -f 1-9 > "$TMPD/obs"
lines "FCC1MK2ACXX:2:2110:4301:28831# 99 chr1 10004 0 100M = 10047 140"
check "$TMPD/obs" "$TMPD/exp"

finish
