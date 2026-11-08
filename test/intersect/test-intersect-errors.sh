. ./helpers-intersect.sh errors
set -e

# begin_n: begin with the next sequential test id.
NEXT=0
begin_n()
{
	flush_pending
	NEXT=$((NEXT + 1))
	begin "$(printf '%03d' $NEXT)"
}

# Local helpers.  Every fixture is a hand-written file in $TMPD.

# tabs: printf wrapper so the fixtures below read like BED.
mk()
{
	local f=$1
	shift
	: > "$TMPD/$f"
	local l
	for l in "$@"; do
		printf '%s\n' "$l" | tr ' ' '\t' >> "$TMPD/$f"
	done
}

# fail1 <substring> <command...>: exit status exactly 1, nothing on stdout,
# <substring> on stderr.  With USAGE=1 the usage line must be on stderr too
# (bedtools prints the help before the error for argument errors).
USAGE=0
ALSO=""
fail1()
{
	local want=$1 ok=1
	shift
	set +e
	"$@" > "$TMPD/f1.stdout" 2> "$TMPD/f1.stderr"
	local status=$?
	set -e
	[ $status -eq 1 ] || ok=0
	[ ! -s "$TMPD/f1.stdout" ] || ok=0
	grep -F -e "$want" "$TMPD/f1.stderr" > /dev/null || ok=0
	if [ -n "$ALSO" ]; then
		grep -F -e "$ALSO" "$TMPD/f1.stderr" > /dev/null || ok=0
	fi
	if [ $USAGE -eq 1 ]; then
		grep -F -e "Usage:   bedtools intersect" "$TMPD/f1.stderr" > /dev/null || ok=0
	fi
	if [ $ok -eq 1 ]; then
		echo ok
	else
		FAILURES=$((FAILURES + 1))
		echo fail
		echo "      exit $status; stdout bytes: $(wc -c < ""$TMPD/f1.stdout""); stderr:"
		sed 's/^/        /' "$TMPD/f1.stderr" | head -20
	fi
}

# argerr <substring> <command...>: fail1 for errors that come back through
# the context (message plus usage text).
argerr()
{
	USAGE=1
	fail1 "$@"
	USAGE=0
}

# runerr <substring> <command...>: fail1 for errors that exit from inside
# the tool (no usage text expected, but exit status is still 1).
runerr()
{
	USAGE=0
	fail1 "$@"
}

# A test that has several assertions decides pass/fail before printing:
# pcheck/pfail/warned/quiet record into PEND_* and flush_pending prints one
# ok or fail line.  begin_n and finish_all flush the previous test.
PEND=0
PEND_FAIL=0
PEND_MSG=""
pend_start() { PEND=1; PEND_FAIL=0; PEND_MSG=""; }
pfail() { PEND_FAIL=1; PEND_MSG="$PEND_MSG
      $1"; }
flush_pending()
{
	if [ $PEND -eq 1 ]; then
		if [ $PEND_FAIL -eq 0 ]; then
			echo ok
		else
			FAILURES=$((FAILURES + 1))
			echo fail
			printf '%s\n' "$PEND_MSG"
		fi
		PEND=0
	fi
}
finish_all() { flush_pending; finish; }

# pcheck <observed file> <expected file>: byte comparison, recorded.
pcheck()
{
	if ! diff "$1" "$2" > "$TMPD/diff.out" 2>&1; then
		pfail "output differs: $(cat ""$TMPD/diff.out"")"
	fi
}

# okrun <expected file> <command...>: status 0 and stdout equals the file,
# stderr is kept in $TMPD/ok.stderr.  Followed by optional warned/quiet.
okrun()
{
	local exp=$1
	shift
	pend_start
	set +e
	"$@" > "$TMPD/ok.stdout" 2> "$TMPD/ok.stderr"
	local status=$?
	set -e
	if [ $status -ne 0 ]; then
		pfail "exit $status; stderr: $(head -5 ""$TMPD/ok.stderr"")"
		return 0
	fi
	pcheck "$TMPD/ok.stdout" "$exp"
}

# warned <substring>: stderr of the last run must contain <substring>.
warned()
{
	if ! grep -F -e "$1" "$TMPD/ok.stderr" > /dev/null; then
		pfail "missing on stderr: $1 (stderr: $(head -5 ""$TMPD/ok.stderr""))"
	fi
}

# quiet: stderr of the last run must be empty.
quiet()
{
	if [ -s "$TMPD/ok.stderr" ]; then
		pfail "unexpected stderr: $(head -5 ""$TMPD/ok.stderr"")"
	fi
}

# Small shared fixtures.
mk a3.bed "chr1 10 20" "chr1 30 40"
mk b3.bed "chr1 15 25"
mk a6.bed "chr1 10 20 x 1 +" "chr1 30 40 y 1 -"
mk b6.bed "chr1 15 25 z 1 +"
mk a5.bed "chr1 10 20 x 1"
mk one.bed "chr1 15 20"
mk e_one.bed "chr1 10 20"
: > "$TMPD/empty.bed"
mk exp_overlap "chr1 15 20"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t12\t.\tA\tC\t.\t.\t.\n' > "$TMPD/v.vcf"

# SAM -> BAM: r1 is chr1:10-20 (+), r2 is chr1:100-110 (-).
printf '@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nr1\t0\tchr1\t11\t40\t10M\t*\t0\t0\tAAAAAAAAAA\tIIIIIIIIII\nr2\t16\tchr1\t101\t40\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII\n' > "$TMPD/q.sam"
$htsutil samtobam "$TMPD/q.sam" "$TMPD/q.bam"
mk bq.bed "chr1 15 25 b1"

# ======================================================================
# Argument errors: message on stderr, usage printed, exit status 1.
# ======================================================================

# unrecognized option is rejected with its name.
begin_n
argerr "ERROR: Unrecognized parameter: -bogus" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -bogus

# every unrecognized option is listed, not only the first.
begin_n
ALSO="Unrecognized parameter: -bogus"
argerr "Unrecognized parameter: -zzz" "$BT" intersect -a "$TMPD/a3.bed" -bogus -b "$TMPD/b3.bed" -zzz
ALSO=""

# -a as the last argument has no value.
begin_n
argerr "ERROR: -a option given, but no query file specified." "$BT" intersect -b "$TMPD/b3.bed" -a

# -b as the last argument has no value.
begin_n
argerr "ERROR: -b option given, but no database file specified." "$BT" intersect -a "$TMPD/a3.bed" -b

# -g as the last argument has no value.
begin_n
argerr "ERROR: -g option given, but no genome file specified." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -g

# -iobuf as the last argument has no value.
begin_n
argerr "ERROR: -iobuf option given, but size of input buffer not specified." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf

# -names as the last argument has no value.
begin_n
argerr "but no database names specified." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -names

# the obsolete -abam as the last argument has no value.
begin_n
argerr "ERROR: -abam option given, but no query BAM file specified." "$BT" intersect -b "$TMPD/b3.bed" -abam

# -iobuf below the 8 byte minimum.
begin_n
argerr "ERROR: specified buffer size is too small." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf 7

# -iobuf with a suffix other than K, M, G names the suffix.
begin_n
argerr "Unrecognized memory buffer size suffix 'x' given." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf 12x

# -iobuf with a valid suffix but a non-numeric body.
begin_n
argerr "ERROR: argument passed to -iobuf is not numeric." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf xK

# -iobuf 8 (the minimum) and a K suffix leave the result unchanged.
begin_n
okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf 8
begin_n
okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -iobuf 1K

# -f must be in (0, 1]: zero, above one, and a non-number (read as 0).
begin_n
argerr "ERROR: -f must be in the range (0.0, 1.0]." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -f 0
begin_n
argerr "ERROR: -f must be in the range (0.0, 1.0]." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -f 1.5
begin_n
argerr "ERROR: -f must be in the range (0.0, 1.0]." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -f abc

# -F must be in (0, 1].
begin_n
argerr "ERROR: -F must be in the range (0.0, 1.0]." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -F 0
begin_n
argerr "ERROR: -F must be in the range (0.0, 1.0]." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -F 2

# -f as the last argument has no value, so it is left as an unrecognized option.
begin_n
argerr "Unrecognized parameter: -f" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -f

# -f 1 is accepted: A chr1:10-20 inside B chr1:5-30 is reported, and when B
# stops at 19 (9 of 10 bases) it is not.
begin_n
mk fa "chr1 10 20"
mk fb_full "chr1 5 30"
mk fb_part "chr1 5 19"
okrun "$TMPD/fa" "$BT" intersect -a "$TMPD/fa" -b "$TMPD/fb_full" -f 1
begin_n
okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/fa" -b "$TMPD/fb_part" -f 1

# -F 1 is accepted: B chr1:12-18 inside A chr1:10-20 is reported (as the
# overlap 12-18); B chr1:5-30 is not covered by A.
begin_n
mk exp_f1 "chr1 12 18"
mk fb_in "chr1 12 18"
okrun "$TMPD/exp_f1" "$BT" intersect -a "$TMPD/fa" -b "$TMPD/fb_in" -F 1
begin_n
okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/fa" -b "$TMPD/fb_full" -F 1

# -e with -F but without -f is an error.
begin_n
argerr "ERROR: -e must be used with -f." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -F 0.5 -e

# -r with -F is an error.
begin_n
argerr "ERROR: -r must be used solely with -f." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -F 0.5 -r

# -r with both -f and -F is an error.
begin_n
argerr "ERROR: -r must be used solely with -f." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -f 0.5 -F 0.5 -r

# -u with -v.
begin_n
argerr "ERROR: request either -u for anyHit OR -v for noHit, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -u -v

# -c with -u.
begin_n
argerr "ERROR: request either -c for writeCount OR -u for anyHit, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -c -u

# -c with -wb.
begin_n
argerr "ERROR: request either -c for writeCount OR -wb for writeB, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -c -wb

# -c with -filenames and with -names.
begin_n
argerr "ERROR: -c cannot be used with -filenames or -names. Use -C." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -c -filenames
begin_n
argerr "ERROR: -c cannot be used with -filenames or -names. Use -C." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -c -names x

# -C with two databases and names: one line per (A record, database); A
# chr1:10-20 overlaps b3 (15-25) once and b_far never, chr1:30-40 neither.
begin_n
mk b_far "chr1 500 600"
mk exp_C "chr1 10 20 x 1" "chr1 10 20 y 0" "chr1 30 40 x 0" "chr1 30 40 y 0"
okrun "$TMPD/exp_C" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" "$TMPD/b_far" -C -names x y

# -wo with -wa, -wb, -c and -u.
begin_n
argerr "ERROR: request either -wa for writeA OR -wo for writeOverlap, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -wo -wa
begin_n
argerr "ERROR: request either -wb for writeB OR -wo for writeOverlap, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -wo -wb
begin_n
argerr "ERROR: request either -c for writeCount OR -wo for writeOverlap, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -wo -c
begin_n
argerr "ERROR: request either -u for anyHit OR -wo for writeOverlap, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -wo -u

# -wao is a form of -wo internally, so -wao with -wa reports the -wo conflict.
begin_n
argerr "ERROR: request either -wa for writeA OR -wo for writeOverlap, not both." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -wao -wa

# -s with -S.
begin_n
argerr "ERROR: request -s for sameStrand, or -S for diffStrand, not both." "$BT" intersect -a "$TMPD/a6.bed" -b "$TMPD/b6.bed" -s -S

# stranded request on a BED5 file (no strand column).
begin_n
argerr "stranded intersect requested, but input file $TMPD/a5.bed does not have strands." "$BT" intersect -a "$TMPD/a5.bed" -b "$TMPD/b6.bed" -s
begin_n
argerr "does not have strands." "$BT" intersect -a "$TMPD/a6.bed" -b "$TMPD/a5.bed" -S

# stranded request with a VCF file, as B and as A.
begin_n
argerr "stranded intersect not supported for VCF file $TMPD/v.vcf." "$BT" intersect -a "$TMPD/a6.bed" -b "$TMPD/v.vcf" -s
begin_n
argerr "stranded intersect not supported for VCF file $TMPD/v.vcf." "$BT" intersect -a "$TMPD/v.vcf" -b "$TMPD/b6.bed" -S

# -names count mismatches: more names than databases, fewer names than databases.
begin_n
argerr "ERROR: Number of database name tags given does not match number of databases." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -names x y
begin_n
argerr "ERROR: Number of database name tags given does not match number of databases." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" "$TMPD/b3.bed" -names x

# no -b, and no arguments at all.  (A run with -b and no -a is left out: bedtools crashes on it.)
begin_n
argerr "ERROR: query and database files not specified." "$BT" intersect -a "$TMPD/a3.bed"
begin_n
argerr "ERROR: No input file given. Exiting." "$BT" intersect

# nonexistent inputs exit with status 1 and name the file.
begin_n
runerr "Error: Unable to open file $TMPD/nonexistent.bed" "$BT" intersect -a "$TMPD/nonexistent.bed" -b "$TMPD/b3.bed"
begin_n
runerr "Error: Unable to open file $TMPD/nonexistent.bed" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/nonexistent.bed"

# genome file problems exit with status 1 as soon as -g is parsed (with -sorted).
begin_n
runerr "Error: Can't open genome file $TMPD/nonexistent.txt" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -g "$TMPD/nonexistent.txt" -sorted
begin_n
: > "$TMPD/g_empty.txt"
runerr "has no valid entries" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -g "$TMPD/g_empty.txt" -sorted
begin_n
mk g_zero "chr1 0"
runerr "Error: chr1 has length equal to 0." "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -g "$TMPD/g_zero" -sorted

# ======================================================================
# BAM query restrictions and warnings.
# ======================================================================

# -c with a BAM query needs -bed.
begin_n
argerr "ERROR: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option." "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -c

# -wo and -wao with a BAM query need -bed.
begin_n
argerr "ERROR: writeAllOverlap option is not valid with BAM query input, unless bed output is specified with -bed option." "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -wo
begin_n
argerr "writeAllOverlap option is not valid with BAM query input" "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -wao

# bamrun <expected rows file or ""> <command...>: BAM output; status must be
# 0; the records (via htsutil viewbamrecords) must equal the file unless the
# file name is empty.  stderr is kept for warned/quiet.
bamrun()
{
	local exp=$1
	shift
	pend_start
	set +e
	"$@" > "$TMPD/o.bam" 2> "$TMPD/ok.stderr"
	local status=$?
	set -e
	if [ $status -ne 0 ]; then
		pfail "exit $status; stderr: $(head -5 ""$TMPD/ok.stderr"")"
		return 0
	fi
	if [ -n "$exp" ]; then
		$htsutil viewbamrecords "$TMPD/o.bam" > "$TMPD/ok.stdout"
		pcheck "$TMPD/ok.stdout" "$exp"
	fi
}

# -c -bed with a BAM query: r1 overlaps bq.bed once, r2 does not.  The first
# six BED columns are the alignment's; the last column is the count.
begin_n
mk exp_bamc "chr1 10 20 r1 40 + 1" "chr1 100 110 r2 40 - 0"
pend_start
set +e
"$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -c -bed > "$TMPD/o.txt" 2> "$TMPD/ok.stderr"
[ $? -eq 0 ] || pfail "nonzero exit"
set -e
cut -f1-6,13 "$TMPD/o.txt" > "$TMPD/obs"
pcheck "$TMPD/obs" "$TMPD/exp_bamc"

# -wo -bed with a BAM query: alignment columns, B columns, overlap 5 (15-20).
begin_n
mk exp_bamwo "chr1 10 20 r1 chr1 15 25 b1 5"
pend_start
set +e
"$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -wo -bed > "$TMPD/o.txt" 2> "$TMPD/ok.stderr"
[ $? -eq 0 ] || pfail "nonzero exit"
set -e
cut -f1-4,13-17 "$TMPD/o.txt" > "$TMPD/obs"
pcheck "$TMPD/obs" "$TMPD/exp_bamwo"

# -wb with BAM output warns and is ignored: only r1 overlaps and is written.
begin_n
mk exp_r1 "r1 0 chr1 11 40 10M * 0 0 AAAAAAAAAA IIIIIIIIII"
bamrun "$TMPD/exp_r1" "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -wb
warned "WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option."

# -loj with BAM output warns; only the exit status and the warning are asserted.
begin_n
bamrun "" "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -loj
warned "WARNING: -wb and -loj are ignored with bam input"

# -wb with -bed does not warn, and writes A (as the overlap 15-20 without -wa)
# followed by the B record.
begin_n
mk exp_wbbed "chr1 15 20 r1 chr1 15 25 b1"
pend_start
set +e
"$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -wb -bed > "$TMPD/o.txt" 2> "$TMPD/ok.stderr"
[ $? -eq 0 ] || pfail "nonzero exit"
set -e
cut -f1-4,13-16 "$TMPD/o.txt" > "$TMPD/obs"
pcheck "$TMPD/obs" "$TMPD/exp_wbbed"
if grep -F WARNING "$TMPD/ok.stderr" > /dev/null; then pfail "unexpected warning"; fi

# -header with a BAM query warns; the BAM output is still written.
begin_n
bamrun "$TMPD/exp_r1" "$BT" intersect -a "$TMPD/q.bam" -b "$TMPD/bq.bed" -header
warned "WARNING: -header option is not valid for BAM input."

# ======================================================================
# Help.
# ======================================================================

# -h prints the usage to stderr, nothing to stdout, and exits 1.
begin_n
argerr "Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>" "$BT" intersect -h

# --help is the same.
begin_n
argerr "Summary: Report overlaps between two feature files." "$BT" intersect --help

# bare "bedtools intersect" prints the usage too.
begin_n
argerr "Tool:    bedtools intersect (aka intersectBed)" "$BT" intersect

# -h after real arguments still gives help rather than running.
begin_n
argerr "Usage:   bedtools intersect" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -h

# the usage names every documented option of each help block.
help_has()
{
	local opt missing=""
	"$BT" intersect -h > /dev/null 2> "$TMPD/help.txt" || true
	for opt in "$@"; do
		grep -F -e "$(printf '\t%s\t' "$opt")" "$TMPD/help.txt" > /dev/null || missing="$missing $opt"
	done
	if [ -z "$missing" ]; then
		echo ok
	else
		FAILURES=$((FAILURES + 1)); echo fail; echo "      missing from help:$missing"
	fi
}

begin_n
help_has -wa -wb -loj -wo -wao -u -c -C -v -ubam
begin_n
help_has -s -S -f -F -r -e -split -g -nonamecheck
begin_n
help_has -sorted
begin_n
help_has -names -filenames -sortout
begin_n
help_has -bed -header -nobuf -iobuf

# the help text ends with the BAM note.
begin_n
argerr "(1) When a BAM file is used for the A file, the alignment is retained" "$BT" intersect -h

# the alias prints the same usage and exits 1.
begin_n
argerr "Usage:   bedtools intersect [OPTIONS]" "$IBT" -h

# ======================================================================
# -sorted with sorted input: hand-worked results, checked in both modes.
# ======================================================================

mk sa "chr1 10 20" "chr1 30 40" "chr2 5 15" "chr3 1 10"
mk sb "chr1 15 35" "chr2 10 12" "chr2 100 200" "chr4 1 5"

# A1 (10-20) and A2 (30-40) overlap B1 (15-35) in 15-20 and 30-35; A3 (5-15)
# overlaps B2 (10-12); A4 on chr3 overlaps nothing.
mk exp_s_int "chr1 15 20" "chr1 30 35" "chr2 10 12"
mk exp_s_wawb "chr1 10 20 chr1 15 35" "chr1 30 40 chr1 15 35" "chr2 5 15 chr2 10 12"
mk exp_s_wa "chr1 10 20" "chr1 30 40" "chr2 5 15"
mk exp_s_c "chr1 10 20 1" "chr1 30 40 1" "chr2 5 15 1" "chr3 1 10 0"
mk exp_s_v "chr3 1 10"
mk exp_s_loj "chr1 10 20 chr1 15 35" "chr1 30 40 chr1 15 35" "chr2 5 15 chr2 10 12" "chr3 1 10 . -1 -1"
mk exp_s_wo "chr1 10 20 chr1 15 35 5" "chr1 30 40 chr1 15 35 5" "chr2 5 15 chr2 10 12 2"

for mode in "" "-sorted"; do
	begin_n
	okrun "$TMPD/exp_s_int" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" $mode
	begin_n
	okrun "$TMPD/exp_s_wawb" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -wa -wb $mode
	begin_n
	okrun "$TMPD/exp_s_wa" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -u $mode
	begin_n
	okrun "$TMPD/exp_s_c" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -c $mode
	begin_n
	okrun "$TMPD/exp_s_v" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -v $mode
	begin_n
	okrun "$TMPD/exp_s_loj" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -loj $mode
	begin_n
	okrun "$TMPD/exp_s_wo" "$BT" intersect -a "$TMPD/sa" -b "$TMPD/sb" -wo $mode
done

# Swapping A and B: B1 (15-35) overlaps A1 in 15-20 and A2 in 30-35, B2
# (10-12) overlaps A3 (5-15) in 10-12; chr3 and chr4 are each on one side only.
begin_n
okrun "$TMPD/exp_s_int" "$BT" intersect -a "$TMPD/sb" -b "$TMPD/sa" -sorted

# ======================================================================
# -sorted order checks on a single file.
# ======================================================================

mk b_ok "chr1 15 25"

# A record whose start is lower than the previous one on the same chromosome.
mk a_unsorted "chr1 100 200" "chr1 10 20"
begin_n
runerr "Error: Sorted input specified, but the file $TMPD/a_unsorted has the following out of order record" "$BT" intersect -a "$TMPD/a_unsorted" -b "$TMPD/b_ok" -sorted
# The same file as B.
begin_n
runerr "Error: Sorted input specified, but the file $TMPD/a_unsorted has the following out of order record" "$BT" intersect -a "$TMPD/b_ok" -b "$TMPD/a_unsorted" -sorted
# A chromosome that comes back after another one.
mk a_revisit "chr1 10 20" "chr2 10 20" "chr1 30 40"
begin_n
runerr "has the following out of order record" "$BT" intersect -a "$TMPD/a_revisit" -b "$TMPD/b_ok" -sorted
begin_n
runerr "has the following out of order record" "$BT" intersect -a "$TMPD/b_ok" -b "$TMPD/a_revisit" -sorted

# Equal starts are in order.
begin_n
mk a_tie "chr1 10 20" "chr1 10 30"
mk exp_tie "chr1 15 20" "chr1 15 25"
okrun "$TMPD/exp_tie" "$BT" intersect -a "$TMPD/a_tie" -b "$TMPD/b_ok" -sorted

# Without -sorted the same unsorted files are accepted: A 100-200 does not
# overlap B 15-25, A 10-20 does.
begin_n
okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/a_unsorted" -b "$TMPD/b_ok"
begin_n
okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/a_revisit" -b "$TMPD/b_ok"

# ======================================================================
# Chromosome order across files (no genome file).
# ======================================================================

# anum is in numeric chromosome order, blex in lexicographic order.
mk anum "chr1 10 20" "chr2 10 20" "chr10 10 20"
mk blex "chr1 10 20" "chr10 10 20" "chr2 10 20"

begin_n
runerr "ERROR: chromosome sort ordering for file " "$BT" intersect -a "$TMPD/anum" -b "$TMPD/blex" -sorted
begin_n
runerr "is inconsistent with other files. Record was:" "$BT" intersect -a "$TMPD/blex" -b "$TMPD/anum" -sorted

# Two files in the same numeric order agree with each other.
begin_n
okrun "$TMPD/anum" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/anum" -sorted
begin_n
okrun "$TMPD/blex" "$BT" intersect -a "$TMPD/blex" -b "$TMPD/blex" -sorted

# Unsorted mode does not care about chromosome order.
begin_n
okrun "$TMPD/anum" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/blex"
begin_n
mk exp_blex_a "chr1 10 20" "chr10 10 20" "chr2 10 20"
okrun "$TMPD/exp_blex_a" "$BT" intersect -a "$TMPD/blex" -b "$TMPD/anum"

# Once numeric order is established, a database chromosome the query does not
# have is an error without -g.
mk bm "chrM 1 2"
begin_n
runerr "ERROR: Database file $TMPD/bm contains chromosome chrM, but the query file does not." "$BT" intersect -a "$TMPD/anum" -b "$TMPD/bm" -sorted
begin_n
runerr "Please rerun with the -g option for a genome file." "$BT" intersect -a "$TMPD/anum" -b "$TMPD/bm" -sorted

# A query that starts at chrM (lexicographically after chr10) followed by a
# database in numeric order: the sweep assumed lexicographic order and finds
# the database is not.
begin_n
runerr "ERROR: Sort order was unspecified, and file $TMPD/anum is not sorted lexicographically." "$BT" intersect -a "$TMPD/bm" -b "$TMPD/anum" -sorted

# Chromosomes missing from one side with lexicographic order are fine:
# only chr2 is shared.
mk a12 "chr1 10 20" "chr2 10 20"
mk b23 "chr2 15 25" "chr3 10 20"
mk exp_chr2 "chr2 15 20"
begin_n
okrun "$TMPD/exp_chr2" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/b23" -sorted
begin_n
okrun "$TMPD/exp_chr2" "$BT" intersect -a "$TMPD/b23" -b "$TMPD/a12" -sorted
# -v keeps A's chr1 record, -loj keeps it with a NULL B.
begin_n
mk exp_a12_v "chr1 10 20"
okrun "$TMPD/exp_a12_v" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/b23" -sorted -v
begin_n
mk exp_a12_loj "chr1 10 20 . -1 -1" "chr2 10 20 chr2 15 25"
okrun "$TMPD/exp_a12_loj" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/b23" -sorted -loj

# ======================================================================
# Genome files with -sorted.
# ======================================================================

mk gnum "chr1 1000" "chr2 1000" "chr10 1000"
mk glex "chr1 1000" "chr10 1000" "chr2 1000"
mk gmiss "chr1 1000" "chr2 1000"
mk gall "chr1 1000" "chr2 1000" "chr10 1000" "chrM 1000"

# A genome file in the numeric order the files use.
begin_n
okrun "$TMPD/anum" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/anum" -sorted -g "$TMPD/gnum"
# Lexicographic files with a lexicographic genome file.
begin_n
okrun "$TMPD/blex" "$BT" intersect -a "$TMPD/blex" -b "$TMPD/blex" -sorted -g "$TMPD/glex"
# A genome file whose order disagrees with the files.
begin_n
runerr "has the following record with a different sort order than the genomeFile $TMPD/glex" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/anum" -sorted -g "$TMPD/glex"
begin_n
runerr "has the following record with a different sort order than the genomeFile $TMPD/gnum" "$BT" intersect -a "$TMPD/blex" -b "$TMPD/blex" -sorted -g "$TMPD/gnum"
# Files with chromosomes the genome file does not list.
begin_n
runerr "Error: requested chromosome chr10 does not exist in the genome file $TMPD/gmiss. Exiting." "$BT" intersect -a "$TMPD/anum" -b "$TMPD/anum" -sorted -g "$TMPD/gmiss"
begin_n
runerr "does not exist in the genome file" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/anum" -sorted -g "$TMPD/gmiss"
# A genome file with extra chromosomes is fine.
begin_n
okrun "$TMPD/anum" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/anum" -sorted -g "$TMPD/gall"

# -g supplies the order the sweep cannot guess: chrM after chr10 in numeric
# order.  A has chrM only; B has the numeric-order file.  Without -g this is
# the error of the chrM test above (no -g); with the genome file the two files are consistent and
# chrM 1-2 overlaps nothing (output empty), -v lists it.
begin_n
okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/bm" -b "$TMPD/anum" -sorted -g "$TMPD/gall"
begin_n
okrun "$TMPD/bm" "$BT" intersect -a "$TMPD/bm" -b "$TMPD/anum" -sorted -g "$TMPD/gall" -v

# Database chromosome the query lacks, numeric order, with -g: no error.
begin_n
okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/anum" -b "$TMPD/bm" -sorted -g "$TMPD/gall"

# Overlaps hand-worked under a non-lexicographic genome order (chr2, chr10, chr1).
mk gweird "chr2 1000" "chr10 1000" "chr1 1000"
mk wa "chr2 10 20" "chr10 10 20" "chr1 10 20" "chr1 50 60"
mk wb "chr2 15 25" "chr1 55 70"
begin_n
mk exp_w_int "chr2 15 20" "chr1 55 60"
okrun "$TMPD/exp_w_int" "$BT" intersect -a "$TMPD/wa" -b "$TMPD/wb" -sorted -g "$TMPD/gweird"
begin_n
mk exp_w_c "chr2 10 20 1" "chr10 10 20 0" "chr1 10 20 0" "chr1 50 60 1"
okrun "$TMPD/exp_w_c" "$BT" intersect -a "$TMPD/wa" -b "$TMPD/wb" -sorted -g "$TMPD/gweird" -c
begin_n
mk exp_w_v "chr10 10 20" "chr1 10 20"
okrun "$TMPD/exp_w_v" "$BT" intersect -a "$TMPD/wa" -b "$TMPD/wb" -sorted -g "$TMPD/gweird" -v
begin_n
mk exp_w_wo "chr2 10 20 chr2 15 25 5" "chr1 50 60 chr1 55 70 5"
okrun "$TMPD/exp_w_wo" "$BT" intersect -a "$TMPD/wa" -b "$TMPD/wb" -sorted -g "$TMPD/gweird" -wo
# -g with a record out of order inside a chromosome still reports the plain
# out-of-order error.
begin_n
runerr "has the following out of order record" "$BT" intersect -a "$TMPD/a_unsorted" -b "$TMPD/b_ok" -sorted -g "$TMPD/gnum"

# ======================================================================
# Empty and header-only inputs, both modes.
# ======================================================================

mk hdr_only "#header line"
printf 'track name=t\n' > "$TMPD/track_only"
mk exp_a12 "chr1 10 20" "chr2 10 20"
mk exp_a12_null "chr1 10 20 . -1 -1" "chr2 10 20 . -1 -1"
mk exp_a12_wao "chr1 10 20 . -1 -1 0" "chr2 10 20 . -1 -1 0"
mk exp_a12_c "chr1 10 20 0" "chr2 10 20 0"

for mode in "" "-sorted"; do
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/empty.bed" -b "$TMPD/a12" $mode
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/empty.bed" -b "$TMPD/empty.bed" $mode
	begin_n
	okrun "$TMPD/exp_a12" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode -v
	begin_n
	okrun "$TMPD/exp_a12_null" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode -loj
	begin_n
	okrun "$TMPD/exp_a12_wao" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode -wao
	begin_n
	okrun "$TMPD/exp_a12_c" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode -c
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/empty.bed" $mode -u
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/empty.bed" -b "$TMPD/a12" $mode -c
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/empty.bed" -b "$TMPD/a12" $mode -v
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/hdr_only" -b "$TMPD/a12" $mode
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/hdr_only" $mode
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/hdr_only" -b "$TMPD/hdr_only" $mode
	begin_n
	okrun "$TMPD/exp_a12" "$BT" intersect -a "$TMPD/a12" -b "$TMPD/track_only" $mode -v
	begin_n
	okrun "$TMPD/hdr_only" "$BT" intersect -a "$TMPD/hdr_only" -b "$TMPD/a12" $mode -header
done

# ======================================================================
# Name conventions.  Warnings go to stderr; the exit status stays 0.
# ======================================================================

mk nc_chr "chr1 10 20"
mk nc_bare "1 15 25"
mk nc_zero "chr01 15 25"
mk nc_mixed "1 15 25" "chr1 15 25"
mk nc_same "chr1 15 25"

for mode in "" "-sorted"; do
	# chr1 against 1: warning, no overlap because the names differ.
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_bare" $mode
	warned "***** WARNING: File $TMPD/"
	warned "has inconsistent naming convention for record:"
	# -nonamecheck silences it and leaves the (empty) result alone.
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_bare" $mode -nonamecheck
	quiet
	# chr1 against chr01: the leading-zero warning.
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_zero" $mode
	warned "has a record where naming convention (leading zero) is inconsistent with other files:"
	begin_n
	okrun "$TMPD/empty.bed" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_zero" $mode -nonamecheck
	quiet
	# A database mixing 1 and chr1: the chr1 record still overlaps A, with a warning.
	begin_n
	okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_mixed" $mode
	warned "has inconsistent naming convention for record:"
	begin_n
	okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_mixed" $mode -nonamecheck
	quiet
	# Consistent names: no warning.
	begin_n
	okrun "$TMPD/exp_overlap" "$BT" intersect -a "$TMPD/nc_chr" -b "$TMPD/nc_same" $mode
	quiet
done

# ======================================================================
# The intersectBed alias gives the same results as "bedtools intersect".
# ======================================================================

begin_n
okrun "$TMPD/exp_s_int" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb"
begin_n
okrun "$TMPD/exp_s_wawb" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb" -wa -wb
begin_n
okrun "$TMPD/exp_s_c" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb" -c -sorted
begin_n
okrun "$TMPD/exp_s_v" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb" -v
begin_n
okrun "$TMPD/exp_s_loj" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb" -loj -sorted
begin_n
okrun "$TMPD/exp_s_wo" "$IBT" -a "$TMPD/sa" -b "$TMPD/sb" -wo
# Errors and their exit status come through the alias unchanged.
begin_n
argerr "ERROR: request either -u for anyHit OR -v for noHit, not both." "$IBT" -a "$TMPD/a3.bed" -b "$TMPD/b3.bed" -u -v
begin_n
runerr "has the following out of order record" "$IBT" -a "$TMPD/a_unsorted" -b "$TMPD/b_ok" -sorted

# ======================================================================
# Records the parsers reject (message on stderr, exit status 1).
# ======================================================================

# Start greater than end on the first data line: the type check rejects the file.
mk bad_se "chr1 20 10"
begin_n
runerr "Error: unable to open file or unable to determine types for file $TMPD/bad_se" "$BT" intersect -a "$TMPD/bad_se" -b "$TMPD/b3.bed"
begin_n
runerr "Error: unable to open file or unable to determine types for file $TMPD/bad_se" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/bad_se"

# Negative start.
mk bad_neg "chr1 -5 10"
begin_n
runerr "Error: Invalid record in file $TMPD/bad_neg" "$BT" intersect -a "$TMPD/bad_neg" -b "$TMPD/b3.bed"
begin_n
runerr "Error: Invalid record in file $TMPD/bad_neg" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/bad_neg"

# Non-numeric start on the first data line.
mk bad_nn "chr1 x 10"
begin_n
runerr "unable to determine types for file $TMPD/bad_nn" "$BT" intersect -a "$TMPD/bad_nn" -b "$TMPD/b3.bed"

# Too few columns.
mk bad_two "chr1 10"
begin_n
runerr "unable to determine types for file $TMPD/bad_two" "$BT" intersect -a "$TMPD/bad_two" -b "$TMPD/b3.bed"
begin_n
runerr "unable to determine types for file $TMPD/bad_two" "$BT" intersect -a "$TMPD/a3.bed" -b "$TMPD/bad_two"

# Different column counts in the first two lines.
mk bad_cols "chr1 10 20" "chr1 30 40 5"
begin_n
runerr "Error: Type checker found wrong number of fields while tokenizing data line." "$BT" intersect -a "$TMPD/bad_cols" -b "$TMPD/b3.bed"

# Records further into a long file: a bad line 5001 (extra column, bad
# interval, bad number).
awk 'BEGIN { for (i = 0; i < 5000; i++) printf "chr1\t%d\t%d\n", i * 10, i * 10 + 5 }' > "$TMPD/long_ok"
cp "$TMPD/long_ok" "$TMPD/long_cols"
printf 'chr1\t60000\t60005\tx\n' >> "$TMPD/long_cols"
begin_n
runerr "Error: line number 5001 of file $TMPD/long_cols has 4 fields, but 3 were expected." "$BT" intersect -a "$TMPD/long_cols" -b "$TMPD/b3.bed"
begin_n
runerr "Error: line number 5001 of file $TMPD/long_cols has 4 fields, but 3 were expected." "$BT" intersect -a "$TMPD/b3.bed" -b "$TMPD/long_cols"
cp "$TMPD/long_ok" "$TMPD/long_se"
printf 'chr1\t60000\t59999\n' >> "$TMPD/long_se"
begin_n
runerr "Error: Invalid record in file $TMPD/long_se" "$BT" intersect -a "$TMPD/long_se" -b "$TMPD/b3.bed"
cp "$TMPD/long_ok" "$TMPD/long_nn"
printf 'chr1\t60000\tfoo\n' >> "$TMPD/long_nn"
begin_n
runerr "***** ERROR: illegal number \"foo\". Exiting..." "$BT" intersect -a "$TMPD/long_nn" -b "$TMPD/b3.bed"

# The long file without the bad line is read without complaint.  A record
# i is chr1:10i-(10i+5), so only i=1 (10-15) overlaps B chr1:10-20.
begin_n
mk exp_long "chr1 10 15"
okrun "$TMPD/exp_long" "$BT" intersect -a "$TMPD/long_ok" -b "$TMPD/e_one.bed"

finish_all
