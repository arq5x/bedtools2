set -e;
echo -e \
"\n###########################################################
#
#  INTERSECT OPTIONS
#
###########################################################\n"

BT=${BT-../../../bin/bedtools}

FAILURES=0;

check()
{
	if diff $1 $2; then
		echo ok
	else
		FAILURES=$(expr $FAILURES + 1);
		echo fail
	fi
}


###########################################################
###########################################################
#  ARGUMENT VALIDATION
#
#  Usage errors print the help text and then the error on
#  stderr, and exit 1. These tests keep the exit code, stdout,
#  the help's "Tool:" line (line 2, so the Version line is
#  never compared) and the last line, which holds the error.
###########################################################
###########################################################

###########################################################
#  intersect.opt.t001: no arguments at all is an error
#  No files were named, so openFiles() reports that no input
#  was given.
###########################################################
echo -e "    intersect.opt.t001...\c"
$BT intersect > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: No input file given. Exiting. *****\n' > exp
check obs exp
rm obs exp out err

###########################################################
#  intersect.opt.t002: -a without -b is an error
###########################################################
echo -e "    intersect.opt.t002...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
$BT intersect -a a.tmp.bed > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: query and database files not specified. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed

###########################################################
#  intersect.opt.t003: -a as the last argument has no file
###########################################################
echo -e "    intersect.opt.t003...\c"
printf 'chr1\t100\t200\n' > b.tmp.bed
$BT intersect -b b.tmp.bed -a > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -a option given, but no query file specified. *****\n' > exp
check obs exp
rm obs exp out err b.tmp.bed

###########################################################
#  intersect.opt.t004: -b as the last argument has no file
###########################################################
echo -e "    intersect.opt.t004...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
$BT intersect -a a.tmp.bed -b > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -b option given, but no database file specified. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed

###########################################################
#  intersect.opt.t005: -u and -v together are an error
###########################################################
echo -e "    intersect.opt.t005...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -u -v > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -u for anyHit OR -v for noHit, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t006: -c and -u together are an error
###########################################################
echo -e "    intersect.opt.t006...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -c -u > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -c for writeCount OR -u for anyHit, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t007: -c and -wb together are an error
###########################################################
echo -e "    intersect.opt.t007...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -c -wb > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -c for writeCount OR -wb for writeB, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t008: a negative -f is outside (0.0, 1.0]
###########################################################
echo -e "    intersect.opt.t008...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -f -0.5 > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -f must be in the range (0.0, 1.0]. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t009: a non-numeric -f reads as 0 and is
#  rejected by the range check
###########################################################
echo -e "    intersect.opt.t009...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -f half > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -f must be in the range (0.0, 1.0]. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t010: -f as the last argument has no value
#  and is left unconsumed, so it is reported as unrecognized
###########################################################
echo -e "    intersect.opt.t010...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -f > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Unrecognized parameter: -f *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t011: -F 0 is outside (0.0, 1.0]
###########################################################
echo -e "    intersect.opt.t011...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0 > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -F must be in the range (0.0, 1.0]. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t012: -F above 1.0 is outside (0.0, 1.0]
###########################################################
echo -e "    intersect.opt.t012...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 1.5 > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -F must be in the range (0.0, 1.0]. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t013: -e with -F but without -f is an error
###########################################################
echo -e "    intersect.opt.t013...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.5 -e > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -e must be used with -f. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t014: -r with -F (and no -f) is an error
###########################################################
echo -e "    intersect.opt.t014...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.5 -r > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -r must be used solely with -f. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t015: -r with both -f and -F is an error
###########################################################
echo -e "    intersect.opt.t015...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 -r > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -r must be used solely with -f. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t016: -c with -names is an error (use -C)
###########################################################
echo -e "    intersect.opt.t016...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
printf 'chr1\t160\t260\n' > c.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -names B C -c > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -c cannot be used with -filenames or -names. Use -C. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t017: -c with -filenames is an error (use -C)
###########################################################
echo -e "    intersect.opt.t017...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
printf 'chr1\t160\t260\n' > c.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -filenames -c > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -c cannot be used with -filenames or -names. Use -C. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t018: fewer -names than -b files is an error
###########################################################
echo -e "    intersect.opt.t018...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
printf 'chr1\t160\t260\n' > c.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -names B > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Number of database name tags given does not match number of databases. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t019: more -names than -b files is an error,
#  also when there is only one -b file
###########################################################
echo -e "    intersect.opt.t019...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -names B C > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Number of database name tags given does not match number of databases. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t020: -wo with -wa is an error
###########################################################
echo -e "    intersect.opt.t020...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -wa > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -wa for writeA OR -wo for writeOverlap, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t021: -wo with -wb is an error
###########################################################
echo -e "    intersect.opt.t021...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -wb > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -wb for writeB OR -wo for writeOverlap, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t022: -wo with -c is an error
###########################################################
echo -e "    intersect.opt.t022...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -c > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -c for writeCount OR -wo for writeOverlap, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t023: -wao with -u is an error; -wao turns
#  on writeOverlap, so the message names -wo
###########################################################
echo -e "    intersect.opt.t023...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao -u > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request either -u for anyHit OR -wo for writeOverlap, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t024: -s and -S together are an error
###########################################################
echo -e "    intersect.opt.t024...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t250\tb\t0\t+\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -s -S > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t025: -s with a BED5 A file, which has no
#  strand column, is an error naming that file
###########################################################
echo -e "    intersect.opt.t025...\c"
printf 'chr1\t100\t200\ta\t0\n' > a.tmp.bed
printf 'chr1\t150\t250\tb\t0\t+\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -s > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: stranded intersect requested, but input file a.tmp.bed does not have strands. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t026: -S with a bedGraph B file (numeric
#  fourth column, no strand) is an error naming that file
###########################################################
echo -e "    intersect.opt.t026...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t250\t3.5\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -S > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: stranded intersect requested, but input file b.tmp.bed does not have strands. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t027: every unrecognized parameter gets its
#  own error line, in command-line order
###########################################################
echo -e "    intersect.opt.t027...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -wq -b b.tmp.bed -zz > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p' err >> obs
tail -n 2 err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Unrecognized parameter: -wq *****\n***** ERROR: Unrecognized parameter: -zz *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t028: -g as the last argument has no file
###########################################################
echo -e "    intersect.opt.t028...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -g > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -g option given, but no genome file specified. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t029: a -g file that does not exist stops the
#  run while the arguments are read, without the help text
###########################################################
echo -e "    intersect.opt.t029...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -g missing.tmp.genome > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
printf 'Error: Can'"'"'t open genome file missing.tmp.genome. Exiting...\nexit 1\n' > exp
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t030: -iobuf as the last argument has no size
###########################################################
echo -e "    intersect.opt.t030...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: -iobuf option given, but size of input buffer not specified. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t031: -iobuf below the 8-byte minimum is an
#  error
###########################################################
echo -e "    intersect.opt.t031...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 7 > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: specified buffer size is too small. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t032: -iobuf suffixes are upper case K/M/G
#  only; a lower-case k is rejected
###########################################################
echo -e "    intersect.opt.t032...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 64k > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Unrecognized memory buffer size suffix '"'"'k'"'"' given. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t033: a -iobuf value ending in a digit but
#  not numeric is an error
###########################################################
echo -e "    intersect.opt.t033...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 1x6 > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: argument passed to -iobuf is not numeric. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t034: -delim belongs to tools with column
#  operations and is an error for intersect
###########################################################
echo -e "    intersect.opt.t034...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -delim , > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;$p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n***** ERROR: Can'"'"'t set delimiter for tools without column operations. Exiting. *****\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t035: -h prints the help to stderr only and
#  exits 1; line 6 is the usage line
###########################################################
echo -e "    intersect.opt.t035...\c"
$BT intersect -h > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p;6p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\nUsage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf/bam> -b <bed/gff/vcf/bam>\n' > exp
check obs exp
rm obs exp out err

###########################################################
#  intersect.opt.t036: -h wins over otherwise valid arguments:
#  no intersection is run
###########################################################
echo -e "    intersect.opt.t036...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed --help > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat out >> obs
sed -n '2p' err >> obs
printf 'exit 1\nTool:    bedtools intersect (aka intersectBed)\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  REPORTING MODES, ONE -b FILE
#
#  Inputs shared by this section (removed at its end):
#    A (BED6)                  B (BED6)
#    chr1 100 200 a1 1 +       chr1 150 250 b1 5 +
#    chr1 300 400 a2 2 -       chr1 180 190 b2 6 -
#    chr1 500 600 a3 3 +       chr1 390 450 b3 7 -
#    chr2 100 200 a4 4 -       chr3 100 200 b4 8 +
#  a1 overlaps b1 on 150-200 (50 bp) and b2 on 180-190
#  (10 bp); a2 overlaps b3 on 390-400 (10 bp); a3 has no
#  hit; chr2 is only in A and chr3 only in B.
#  b1 and b2 sit in the same 16 kb bin, so they come back in
#  B's file order (pins the current output: hit order).
###########################################################
###########################################################
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\nchr1\t500\t600\ta3\t3\t+\nchr2\t100\t200\ta4\t4\t-\n' > modes_a.tmp.bed
printf 'chr1\t150\t250\tb1\t5\t+\nchr1\t180\t190\tb2\t6\t-\nchr1\t390\t450\tb3\t7\t-\nchr3\t100\t200\tb4\t8\t+\n' > modes_b.tmp.bed

###########################################################
#  intersect.opt.t037: default mode prints the overlapping
#  part of A, one line per hit, with A's other fields
###########################################################
echo -e "    intersect.opt.t037...\c"
printf 'chr1\t150\t200\ta1\t1\t+\nchr1\t180\t190\ta1\t1\t+\nchr1\t390\t400\ta2\t2\t-\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t038: -wa prints the whole A record once per
#  hit (a1 twice)
###########################################################
echo -e "    intersect.opt.t038...\c"
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -wa > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t039: -wb alone prints the overlapping part
#  of A followed by the whole B record
###########################################################
echo -e "    intersect.opt.t039...\c"
printf 'chr1\t150\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t180\t190\ta1\t1\t+\tchr1\t180\t190\tb2\t6\t-\n' >> exp
printf 'chr1\t390\t400\ta2\t2\t-\tchr1\t390\t450\tb3\t7\t-\n' >> exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t040: -wa -wb prints whole A and whole B
###########################################################
echo -e "    intersect.opt.t040...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t180\t190\tb2\t6\t-\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\tb3\t7\t-\n' >> exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t041: -loj prints A and B like -wa -wb and
#  adds a NULL BED6 record for A records without a hit, in A
#  order (a3, and a4 whose chromosome B lacks)
###########################################################
echo -e "    intersect.opt.t041...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t180\t190\tb2\t6\t-\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\tb3\t7\t-\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t.\t-1\t-1\t.\t-1\t.\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t-1\t-1\t.\t-1\t.\n' >> exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -loj > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t042: -wo adds the overlap in bp and skips A
#  records without a hit
###########################################################
echo -e "    intersect.opt.t042...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t180\t190\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\tb3\t7\t-\t10\n' >> exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -wo > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t043: -wao is -wo plus a NULL BED6 record
#  and an overlap of 0 for A records without a hit
###########################################################
echo -e "    intersect.opt.t043...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t180\t190\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\tb3\t7\t-\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t044: -u prints each A record with a hit
#  once, however many hits it has
###########################################################
echo -e "    intersect.opt.t044...\c"
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -u > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t045: -c appends the number of hits to every
#  A record, 0 included
###########################################################
echo -e "    intersect.opt.t045...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t2\nchr1\t300\t400\ta2\t2\t-\t1\nchr1\t500\t600\ta3\t3\t+\t0\nchr2\t100\t200\ta4\t4\t-\t0\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -c > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t046: -C with one -b file prints one line per
#  A record with its count and no file column
###########################################################
echo -e "    intersect.opt.t046...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t2\nchr1\t300\t400\ta2\t2\t-\t1\nchr1\t500\t600\ta3\t3\t+\t0\nchr2\t100\t200\ta4\t4\t-\t0\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -C > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t047: -v prints only the A records without a
#  hit
###########################################################
echo -e "    intersect.opt.t047...\c"
printf 'chr1\t500\t600\ta3\t3\t+\nchr2\t100\t200\ta4\t4\t-\n' > exp
$BT intersect -a modes_a.tmp.bed -b modes_b.tmp.bed -v > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t048: swapping the files swaps the roles: B
#  records are now reported, with the A overlaps
###########################################################
echo -e "    intersect.opt.t048...\c"
printf 'chr1\t150\t250\tb1\t5\t+\tchr1\t100\t200\ta1\t1\t+\t50\n' > exp
printf 'chr1\t180\t190\tb2\t6\t-\tchr1\t100\t200\ta1\t1\t+\t10\n' >> exp
printf 'chr1\t390\t450\tb3\t7\t-\tchr1\t300\t400\ta2\t2\t-\t10\n' >> exp
printf 'chr3\t100\t200\tb4\t8\t+\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
$BT intersect -a modes_b.tmp.bed -b modes_a.tmp.bed -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t049: with a BED3 A file and a BED6 B file,
#  -loj prints each record in its own width
###########################################################
echo -e "    intersect.opt.t049...\c"
printf 'chr1\t100\t200\nchr1\t500\t600\n' > a.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t100\t200\tchr1\t180\t190\tb2\t6\t-\n' >> exp
printf 'chr1\t500\t600\t.\t-1\t-1\t.\t-1\t.\n' >> exp
$BT intersect -a a.tmp.bed -b modes_b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.opt.t050: -loj with a BED3 B file prints a BED3
#  NULL record
###########################################################
echo -e "    intersect.opt.t050...\c"
printf 'chr1\t150\t250\nchr1\t390\t450\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\n' > exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t.\t-1\t-1\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t-1\t-1\n' >> exp
$BT intersect -a modes_a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp b.tmp.bed

rm modes_a.tmp.bed modes_b.tmp.bed

###########################################################
###########################################################
#  REPORTING MODES, SEVERAL -b FILES
#
#  Inputs shared by this section (removed at its end):
#    A (BED6)                B (BED6)               C (BED6)
#    chr1 100 200 a1 1 +     chr1 150 250 b1 5 +    chr1 120 130 c1 7 -
#    chr1 300 400 a2 2 -     chr1 390 450 b2 6 -    chr1 550 560 c2 8 +
#    chr1 500 600 a3 3 +
#    chr2 100 200 a4 4 -
#  a1 hits b1 (50 bp) and c1 (10 bp); a2 hits b2 (10 bp); a3
#  hits c2 (10 bp); a4 has no hit.
#  With "-a A -b B C" the file ids are A=0, B=1, C=2. All four
#  B/C records sit in one 16 kb bin and B is loaded first, so
#  a1's hits come back as b1 then c1 (pins the current output:
#  hit order across -b files).
###########################################################
###########################################################
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\nchr1\t500\t600\ta3\t3\t+\nchr2\t100\t200\ta4\t4\t-\n' > multi_a.tmp.bed
printf 'chr1\t150\t250\tb1\t5\t+\nchr1\t390\t450\tb2\t6\t-\n' > multi_b.tmp.bed
printf 'chr1\t120\t130\tc1\t7\t-\nchr1\t550\t560\tc2\t8\t+\n' > multi_c.tmp.bed

###########################################################
#  intersect.opt.t051: default mode with two -b files prints
#  only the A part, with no file column
###########################################################
echo -e "    intersect.opt.t051...\c"
printf 'chr1\t150\t200\ta1\t1\t+\nchr1\t120\t130\ta1\t1\t+\nchr1\t390\t400\ta2\t2\t-\nchr1\t550\t560\ta3\t3\t+\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t052: -wa with two -b files repeats A once
#  per hit in either file
###########################################################
echo -e "    intersect.opt.t052...\c"
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\nchr1\t500\t600\ta3\t3\t+\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wa > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t053: -wb with two -b files puts the file id
#  between the A part and the B record
###########################################################
echo -e "    intersect.opt.t053...\c"
printf 'chr1\t150\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t120\t130\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\n' >> exp
printf 'chr1\t390\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\n' >> exp
printf 'chr1\t550\t560\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t054: -wa -wb with two -b files puts the file
#  id between the two whole records
###########################################################
echo -e "    intersect.opt.t054...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t055: -wo with two -b files: A, file id, B,
#  overlap
###########################################################
echo -e "    intersect.opt.t055...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\t10\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wo > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t056: -wao with two -b files gives a4 a "."
#  file id, a NULL BED6 record and 0
###########################################################
echo -e "    intersect.opt.t056...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\t10\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t057: -u with two -b files reports A records
#  with a hit in any file, once
###########################################################
echo -e "    intersect.opt.t057...\c"
printf 'chr1\t100\t200\ta1\t1\t+\nchr1\t300\t400\ta2\t2\t-\nchr1\t500\t600\ta3\t3\t+\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -u > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t058: -c with two -b files counts hits in all
#  files together
###########################################################
echo -e "    intersect.opt.t058...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t2\nchr1\t300\t400\ta2\t2\t-\t1\nchr1\t500\t600\ta3\t3\t+\t1\nchr2\t100\t200\ta4\t4\t-\t0\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -c > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t059: -C with two -b files prints one line
#  per A record and file, in file-id order, zeros included
###########################################################
echo -e "    intersect.opt.t059...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\t1\nchr1\t100\t200\ta1\t1\t+\t2\t1\n' > exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\t1\nchr1\t300\t400\ta2\t2\t-\t2\t0\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t1\t0\nchr1\t500\t600\ta3\t3\t+\t2\t1\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t1\t0\nchr2\t100\t200\ta4\t4\t-\t2\t0\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -C > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t060: -v with two -b files reports A records
#  with no hit in either file
###########################################################
echo -e "    intersect.opt.t060...\c"
printf 'chr2\t100\t200\ta4\t4\t-\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -v > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t061: -names replaces the file id in -wo
#  output; the NULL row of -wao keeps "."
###########################################################
echo -e "    intersect.opt.t061...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tbee\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\tsea\tchr1\t120\t130\tc1\t7\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tbee\tchr1\t390\t450\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\tsea\tchr1\t550\t560\tc2\t8\t+\t10\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -names bee sea -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t062: -names replaces the file id in -wb
#  output, after the clipped A part
###########################################################
echo -e "    intersect.opt.t062...\c"
printf 'chr1\t150\t200\ta1\t1\t+\tbee\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t120\t130\ta1\t1\t+\tsea\tchr1\t120\t130\tc1\t7\t-\n' >> exp
printf 'chr1\t390\t400\ta2\t2\t-\tbee\tchr1\t390\t450\tb2\t6\t-\n' >> exp
printf 'chr1\t550\t560\ta3\t3\t+\tsea\tchr1\t550\t560\tc2\t8\t+\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wb -names bee sea > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t063: -filenames prints each -b file name as
#  it was typed (here with a ./ prefix); the -loj NULL row
#  keeps "."
###########################################################
echo -e "    intersect.opt.t063...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tmulti_b.tmp.bed\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t./multi_c.tmp.bed\tchr1\t120\t130\tc1\t7\t-\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\tmulti_b.tmp.bed\tchr1\t390\t450\tb2\t6\t-\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t./multi_c.tmp.bed\tchr1\t550\t560\tc2\t8\t+\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t.\t-1\t-1\t.\t-1\t.\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed ./multi_c.tmp.bed -filenames -loj > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t066: two separate -b options act like one -b
#  followed by both files (same rows as the -wo test above)
###########################################################
echo -e "    intersect.opt.t066...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\t10\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed -b multi_c.tmp.bed -wo > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t067: with one -b file, -names and -filenames
#  add no column
###########################################################
echo -e "    intersect.opt.t067...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tchr1\t150\t250\tb1\t5\t+\n' > exp
printf 'chr1\t300\t400\ta2\t2\t-\tchr1\t390\t450\tb2\t6\t-\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed -names bee -wa -wb > obs
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed -filenames -wa -wb >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t068: the NULL record of -wao takes its shape
#  from the first -b file: BED3 here, though C is BED6
#  (pins the current output: NULL record type)
###########################################################
echo -e "    intersect.opt.t068...\c"
printf 'chr1\t150\t250\nchr1\t390\t450\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\t50\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\t10\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\t10\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t.\t.\t-1\t-1\t0\n' >> exp
$BT intersect -a multi_a.tmp.bed -b b.tmp.bed multi_c.tmp.bed -wao > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t069: -sortout orders each A record's hits by
#  position across files: a1 gets c1 (120) before b1 (150)
###########################################################
echo -e "    intersect.opt.t069...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t2\tchr1\t120\t130\tc1\t7\t-\n' > exp
printf 'chr1\t100\t200\ta1\t1\t+\t1\tchr1\t150\t250\tb1\t5\t+\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\tchr1\t390\t450\tb2\t6\t-\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t2\tchr1\t550\t560\tc2\t8\t+\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -wa -wb -sortout > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t070: -S with two -b files keeps only
#  opposite-strand hits: a1(+) keeps c1(-), a3(+) loses c2(+),
#  a2(-) loses b2(-)
###########################################################
echo -e "    intersect.opt.t070...\c"
printf 'chr1\t100\t200\ta1\t1\t+\tsea\tchr1\t120\t130\tc1\t7\t-\t10\n' > exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -names bee sea -S -wo > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t071: -C counts per file after -s and after
#  -f are applied
#    -s: a1(+) keeps b1(+) not c1(-); a2(-) keeps b2(-); a3(+)
#        keeps c2(+)
#    -f 0.5: only b1 (50 of a1's 100 bp) passes
###########################################################
echo -e "    intersect.opt.t071...\c"
printf 'chr1\t100\t200\ta1\t1\t+\t1\t1\nchr1\t100\t200\ta1\t1\t+\t2\t0\n' > exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\t1\nchr1\t300\t400\ta2\t2\t-\t2\t0\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t1\t0\nchr1\t500\t600\ta3\t3\t+\t2\t1\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t1\t0\nchr2\t100\t200\ta4\t4\t-\t2\t0\n' >> exp
printf 'chr1\t100\t200\ta1\t1\t+\t1\t1\nchr1\t100\t200\ta1\t1\t+\t2\t0\n' >> exp
printf 'chr1\t300\t400\ta2\t2\t-\t1\t0\nchr1\t300\t400\ta2\t2\t-\t2\t0\n' >> exp
printf 'chr1\t500\t600\ta3\t3\t+\t1\t0\nchr1\t500\t600\ta3\t3\t+\t2\t0\n' >> exp
printf 'chr2\t100\t200\ta4\t4\t-\t1\t0\nchr2\t100\t200\ta4\t4\t-\t2\t0\n' >> exp
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -C -s > obs
$BT intersect -a multi_a.tmp.bed -b multi_b.tmp.bed multi_c.tmp.bed -C -f 0.5 >> obs
check obs exp
rm obs exp

rm multi_a.tmp.bed multi_b.tmp.bed multi_c.tmp.bed

###########################################################
###########################################################
#  HIT ORDER AND -sortout
#
#  Unsorted intersect loads B into a bin index (8 levels,
#  16 kb bins at the finest level, 8x larger per level) and
#  reports an A record's hits level by level from the finest,
#  bin by bin within a level, and in load order within a bin.
#  None of this is documented, so the tests below say
#  "pins the current output: hit order".
###########################################################
###########################################################

###########################################################
#  intersect.opt.t072: an A record with hits at several bin
#  levels (pins the current output: hit order)
#  A is chr1:0-200000. B, in file order:
#    r1 16000-17000    crosses 16384, so a 128 kb (level 1) bin
#    r2 100-200        16 kb bin 0
#    r3 20000-20100    16 kb bin 1
#    r4 5000-6000      16 kb bin 0
#    r5 130000-140000  crosses 131072, so a 1 Mb (level 2) bin
#  Order: bin 0 (r2, r4 in file order), bin 1 (r3), the
#  level-1 bin (r1), the level-2 bin (r5).
###########################################################
echo -e "    intersect.opt.t072...\c"
printf 'chr1\t0\t200000\tA\n' > a.tmp.bed
printf 'chr1\t16000\t17000\tr1\nchr1\t100\t200\tr2\nchr1\t20000\t20100\tr3\nchr1\t5000\t6000\tr4\nchr1\t130000\t140000\tr5\n' > b.tmp.bed
printf 'chr1\t100\t200\tA\tchr1\t100\t200\tr2\n' > exp
printf 'chr1\t5000\t6000\tA\tchr1\t5000\t6000\tr4\n' >> exp
printf 'chr1\t20000\t20100\tA\tchr1\t20000\t20100\tr3\n' >> exp
printf 'chr1\t16000\t17000\tA\tchr1\t16000\t17000\tr1\n' >> exp
printf 'chr1\t130000\t140000\tA\tchr1\t130000\t140000\tr5\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t073: -sortout with one -b file orders the
#  same five hits by start: r2, r4, r1, r3, r5
###########################################################
echo -e "    intersect.opt.t073...\c"
printf 'chr1\t0\t200000\tA\n' > a.tmp.bed
printf 'chr1\t16000\t17000\tr1\nchr1\t100\t200\tr2\nchr1\t20000\t20100\tr3\nchr1\t5000\t6000\tr4\nchr1\t130000\t140000\tr5\n' > b.tmp.bed
printf 'chr1\t100\t200\tA\tchr1\t100\t200\tr2\n' > exp
printf 'chr1\t5000\t6000\tA\tchr1\t5000\t6000\tr4\n' >> exp
printf 'chr1\t16000\t17000\tA\tchr1\t16000\t17000\tr1\n' >> exp
printf 'chr1\t20000\t20100\tA\tchr1\t20000\t20100\tr3\n' >> exp
printf 'chr1\t130000\t140000\tA\tchr1\t130000\t140000\tr5\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb -sortout > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t074: -sorted reports the hits in B's (sorted)
#  file order, the same order -sortout gives above
###########################################################
echo -e "    intersect.opt.t074...\c"
printf 'chr1\t0\t200000\tA\n' > a.tmp.bed
printf 'chr1\t100\t200\tr2\nchr1\t5000\t6000\tr4\nchr1\t16000\t17000\tr1\nchr1\t20000\t20100\tr3\nchr1\t130000\t140000\tr5\n' > b.tmp.bed
printf 'chr1\t100\t200\tA\tchr1\t100\t200\tr2\n' > exp
printf 'chr1\t5000\t6000\tA\tchr1\t5000\t6000\tr4\n' >> exp
printf 'chr1\t16000\t17000\tA\tchr1\t16000\t17000\tr1\n' >> exp
printf 'chr1\t20000\t20100\tA\tchr1\t20000\t20100\tr3\n' >> exp
printf 'chr1\t130000\t140000\tA\tchr1\t130000\t140000\tr5\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t075: hits in one 16 kb bin keep B's file
#  order even when that is descending (pins the current
#  output: hit order)
###########################################################
echo -e "    intersect.opt.t075...\c"
printf 'chr1\t0\t1000\tA\n' > a.tmp.bed
printf 'chr1\t300\t310\tx\nchr1\t200\t210\ty\nchr1\t100\t110\tz\n' > b.tmp.bed
printf 'chr1\t0\t1000\tA\tchr1\t300\t310\tx\n' > exp
printf 'chr1\t0\t1000\tA\tchr1\t200\t210\ty\n' >> exp
printf 'chr1\t0\t1000\tA\tchr1\t100\t110\tz\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t076: hits in different 16 kb bins come back
#  by bin, so ascending, though B lists them descending
#  (pins the current output: hit order)
###########################################################
echo -e "    intersect.opt.t076...\c"
printf 'chr1\t0\t100000\tA\n' > a.tmp.bed
printf 'chr1\t80000\t80010\tp\nchr1\t60000\t60010\tq\nchr1\t40000\t40010\tr\nchr1\t20000\t20010\ts\nchr1\t0\t10\tt\n' > b.tmp.bed
printf 'chr1\t0\t10\tA\nchr1\t20000\t20010\tA\nchr1\t40000\t40010\tA\nchr1\t60000\t60010\tA\nchr1\t80000\t80010\tA\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t077: -sortout with -wo sorts the descending
#  same-bin hits of the test before last into z, y, x
###########################################################
echo -e "    intersect.opt.t077...\c"
printf 'chr1\t0\t1000\tA\n' > a.tmp.bed
printf 'chr1\t300\t310\tx\nchr1\t200\t210\ty\nchr1\t100\t110\tz\n' > b.tmp.bed
printf 'chr1\t0\t1000\tA\tchr1\t100\t110\tz\t10\n' > exp
printf 'chr1\t0\t1000\tA\tchr1\t200\t210\ty\t10\n' >> exp
printf 'chr1\t0\t1000\tA\tchr1\t300\t310\tx\t10\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -sortout > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t078: -sortout sorts by start, then end
#  B: 100-300, 100-200, 50-400 -> 50-400, 100-200, 100-300
###########################################################
echo -e "    intersect.opt.t078...\c"
printf 'chr1\t0\t1000\tA\n' > a.tmp.bed
printf 'chr1\t100\t300\tu\nchr1\t100\t200\tv\nchr1\t50\t400\tw\n' > b.tmp.bed
printf 'chr1\t50\t400\tA\nchr1\t100\t200\tA\nchr1\t100\t300\tA\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -sortout > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t079: -sortout leaves hits with equal
#  coordinates in their load order (pins the current output:
#  order of equal hits). This holds only for small hit sets:
#  -sortout uses std::sort (RecordKeyVector::sortVector),
#  which keeps tied hits in load order only while the A
#  record's hits fall under the library's insertion-sort
#  threshold (16 hits in libstdc++; 16 tied hits plus one
#  other are already reordered), and libstdc++ (Linux) and
#  libc++ (macOS) differ above it. Tests with tied hits under
#  -sortout keep to 16 or fewer hits per A record and never
#  assert the order of 17 or more.
###########################################################
echo -e "    intersect.opt.t079...\c"
printf 'chr1\t0\t1000\tA\n' > a.tmp.bed
printf 'chr1\t100\t200\tfirst\nchr1\t100\t200\tsecond\nchr1\t50\t60\tearly\n' > b.tmp.bed
printf 'chr1\t0\t1000\tA\tchr1\t50\t60\tearly\n' > exp
printf 'chr1\t0\t1000\tA\tchr1\t100\t200\tfirst\n' >> exp
printf 'chr1\t0\t1000\tA\tchr1\t100\t200\tsecond\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb -sortout > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t080: bin order outranks -b order: the second
#  file's record in a 16 kb bin comes before the first file's
#  record in a 128 kb bin (pins the current output: hit order)
###########################################################
echo -e "    intersect.opt.t080...\c"
printf 'chr1\t0\t200000\tA\n' > a.tmp.bed
printf 'chr1\t16000\t17000\tbig\n' > b.tmp.bed
printf 'chr1\t100\t200\tsmall\n' > c.tmp.bed
printf 'chr1\t0\t200000\tA\t2\tchr1\t100\t200\tsmall\n' > exp
printf 'chr1\t0\t200000\tA\t1\tchr1\t16000\t17000\tbig\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t081: within one bin, hits from three -b
#  files come back in -b order, whatever their positions
#  (pins the current output: hit order)
###########################################################
echo -e "    intersect.opt.t081...\c"
printf 'chr1\t0\t1000\tA\n' > a.tmp.bed
printf 'chr1\t300\t310\tb\n' > b.tmp.bed
printf 'chr1\t200\t210\tc\n' > c.tmp.bed
printf 'chr1\t100\t110\td\n' > d.tmp.bed
printf 'chr1\t0\t1000\tA\t1\tchr1\t300\t310\tb\t10\n' > exp
printf 'chr1\t0\t1000\tA\t2\tchr1\t200\t210\tc\t10\n' >> exp
printf 'chr1\t0\t1000\tA\t3\tchr1\t100\t110\td\t10\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed d.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed d.tmp.bed

###########################################################
###########################################################
#  STRAND: -s AND -S
#
#  Inputs shared by this section (removed at its end):
#    A (BED6)                B (BED6), by A record
#    chr1 100 200 p 0 +      150-160 bp +, 170-180 bm -
#    chr1 300 400 m 0 -      350-360 bp2 +, 370-380 bm2 -
#    chr1 500 600 d 0 .      550-560 bd ., 570-580 bp3 +
#    chr1 700 800 s 0 *      750-760 bs *, 770-780 bm3 -
#    chr1 900 1000 q 0 +     950-960 bq +
#  Only + and - are strands; '.' and '*' are unknown and can
#  neither match nor differ. Same strand: p-bp, m-bm2, q-bq.
#  Opposite strand: p-bm, m-bp2.
###########################################################
###########################################################
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t300\t400\tm\t0\t-\nchr1\t500\t600\td\t0\t.\nchr1\t700\t800\ts\t0\t*\nchr1\t900\t1000\tq\t0\t+\n' > str_a.tmp.bed
printf 'chr1\t150\t160\tbp\t0\t+\nchr1\t170\t180\tbm\t0\t-\nchr1\t350\t360\tbp2\t0\t+\nchr1\t370\t380\tbm2\t0\t-\n' > str_b.tmp.bed
printf 'chr1\t550\t560\tbd\t0\t.\nchr1\t570\t580\tbp3\t0\t+\nchr1\t750\t760\tbs\t0\t*\nchr1\t770\t780\tbm3\t0\t-\nchr1\t950\t960\tbq\t0\t+\n' >> str_b.tmp.bed

###########################################################
#  intersect.opt.t082: -s keeps same-strand hits only; '.' and
#  '*' records get none
###########################################################
echo -e "    intersect.opt.t082...\c"
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t150\t160\tbp\t0\t+\n' > exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t370\t380\tbm2\t0\t-\n' >> exp
printf 'chr1\t900\t1000\tq\t0\t+\tchr1\t950\t960\tbq\t0\t+\n' >> exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t083: -S keeps opposite-strand hits only;
#  '.' and '*' are not opposite to anything
###########################################################
echo -e "    intersect.opt.t083...\c"
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t170\t180\tbm\t0\t-\n' > exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t350\t360\tbp2\t0\t+\n' >> exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -S -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t084: -s in default mode prints the clipped A
#  part of each same-strand hit
###########################################################
echo -e "    intersect.opt.t084...\c"
printf 'chr1\t150\t160\tp\t0\t+\nchr1\t370\t380\tm\t0\t-\nchr1\t950\t960\tq\t0\t+\n' > exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t085: -s -v reports the A records with no
#  same-strand hit: the '.' and '*' records
###########################################################
echo -e "    intersect.opt.t085...\c"
printf 'chr1\t500\t600\td\t0\t.\nchr1\t700\t800\ts\t0\t*\n' > exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s -v > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t086: -S -v also reports q, whose only hit is
#  on its own strand
###########################################################
echo -e "    intersect.opt.t086...\c"
printf 'chr1\t500\t600\td\t0\t.\nchr1\t700\t800\ts\t0\t*\nchr1\t900\t1000\tq\t0\t+\n' > exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -S -v > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t087: -s -c counts same-strand hits only
###########################################################
echo -e "    intersect.opt.t087...\c"
printf 'chr1\t100\t200\tp\t0\t+\t1\nchr1\t300\t400\tm\t0\t-\t1\nchr1\t500\t600\td\t0\t.\t0\nchr1\t700\t800\ts\t0\t*\t0\nchr1\t900\t1000\tq\t0\t+\t1\n' > exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s -c > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t088: -S -u reports p and m once each
###########################################################
echo -e "    intersect.opt.t088...\c"
printf 'chr1\t100\t200\tp\t0\t+\nchr1\t300\t400\tm\t0\t-\n' > exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -S -u > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t089: -s -wao gives the NULL row to A records
#  whose hits are all filtered out by strand
###########################################################
echo -e "    intersect.opt.t089...\c"
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t150\t160\tbp\t0\t+\t10\n' > exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t370\t380\tbm2\t0\t-\t10\n' >> exp
printf 'chr1\t500\t600\td\t0\t.\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
printf 'chr1\t700\t800\ts\t0\t*\t.\t-1\t-1\t.\t-1\t.\t0\n' >> exp
printf 'chr1\t900\t1000\tq\t0\t+\tchr1\t950\t960\tbq\t0\t+\t10\n' >> exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t090: -S -loj gives the NULL record to d, s
#  and q
###########################################################
echo -e "    intersect.opt.t090...\c"
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t170\t180\tbm\t0\t-\n' > exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t350\t360\tbp2\t0\t+\n' >> exp
printf 'chr1\t500\t600\td\t0\t.\t.\t-1\t-1\t.\t-1\t.\n' >> exp
printf 'chr1\t700\t800\ts\t0\t*\t.\t-1\t-1\t.\t-1\t.\n' >> exp
printf 'chr1\t900\t1000\tq\t0\t+\t.\t-1\t-1\t.\t-1\t.\n' >> exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -S -loj > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t091: -s and -S give the same rows with
#  -sorted (both inputs are sorted)
###########################################################
echo -e "    intersect.opt.t091...\c"
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t150\t160\tbp\t0\t+\t10\n' > exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t370\t380\tbm2\t0\t-\t10\n' >> exp
printf 'chr1\t900\t1000\tq\t0\t+\tchr1\t950\t960\tbq\t0\t+\t10\n' >> exp
printf 'chr1\t100\t200\tp\t0\t+\tchr1\t170\t180\tbm\t0\t-\t10\n' >> exp
printf 'chr1\t300\t400\tm\t0\t-\tchr1\t350\t360\tbp2\t0\t+\t10\n' >> exp
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -s -wo -sorted > obs
$BT intersect -a str_a.tmp.bed -b str_b.tmp.bed -S -wo -sorted >> obs
check obs exp
rm obs exp

rm str_a.tmp.bed str_b.tmp.bed

###########################################################
#  intersect.opt.t092: -s reads the strand from column 6 of a
#  BED6+ record and keeps its extra column
###########################################################
echo -e "    intersect.opt.t092...\c"
printf 'chr1\t100\t200\tx\t0\t+\textra1\nchr1\t300\t400\ty\t0\t-\textra2\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t0\t+\nchr1\t350\t360\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\tx\t0\t+\textra1\tchr1\t150\t160\tb1\t0\t+\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -s -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t093: -S reads the strand from column 6 of a
#  BED12 record
###########################################################
echo -e "    intersect.opt.t093...\c"
printf 'chr1\t100\t500\ttx\t0\t-\t100\t500\t0\t2\t100,100,\t0,300,\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t0\t+\nchr1\t450\t460\tb2\t0\t-\n' > b.tmp.bed
printf 'chr1\t100\t500\ttx\t0\t-\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t160\tb1\t0\t+\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -S -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t094: -s is applied before -f: A(+) is fully
#  covered by an opposite-strand B, and 10% by a same-strand B
#    -f 0.5       keeps the opposite-strand B (100%)
#    -s -f 0.5    keeps nothing
#    -s -f 0.1    keeps the same-strand B (10/100 = 0.1)
###########################################################
echo -e "    intersect.opt.t094...\c"
printf 'chr1\t100\t200\ta\t0\t+\n' > a.tmp.bed
printf 'chr1\t100\t200\tminus\t0\t-\nchr1\t190\t200\tplus\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\t0\t+\tchr1\t100\t200\tminus\t0\t-\n' > exp
printf 'chr1\t100\t200\ta\t0\t+\tchr1\t190\t200\tplus\t0\t+\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -s -f 0.5 -wa -wb >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -s -f 0.1 -wa -wb >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  OVERLAP FRACTIONS: -f, -F, -r, -e
#
#  A hit passes -f when overlap / len(A) >= f, and -F when
#  overlap / len(B) >= F (both in single precision). With -e
#  either test suffices; -r sets F = f.
###########################################################
###########################################################

###########################################################
#  intersect.opt.t095: -f 0.5 keeps an overlap of exactly half
#  of A; -f 0.51 drops it
#  A 100-200 (100 bp), B 150-300: 50 bp = 0.50 of A.
###########################################################
echo -e "    intersect.opt.t095...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t150\t300\tb\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -wa > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.51 -wa >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t096: -f 1.0 keeps an A record covered end to
#  end and drops one covered 99 of 100 bp
#  in: 100-200 inside B 50-250; part: 300-400 vs B 301-500.
###########################################################
echo -e "    intersect.opt.t096...\c"
printf 'chr1\t100\t200\tin\nchr1\t300\t400\tpart\n' > a.tmp.bed
printf 'chr1\t50\t250\tb1\nchr1\t301\t500\tb2\n' > b.tmp.bed
printf 'chr1\t100\t200\tin\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 1.0 -u > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t097: -f at exact decimal boundaries
#  Three 10 bp A records overlap B by 1, 3 and 7 bp.
#    -f 0.1 keeps all three, -f 0.3 keeps 3 and 7 bp,
#    -f 0.7 keeps 7 bp: the fraction equals f exactly.
###########################################################
echo -e "    intersect.opt.t097...\c"
printf 'chr1\t0\t10\tone\nchr1\t100\t110\tthree\nchr1\t200\t210\tseven\n' > a.tmp.bed
printf 'chr1\t9\t20\nchr1\t107\t120\nchr1\t203\t220\n' > b.tmp.bed
printf 'chr1\t0\t10\tone\nchr1\t100\t110\tthree\nchr1\t200\t210\tseven\n' > exp
printf 'chr1\t100\t110\tthree\nchr1\t200\t210\tseven\n' >> exp
printf 'chr1\t200\t210\tseven\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.1 -u > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.3 -u >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.7 -u >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t098: -f just above those boundaries drops
#  them: -f 0.11 keeps 3 and 7 bp, -f 0.31 keeps 7 bp,
#  -f 0.71 keeps none
###########################################################
echo -e "    intersect.opt.t098...\c"
printf 'chr1\t0\t10\tone\nchr1\t100\t110\tthree\nchr1\t200\t210\tseven\n' > a.tmp.bed
printf 'chr1\t9\t20\nchr1\t107\t120\nchr1\t203\t220\n' > b.tmp.bed
printf 'chr1\t100\t110\tthree\nchr1\t200\t210\tseven\n' > exp
printf 'chr1\t200\t210\tseven\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.11 -u > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.31 -u >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.71 -u >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t099: a 1 bp overlap of a 3 bp A is 1/3, so
#  -f 0.333333 keeps it and -f 0.334 drops it
###########################################################
echo -e "    intersect.opt.t099...\c"
printf 'chr1\t10\t13\ta\n' > a.tmp.bed
printf 'chr1\t12\t50\tb\n' > b.tmp.bed
printf 'chr1\t12\t13\ta\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.333333 > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.334 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t100: -F 0.5 keeps an overlap of exactly half
#  of B; -F 0.51 drops it
#  A 100-200, B 150-250 (100 bp): 50 bp = 0.50 of B.
###########################################################
echo -e "    intersect.opt.t100...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t150\t250\tb\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t150\t250\tb\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.5 -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.51 -wa -wb >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t101: -F 1.0 keeps a B inside A and drops a B
#  reaching 1 bp past A's end
###########################################################
echo -e "    intersect.opt.t101...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t120\t200\tin\nchr1\t150\t201\tout\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t120\t200\tin\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 1.0 -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t102: -f and -F together must both pass
#  A 100-200. big 100-300: A 1.0, B 0.5. small 150-160: A
#  0.1, B 1.0.
#    -f 0.5 -F 0.5 keeps big; -f 0.5 -F 0.6 keeps nothing.
###########################################################
echo -e "    intersect.opt.t102...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t100\t300\tbig\nchr1\t150\t160\tsmall\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t100\t300\tbig\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.5 -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.6 -wa -wb >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t103: -e lets either fraction pass: with
#  -f 0.5 -F 0.6 -e, big passes on A (1.0) and small on B
#  (1.0); same inputs as the test above
###########################################################
echo -e "    intersect.opt.t103...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t100\t300\tbig\nchr1\t150\t160\tsmall\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t100\t300\tbig\n' > exp
printf 'chr1\t100\t200\ta\tchr1\t150\t160\tsmall\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -F 0.6 -e -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t104: -e with -f and no -F keeps any overlap:
#  -F is unset, so the B test passes for any hit. A 1 bp hit
#  fails -f 0.9 alone and passes -f 0.9 -e
###########################################################
echo -e "    intersect.opt.t104...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t199\t300\tb\n' > b.tmp.bed
printf 'chr1\t199\t200\ta\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.9 > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.9 -e >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t105: -r applies -f to B as well
#  A 100-200 (100 bp). b1 100-300: A 1.0, B 100/200 = 0.5.
#  b2 101-300: A 0.99, B 99/199 < 0.5.
#    -f 0.5 keeps both; -f 0.5 -r keeps b1 only.
###########################################################
echo -e "    intersect.opt.t105...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t100\t300\tb1\nchr1\t101\t300\tb2\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t100\t300\tb1\n' > exp
printf 'chr1\t100\t200\ta\tchr1\t101\t300\tb2\n' >> exp
printf 'chr1\t100\t200\ta\tchr1\t100\t300\tb1\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 -r -wa -wb >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t106: -r without -f changes nothing: a 1 bp
#  hit is still reported
###########################################################
echo -e "    intersect.opt.t106...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t199\t300\tb\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t199\t300\tb\t1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -r -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t107: -f 1E-9, the documented default, gives
#  the default rows: a 1 bp hit kept, a touching B not
###########################################################
echo -e "    intersect.opt.t107...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t199\t300\tone\nchr1\t200\t300\ttouch\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t199\t300\tone\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 1E-9 -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  Inputs for the next six tests (removed after them):
#    A: x 100-200, y 300-400
#    B: bx1 190-260 (10 bp of x = 0.1), bx2 100-160 (60 bp of
#       x = 0.6), by 399-500 (1 bp of y = 0.01)
#  With -f 0.5 only bx2 passes; y keeps no hit.
###########################################################
printf 'chr1\t100\t200\tx\nchr1\t300\t400\ty\n' > frac_a.tmp.bed
printf 'chr1\t190\t260\tbx1\nchr1\t100\t160\tbx2\nchr1\t399\t500\tby\n' > frac_b.tmp.bed

###########################################################
#  intersect.opt.t108: -f filters what -c counts
###########################################################
echo -e "    intersect.opt.t108...\c"
printf 'chr1\t100\t200\tx\t1\nchr1\t300\t400\ty\t0\n' > exp
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -f 0.5 -c > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t109: with -f, -u keeps x and -v reports y,
#  whose only hit is too small
###########################################################
echo -e "    intersect.opt.t109...\c"
printf 'chr1\t100\t200\tx\nchr1\t300\t400\ty\n' > exp
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -f 0.5 -u > obs
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -f 0.5 -v >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t110: with -f, -wao shows only passing hits
#  and gives y the NULL row
###########################################################
echo -e "    intersect.opt.t110...\c"
printf 'chr1\t100\t200\tx\tchr1\t100\t160\tbx2\t60\n' > exp
printf 'chr1\t300\t400\ty\t.\t-1\t-1\t.\t0\n' >> exp
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -f 0.5 -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t111: -F filters -loj rows: bx2 covers all of
#  its 60 bp, bx1 only 10 of 70 and by 1 of 101, so y gets
#  the NULL record
###########################################################
echo -e "    intersect.opt.t111...\c"
printf 'chr1\t100\t200\tx\tchr1\t100\t160\tbx2\n' > exp
printf 'chr1\t300\t400\ty\t.\t-1\t-1\t.\n' >> exp
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -F 0.5 -loj > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t112: -f in default mode prints only the
#  passing hit's part of x
###########################################################
echo -e "    intersect.opt.t112...\c"
printf 'chr1\t100\t160\tx\n' > exp
$BT intersect -a frac_a.tmp.bed -b frac_b.tmp.bed -f 0.5 > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t113: -f gives the same rows with -sorted
#  (-wao test above, B sorted first)
###########################################################
echo -e "    intersect.opt.t113...\c"
printf 'chr1\t100\t160\tbx2\nchr1\t190\t260\tbx1\nchr1\t399\t500\tby\n' > b.tmp.bed
printf 'chr1\t100\t200\tx\tchr1\t100\t160\tbx2\t60\n' > exp
printf 'chr1\t300\t400\ty\t.\t-1\t-1\t.\t0\n' >> exp
$BT intersect -a frac_a.tmp.bed -b b.tmp.bed -f 0.5 -wao -sorted > obs
check obs exp
rm obs exp b.tmp.bed

rm frac_a.tmp.bed frac_b.tmp.bed

###########################################################
#  intersect.opt.t114: fractions above 2^31: A 3000000000-
#  3000000100, B 3000000050-3000000200 overlap 50 bp
#  (A 0.5, B 1/3). -f 0.5 keeps, -f 0.6 drops, -F 0.3 keeps,
#  -F 0.34 drops.
###########################################################
echo -e "    intersect.opt.t114...\c"
printf 'chr1\t3000000000\t3000000100\ta\n' > a.tmp.bed
printf 'chr1\t3000000050\t3000000200\tb\n' > b.tmp.bed
printf 'chr1\t3000000050\t3000000100\ta\n' > exp
printf 'chr1\t3000000050\t3000000100\ta\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.5 > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 0.6 >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.3 >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 0.34 >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t115: a zero-length A inside B passes
#  -f 1.0 (pins the current output: a zero-length record is
#  widened by 1 bp on each side, so 2 of 2 bp overlap)
###########################################################
echo -e "    intersect.opt.t115...\c"
printf 'chr1\t150\t150\tz\n' > a.tmp.bed
printf 'chr1\t100\t200\tb\n' > b.tmp.bed
printf 'chr1\t150\t150\tz\tchr1\t100\t200\tb\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 1.0 -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t116: a zero-length B inside A passes
#  -F 1.0 (pins the current output: zero-length widening)
###########################################################
echo -e "    intersect.opt.t116...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t150\t150\tz\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t150\t150\tz\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -F 1.0 -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  -split WITH BED12
#
#  Inputs shared by this section (removed at its end):
#    tx (BED12): chr1 100-500, strand +, blocks 100-200 and
#                400-500 (sizes 100,100  starts 0,300)
#    q (BED6):   q1 150-180 +  inside block 1
#                q2 250-350 +  inside the intron
#                q3 420-460 -  inside block 2
#                q4 0-1000  +  covers everything
#  With -split, q2 overlaps nothing; q4 overlaps both blocks
#  (200 bp). All records share one 16 kb bin, so hits come
#  back in file order (q1, q2, q3, q4), not by position (q4
#  starts first): pins the current output: hit order within
#  one bin follows file order, not position.
###########################################################
###########################################################
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\n' > split_tx.tmp.bed
printf 'chr1\t150\t180\tq1\t0\t+\nchr1\t250\t350\tq2\t0\t+\nchr1\t420\t460\tq3\t0\t-\nchr1\t0\t1000\tq4\t0\t+\n' > split_q.tmp.bed

###########################################################
#  intersect.opt.t117: -split drops the B record that lies in
#  A's intron (pins the current output: hit order within one
#  bin follows file order, not position)
###########################################################
echo -e "    intersect.opt.t117...\c"
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t180\tq1\t0\t+\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t420\t460\tq3\t0\t-\n' >> exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t0\t1000\tq4\t0\t+\n' >> exp
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -split -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t118: -c counts 4 hits without -split and 3
#  with it
###########################################################
echo -e "    intersect.opt.t118...\c"
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\t4\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\t3\n' >> exp
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -c > obs
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -split -c >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t119: -split -wo reports block overlap only:
#  q1 30, q3 40, q4 100 + 100 = 200 (pins the current output:
#  hit order within one bin follows file order, not position)
###########################################################
echo -e "    intersect.opt.t119...\c"
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t180\tq1\t0\t+\t30\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t420\t460\tq3\t0\t-\t40\n' >> exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t0\t1000\tq4\t0\t+\t200\n' >> exp
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -split -wo > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t120: -split with -s keeps q1 and q4 (+);
#  with -S it keeps q3 (-) (pins the current output: hit order
#  within one bin follows file order, not position)
###########################################################
echo -e "    intersect.opt.t120...\c"
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t180\tq1\t0\t+\t30\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t0\t1000\tq4\t0\t+\t200\n' >> exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t420\t460\tq3\t0\t-\t40\n' >> exp
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -split -s -wo > obs
$BT intersect -a split_tx.tmp.bed -b split_q.tmp.bed -split -S -wo >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t121: default mode with -split prints the A
#  part clipped to each B lying inside one block, keeping the
#  BED12 fields as they are
###########################################################
echo -e "    intersect.opt.t121...\c"
printf 'chr1\t150\t180\tq1\t0\t+\nchr1\t420\t460\tq3\t0\t-\n' > b.tmp.bed
printf 'chr1\t150\t180\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\n' > exp
printf 'chr1\t420\t460\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\n' >> exp
$BT intersect -a split_tx.tmp.bed -b b.tmp.bed -split > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t122: with BED12 as B, -split -v reports the
#  A record in the intron and -split -u the other three
###########################################################
echo -e "    intersect.opt.t122...\c"
printf 'chr1\t250\t350\tq2\t0\t+\n' > exp
printf 'chr1\t150\t180\tq1\t0\t+\nchr1\t420\t460\tq3\t0\t-\nchr1\t0\t1000\tq4\t0\t+\n' >> exp
$BT intersect -a split_q.tmp.bed -b split_tx.tmp.bed -split -v > obs
$BT intersect -a split_q.tmp.bed -b split_tx.tmp.bed -split -u >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t123: with BED12 as B, -split -wao gives q2 a
#  NULL BED12 record (12 fields) and 0
###########################################################
echo -e "    intersect.opt.t123...\c"
printf 'chr1\t150\t180\tq1\t0\t+\tchr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\t30\n' > exp
printf 'chr1\t250\t350\tq2\t0\t+\t.\t-1\t-1\t.\t-1\t.\t.\t.\t.\t.\t.\t.\t0\n' >> exp
printf 'chr1\t420\t460\tq3\t0\t-\tchr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\t40\n' >> exp
printf 'chr1\t0\t1000\tq4\t0\t+\tchr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\t200\n' >> exp
$BT intersect -a split_q.tmp.bed -b split_tx.tmp.bed -split -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t124: -split -f / -F measure blocks
#    -f 0.5: q1 30/30 and q3 40/40 pass, q4 200/1000 fails
#    -F 0.6: q4 covers all 200 of tx's block bp (1.0) and
#      passes; by tx's 400 bp extent it would be 0.5 and fail
###########################################################
echo -e "    intersect.opt.t124...\c"
printf 'chr1\t150\t180\tq1\t0\t+\nchr1\t420\t460\tq3\t0\t-\n' > exp
printf 'chr1\t0\t1000\tq4\t0\t+\n' >> exp
$BT intersect -a split_q.tmp.bed -b split_tx.tmp.bed -split -f 0.5 -u > obs
$BT intersect -a split_q.tmp.bed -b split_tx.tmp.bed -split -F 0.6 -u >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t125: BED12 on both sides
#    tx2 150-550, blocks 150-250 and 450-550: with -split it
#      overlaps tx on 150-200 and 450-500 = 100 bp; without,
#      on 150-500 = 350 bp
#    tx3 210-390, blocks 210-260 and 340-390: inside tx's
#      intron, so a hit (180 bp) only without -split
###########################################################
echo -e "    intersect.opt.t125...\c"
printf 'chr1\t150\t550\ttx2\t0\t+\t150\t550\t0\t2\t100,100,\t0,300,\nchr1\t210\t390\ttx3\t0\t+\t210\t390\t0\t2\t50,50,\t0,130,\n' > b.tmp.bed
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t550\ttx2\t0\t+\t150\t550\t0\t2\t100,100,\t0,300,\t350\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t210\t390\ttx3\t0\t+\t210\t390\t0\t2\t50,50,\t0,130,\t180\n' >> exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t550\ttx2\t0\t+\t150\t550\t0\t2\t100,100,\t0,300,\t100\n' >> exp
$BT intersect -a split_tx.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a split_tx.tmp.bed -b b.tmp.bed -split -wo >> obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t126: -split with BED6 on both sides changes
#  nothing (q against a BED6 100-500)
###########################################################
echo -e "    intersect.opt.t126...\c"
printf 'chr1\t100\t500\tr\t0\t+\n' > b.tmp.bed
printf 'chr1\t150\t180\tq1\t0\t+\tchr1\t100\t500\tr\t0\t+\t30\n' > exp
printf 'chr1\t250\t350\tq2\t0\t+\tchr1\t100\t500\tr\t0\t+\t100\n' >> exp
printf 'chr1\t420\t460\tq3\t0\t-\tchr1\t100\t500\tr\t0\t+\t40\n' >> exp
printf 'chr1\t0\t1000\tq4\t0\t+\tchr1\t100\t500\tr\t0\t+\t400\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a split_q.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a split_q.tmp.bed -b b.tmp.bed -split -wo >> obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t127: -split gives the same rows with and
#  without -sorted (q sorted by start: q4, q1, q2, q3)
###########################################################
echo -e "    intersect.opt.t127...\c"
printf 'chr1\t0\t1000\tq4\t0\t+\nchr1\t150\t180\tq1\t0\t+\nchr1\t250\t350\tq2\t0\t+\nchr1\t420\t460\tq3\t0\t-\n' > b.tmp.bed
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t0\t1000\tq4\t0\t+\t200\n' > exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t150\t180\tq1\t0\t+\t30\n' >> exp
printf 'chr1\t100\t500\ttx\t0\t+\t100\t500\t0\t2\t100,100,\t0,300,\tchr1\t420\t460\tq3\t0\t-\t40\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a split_tx.tmp.bed -b b.tmp.bed -split -wo > obs
$BT intersect -a split_tx.tmp.bed -b b.tmp.bed -split -wo -sorted >> obs
check obs exp
rm obs exp b.tmp.bed

rm split_tx.tmp.bed split_q.tmp.bed

###########################################################
###########################################################
#  -sorted (CHROMSWEEP) AND -g
#
#  Inputs shared by the first six tests (removed after them),
#  both sorted:
#    A                    B
#    chr1 100 200 a1      chr1 150 350 b1    a1 50 bp, a2 50 bp
#    chr1 300 400 a2      chr1 380 390 b2    a2 10 bp
#    chr1 500 600 a3      chr15 10 20 b5     chromosome only in B
#    chr2 100 200 a4      chr2 150 250 b3    a4 50 bp
#    chr3 50 60 a5        chr2 300 400 b4    no hit
#  Chromosomes are in lexicographic order (chr1 < chr15 <
#  chr2 < chr3). Each test runs unsorted and then -sorted and
#  expects the same rows twice.
###########################################################
###########################################################
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\nchr1\t500\t600\ta3\nchr2\t100\t200\ta4\nchr3\t50\t60\ta5\n' > srt_a.tmp.bed
printf 'chr1\t150\t350\tb1\nchr1\t380\t390\tb2\nchr15\t10\t20\tb5\nchr2\t150\t250\tb3\nchr2\t300\t400\tb4\n' > srt_b.tmp.bed

###########################################################
#  intersect.opt.t128: default mode, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t128...\c"
printf 'chr1\t150\t200\ta1\nchr1\t300\t350\ta2\nchr1\t380\t390\ta2\nchr2\t150\t200\ta4\n' > exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t129: -wa -wb, unsorted and -sorted; b1 is
#  reported for both a1 and a2
###########################################################
echo -e "    intersect.opt.t129...\c"
printf 'chr1\t100\t200\ta1\tchr1\t150\t350\tb1\n' > exp
printf 'chr1\t300\t400\ta2\tchr1\t150\t350\tb1\n' >> exp
printf 'chr1\t300\t400\ta2\tchr1\t380\t390\tb2\n' >> exp
printf 'chr2\t100\t200\ta4\tchr2\t150\t250\tb3\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -wa -wb > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -wa -wb -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t130: -wao, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t130...\c"
printf 'chr1\t100\t200\ta1\tchr1\t150\t350\tb1\t50\n' > exp
printf 'chr1\t300\t400\ta2\tchr1\t150\t350\tb1\t50\n' >> exp
printf 'chr1\t300\t400\ta2\tchr1\t380\t390\tb2\t10\n' >> exp
printf 'chr1\t500\t600\ta3\t.\t-1\t-1\t.\t0\n' >> exp
printf 'chr2\t100\t200\ta4\tchr2\t150\t250\tb3\t50\n' >> exp
printf 'chr3\t50\t60\ta5\t.\t-1\t-1\t.\t0\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -wao > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -wao -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t131: -loj, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t131...\c"
printf 'chr1\t100\t200\ta1\tchr1\t150\t350\tb1\n' > exp
printf 'chr1\t300\t400\ta2\tchr1\t150\t350\tb1\n' >> exp
printf 'chr1\t300\t400\ta2\tchr1\t380\t390\tb2\n' >> exp
printf 'chr1\t500\t600\ta3\t.\t-1\t-1\t.\n' >> exp
printf 'chr2\t100\t200\ta4\tchr2\t150\t250\tb3\n' >> exp
printf 'chr3\t50\t60\ta5\t.\t-1\t-1\t.\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -loj > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -loj -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t132: -u then -v, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t132...\c"
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\nchr2\t100\t200\ta4\n' > exp
printf 'chr1\t500\t600\ta3\nchr3\t50\t60\ta5\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -u > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -v >> obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -u -sorted >> obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -v -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t133: -c, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t133...\c"
printf 'chr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t2\nchr1\t500\t600\ta3\t0\nchr2\t100\t200\ta4\t1\nchr3\t50\t60\ta5\t0\n' > exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -c > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -c -sorted >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t134: -C, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t134...\c"
printf 'chr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t2\nchr1\t500\t600\ta3\t0\nchr2\t100\t200\ta4\t1\nchr3\t50\t60\ta5\t0\n' > exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -C > obs
$BT intersect -a srt_a.tmp.bed -b srt_b.tmp.bed -C -sorted >> obs
check obs exp
rm obs exp

rm srt_a.tmp.bed srt_b.tmp.bed

###########################################################
#  intersect.opt.t135: -sorted with a chromosome only A has
#  between two shared ones (lexicographic chr1, chr15, chr2):
#  -v reports its record and the chr2 hit is still found
###########################################################
echo -e "    intersect.opt.t135...\c"
printf 'chr1\t10\t20\ta1\nchr15\t10\t20\ta15\nchr2\t10\t20\ta2\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr2\t15\t25\tb2\n' > b.tmp.bed
printf 'chr15\t10\t20\ta15\n' > exp
printf 'chr1\t10\t20\ta1\tchr1\t15\t25\tb1\t5\nchr15\t10\t20\ta15\t.\t-1\t-1\t.\t0\nchr2\t10\t20\ta2\tchr2\t15\t25\tb2\t5\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -v -sorted > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t136: -sorted keeps a long B record for every
#  later A record it overlaps
###########################################################
echo -e "    intersect.opt.t136...\c"
printf 'chr1\t100\t110\tx\nchr1\t200\t210\ty\nchr1\t300\t310\tz\nchr1\t2000\t2010\tw\n' > a.tmp.bed
printf 'chr1\t0\t1000\tlong\n' > b.tmp.bed
printf 'chr1\t100\t110\tx\tchr1\t0\t1000\tlong\n' > exp
printf 'chr1\t200\t210\ty\tchr1\t0\t1000\tlong\n' >> exp
printf 'chr1\t300\t310\tz\tchr1\t0\t1000\tlong\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t137: -sorted finds a hit again for a nested
#  A record after a long A record has passed it
#  big 0-1000 hits all three B; small 100-150 hits only
#  120-130 (50-60 ends before it, 500-600 starts after).
###########################################################
echo -e "    intersect.opt.t137...\c"
printf 'chr1\t0\t1000\tbig\nchr1\t100\t150\tsmall\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\nchr1\t120\t130\tb2\nchr1\t500\t600\tb3\n' > b.tmp.bed
printf 'chr1\t0\t1000\tbig\tchr1\t50\t60\tb1\n' > exp
printf 'chr1\t0\t1000\tbig\tchr1\t120\t130\tb2\n' >> exp
printf 'chr1\t0\t1000\tbig\tchr1\t500\t600\tb3\n' >> exp
printf 'chr1\t100\t150\tsmall\tchr1\t120\t130\tb2\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t138: -sorted -c and -v keep going after B
#  runs out: A records past the last B record, on the same
#  and on a later chromosome, are still reported
###########################################################
echo -e "    intersect.opt.t138...\c"
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\nchr1\t500\t600\ta3\nchr2\t1\t5\ta4\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t0\nchr1\t500\t600\ta3\t0\nchr2\t1\t5\ta4\t0\n' > exp
printf 'chr1\t300\t400\ta2\nchr1\t500\t600\ta3\nchr2\t1\t5\ta4\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c -sorted > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -v -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t139: -sorted -wao keeps going after B runs
#  out (inputs as in the test above)
###########################################################
echo -e "    intersect.opt.t139...\c"
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\nchr1\t500\t600\ta3\nchr2\t1\t5\ta4\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\tchr1\t150\t160\tb1\t10\n' > exp
printf 'chr1\t300\t400\ta2\t.\t-1\t-1\t.\t0\n' >> exp
printf 'chr1\t500\t600\ta3\t.\t-1\t-1\t.\t0\n' >> exp
printf 'chr2\t1\t5\ta4\t.\t-1\t-1\t.\t0\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t140: without -g, -sorted accepts any
#  chromosome order the files share: chr1, chr2, chr10 in
#  both, then chr1, chr10, chr2 in both
###########################################################
echo -e "    intersect.opt.t140...\c"
printf 'chr1\t10\t20\ta1\nchr2\t10\t20\ta2\nchr10\t10\t20\ta10\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr2\t15\t25\tb2\nchr10\t15\t25\tb10\n' > b.tmp.bed
printf 'chr1\t10\t20\ta1\nchr10\t10\t20\ta10\nchr2\t10\t20\ta2\n' > c.tmp.bed
printf 'chr1\t15\t25\tb1\nchr10\t15\t25\tb10\nchr2\t15\t25\tb2\n' > d.tmp.bed
printf 'chr1\t15\t20\ta1\nchr2\t15\t20\ta2\nchr10\t15\t20\ta10\n' > exp
printf 'chr1\t15\t20\ta1\nchr10\t15\t20\ta10\nchr2\t15\t20\ta2\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > obs
$BT intersect -a c.tmp.bed -b d.tmp.bed -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed d.tmp.bed

###########################################################
#  intersect.opt.t141: -sorted -g follows the genome file's
#  order (chr2 before chr1), with a chromosome only B has
#  (chr3) between them
###########################################################
echo -e "    intersect.opt.t141...\c"
printf 'chr2\t10\t20\ta2\nchr1\t10\t20\ta1\n' > a.tmp.bed
printf 'chr2\t15\t25\tb2\nchr3\t1\t5\tb3\nchr1\t15\t25\tb1\n' > b.tmp.bed
printf 'chr2\t100\nchr3\t100\nchr1\t100\n' > g.tmp.genome
printf 'chr2\t10\t20\ta2\tchr2\t15\t25\tb2\t5\nchr1\t10\t20\ta1\tchr1\t15\t25\tb1\t5\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.genome -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.genome

###########################################################
#  intersect.opt.t142: -sorted -g stops on a chromosome the
#  genome file lacks
###########################################################
echo -e "    intersect.opt.t142...\c"
printf 'chr1\t10\t20\ta1\nchr2\t10\t20\ta2\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\n' > b.tmp.bed
printf 'chr1\t100\n' > g.tmp.genome
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.genome > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: requested chromosome chr2 does not exist in the genome file g.tmp.genome. Exiting.\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed g.tmp.genome

###########################################################
#  intersect.opt.t143: -sorted -g stops when a file's order
#  disagrees with the genome file (files chr1, chr10, chr2;
#  genome chr1, chr2, chr10). B is read ahead, so B's chr2
#  record is the one named (pins the current output: which
#  file and record are named)
###########################################################
echo -e "    intersect.opt.t143...\c"
printf 'chr1\t10\t20\ta1\nchr10\t10\t20\ta10\nchr2\t10\t20\ta2\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr10\t15\t25\tb10\nchr2\t15\t25\tb2\n' > b.tmp.bed
printf 'chr1\t100\nchr2\t100\nchr10\t100\n' > g.tmp.genome
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -g g.tmp.genome > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: Sorted input specified, but the file b.tmp.bed has the following record with a different sort order than the genomeFile g.tmp.genome\nchr2\t15\t25\tb2\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed g.tmp.genome

###########################################################
#  intersect.opt.t144: without -g, -sorted stops when the two
#  files order chromosomes differently (A lexicographic, B
#  natural); pins the current output: which file and record
#  are named
#  The error depends on ContextBase::_runToQueryEnd and
#  ContextIntersect::_shouldRunToDbEnd being false, which no
#  constructor guarantees; with MALLOC_PERTURB_ set, a g++ build
#  exits 0 without the error.
###########################################################
echo -e "    intersect.opt.t144...\c"
printf 'chr1\t10\t20\ta1\nchr10\t10\t20\ta10\nchr2\t10\t20\ta2\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr2\t15\t25\tb2\nchr10\t15\t25\tb10\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nERROR: chromosome sort ordering for file b.tmp.bed is inconsistent with other files. Record was:\nchr10\t15\t25\tb10\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t145: without -g, a B chromosome missing from
#  A is an error once the order is known not to be
#  lexicographic (chr2 before chr1); rows printed before the
#  error are not compared
###########################################################
echo -e "    intersect.opt.t145...\c"
printf 'chr2\t10\t20\ta2\nchr1\t10\t20\ta1\n' > a.tmp.bed
printf 'chr2\t15\t25\tb2\nchr3\t1\t5\tb3\nchr1\t15\t25\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nERROR: Database file b.tmp.bed contains chromosome chr3, but the query file does not.\n       Please rerun with the -g option for a genome file.\n       See documentation for details.\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t146: -sorted stops on an A record that
#  starts before the previous one (stdout is not compared:
#  buffered rows are lost on exit)
###########################################################
echo -e "    intersect.opt.t146...\c"
printf 'chr1\t10\t20\ta1\nchr1\t5\t8\ta0\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t5\t8\ta0\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t147: -sorted stops on an out-of-order B
#  record
###########################################################
echo -e "    intersect.opt.t147...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t120\t130\tb2\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: Sorted input specified, but the file b.tmp.bed has the following out of order record\nchr1\t120\t130\tb2\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t148: -sorted stops when A returns to a
#  chromosome it already left (chr1, chr2, chr1)
###########################################################
echo -e "    intersect.opt.t148...\c"
printf 'chr1\t10\t20\ta1\nchr2\t5\t8\ta2\nchr1\t30\t40\ta3\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr1\t35\t38\tb2\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t30\t40\ta3\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t149: -sorted reads the rest of B after A
#  ends and still catches B's out-of-order record there
###########################################################
echo -e "    intersect.opt.t149...\c"
printf 'chr1\t10\t20\ta1\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr2\t1\t5\tb2\nchr2\t0\t3\tb3\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > out 2> err && echo "exit 0" > obs || echo "exit $?" > obs
cat err >> obs
printf 'exit 1\nError: Sorted input specified, but the file b.tmp.bed has the following out of order record\nchr2\t0\t3\tb3\n' > exp
check obs exp
rm obs exp out err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t150: -sorted with two -b files reports all
#  hits from the first file, then the second, while unsorted
#  follows bin order (pins the current output: hit order)
#  A 0-20000; B big 16000-17000 (a 128 kb bin), C small
#  100-200 (a 16 kb bin).
###########################################################
echo -e "    intersect.opt.t150...\c"
printf 'chr1\t0\t20000\tA\n' > a.tmp.bed
printf 'chr1\t16000\t17000\tbig\n' > b.tmp.bed
printf 'chr1\t100\t200\tsmall\n' > c.tmp.bed
printf 'chr1\t0\t20000\tA\t1\tchr1\t16000\t17000\tbig\n' > exp
printf 'chr1\t0\t20000\tA\t2\tchr1\t100\t200\tsmall\n' >> exp
printf 'chr1\t0\t20000\tA\t2\tchr1\t100\t200\tsmall\n' >> exp
printf 'chr1\t0\t20000\tA\t1\tchr1\t16000\t17000\tbig\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -wa -wb -sorted > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -wa -wb >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t151: -sorted -sortout orders hits from both
#  files by position (inputs as in the test above)
###########################################################
echo -e "    intersect.opt.t151...\c"
printf 'chr1\t0\t20000\tA\n' > a.tmp.bed
printf 'chr1\t16000\t17000\tbig\n' > b.tmp.bed
printf 'chr1\t100\t200\tsmall\n' > c.tmp.bed
printf 'chr1\t0\t20000\tA\t2\tchr1\t100\t200\tsmall\n' > exp
printf 'chr1\t0\t20000\tA\t1\tchr1\t16000\t17000\tbig\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -wa -wb -sorted -sortout > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.opt.t152: -g without -sorted changes nothing when
#  the genome file has every chromosome (its order disagrees
#  with the unsorted files)
###########################################################
echo -e "    intersect.opt.t152...\c"
printf 'chr2\t10\t20\ta2\nchr1\t10\t20\ta1\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\nchr2\t15\t25\tb2\n' > b.tmp.bed
printf 'chr1\t100\nchr2\t100\n' > g.tmp.genome
printf 'chr2\t15\t20\ta2\nchr1\t15\t20\ta1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -g g.tmp.genome > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed g.tmp.genome

###########################################################
###########################################################
#  CHROMOSOME NAMING CHECK AND -nonamecheck
#
#  The first record whose name disagrees with the records
#  seen so far ("chr" prefix or not, a leading zero after
#  "chr" or not) triggers one warning. The warning is
#  printed when found and again at exit (pins the current
#  output: printed twice). Unsorted intersect loads B before
#  reading A, so B sets the convention and the warning names
#  the A file (pins the current output: file named).
###########################################################
###########################################################

###########################################################
#  intersect.opt.t153: unsorted, "chr1" in A against "1" in B
#  warns about the first A record only; rows are unaffected
###########################################################
echo -e "    intersect.opt.t153...\c"
printf 'chr1\t100\t200\ta1\nchr2\t100\t200\ta2\n' > a.tmp.bed
printf '1\t150\t160\tb1\n2\t150\t160\tb2\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'chr1\t100\t200\ta1\t.\t-1\t-1\t.\t0\nchr2\t100\t200\ta2\t.\t-1\t-1\t.\t0\nexit 0\n' > exp
printf '***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\ta1\n\n' >> exp
printf '***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t100\t200\ta1\n\n' >> exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t154: -nonamecheck silences that warning and
#  leaves the rows alone
###########################################################
echo -e "    intersect.opt.t154...\c"
printf 'chr1\t100\t200\ta1\nchr2\t100\t200\ta2\n' > a.tmp.bed
printf '1\t150\t160\tb1\n2\t150\t160\tb2\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao -nonamecheck > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'chr1\t100\t200\ta1\t.\t-1\t-1\t.\t0\nchr2\t100\t200\ta2\t.\t-1\t-1\t.\t0\nexit 0\n' > exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t155: unsorted, "chr1" in A against "chr01"
#  in B gives the leading-zero warning
###########################################################
echo -e "    intersect.opt.t155...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'chr01\t150\t160\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'exit 0\n' > exp
printf '***** WARNING: File a.tmp.bed has a record where naming convention (leading zero) is inconsistent with other files:\nchr1\t100\t200\ta1\n\n' >> exp
printf '***** WARNING: File a.tmp.bed has a record where naming convention (leading zero) is inconsistent with other files:\nchr1\t100\t200\ta1\n\n' >> exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t156: -nonamecheck silences the leading-zero
#  warning
###########################################################
echo -e "    intersect.opt.t156...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'chr01\t150\t160\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -nonamecheck > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'exit 0\n' > exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t157: with -sorted the A record is tested
#  first, so the warning names the B file and record (pins
#  the current output: file named)
###########################################################
echo -e "    intersect.opt.t157...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf '1\t150\t160\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'exit 0\n' > exp
printf '***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t150\t160\tb1\n\n' >> exp
printf '***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t150\t160\tb1\n\n' >> exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t158: -sorted -nonamecheck prints no warning
###########################################################
echo -e "    intersect.opt.t158...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf '1\t150\t160\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -nonamecheck -c > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'chr1\t100\t200\ta1\t0\nexit 0\n' > exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t159: names that agree give no warning,
#  whether both lack "chr" or differ only in its case
#  ("chr1" and "Chr1" are still different chromosomes)
###########################################################
echo -e "    intersect.opt.t159...\c"
printf '1\t100\t200\ta1\n' > a.tmp.bed
printf '1\t150\t160\tb1\n' > b.tmp.bed
printf 'chr1\t100\t200\tc1\n' > c.tmp.bed
printf 'Chr1\t150\t160\td1\n' > d.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs 2> err
$BT intersect -a c.tmp.bed -b d.tmp.bed -c >> obs 2>> err
cat err >> obs
printf '1\t100\t200\ta1\t1\nchr1\t100\t200\tc1\t0\n' > exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed c.tmp.bed d.tmp.bed

###########################################################
###########################################################
#  -header
#
#  Inputs shared by this section (removed at its end):
#    A: two "#" header lines (the second has a tab), then
#       chr1 100 200 a1 and chr1 300 400 a2
#    B: one "#" header line, then chr1 150 160 b1
#  a1 overlaps b1 by 10 bp; a2 has no hit.
###########################################################
###########################################################
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\nchr1\t300\t400\ta2\n' > hdr_a.tmp.bed
printf '#bh\nchr1\t150\t160\tb1\n' > hdr_b.tmp.bed

###########################################################
#  intersect.opt.t160: without -header no header is printed
###########################################################
echo -e "    intersect.opt.t160...\c"
printf 'chr1\t150\t160\ta1\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t161: -header prints A's header, never B's,
#  before the rows
###########################################################
echo -e "    intersect.opt.t161...\c"
printf '#h1\n#h2\tx\nchr1\t150\t160\ta1\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t162: -header with -c
###########################################################
echo -e "    intersect.opt.t162...\c"
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t0\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -c > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t163: -header with -v
###########################################################
echo -e "    intersect.opt.t163...\c"
printf '#h1\n#h2\tx\nchr1\t300\t400\ta2\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -v > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t164: -header with -wao
###########################################################
echo -e "    intersect.opt.t164...\c"
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\tchr1\t150\t160\tb1\t10\nchr1\t300\t400\ta2\t.\t-1\t-1\t.\t0\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -wao > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t165: -header with -u and with -C
###########################################################
echo -e "    intersect.opt.t165...\c"
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\n' > exp
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t0\n' >> exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -u > obs
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -C >> obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t166: -header prints the header even when no
#  row follows, unsorted and -sorted
###########################################################
echo -e "    intersect.opt.t166...\c"
printf 'chr9\t1\t2\n' > b.tmp.bed
printf '#h1\n#h2\tx\n#h1\n#h2\tx\n' > exp
$BT intersect -a hdr_a.tmp.bed -b b.tmp.bed -header > obs
$BT intersect -a hdr_a.tmp.bed -b b.tmp.bed -header -sorted >> obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t167: -header -sorted with an empty B file:
#  the sweep ends before any A record is handled, and the
#  header is still printed
###########################################################
echo -e "    intersect.opt.t167...\c"
printf '' > b.tmp.bed
printf '#h1\n#h2\tx\n' > exp
$BT intersect -a hdr_a.tmp.bed -b b.tmp.bed -header -sorted > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.opt.t168: -header with -sorted
###########################################################
echo -e "    intersect.opt.t168...\c"
printf '#h1\n#h2\tx\nchr1\t100\t200\ta1\tchr1\t150\t160\tb1\n' > exp
$BT intersect -a hdr_a.tmp.bed -b hdr_b.tmp.bed -header -sorted -wa -wb > obs
check obs exp
rm obs exp

###########################################################
#  intersect.opt.t169: -header with an A file that has no
#  header adds nothing, even when B has one
###########################################################
echo -e "    intersect.opt.t169...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'chr1\t150\t160\ta1\n' > exp
$BT intersect -a a.tmp.bed -b hdr_b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed

rm hdr_a.tmp.bed hdr_b.tmp.bed

###########################################################
###########################################################
#  -nobuf AND -iobuf
###########################################################
###########################################################

###########################################################
#  intersect.opt.t170: -nobuf gives the same rows as buffered
#  output
###########################################################
echo -e "    intersect.opt.t170...\c"
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\nchr1\t500\t600\ta3\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t550\t700\tb2\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\tchr1\t150\t160\tb1\t10\n' > exp
printf 'chr1\t300\t400\ta2\t.\t-1\t-1\t.\t0\n' >> exp
printf 'chr1\t500\t600\ta3\tchr1\t550\t700\tb2\t50\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao -nobuf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t171: -nobuf writes each row as it is made,
#  so the row for a1 reaches stdout before -sorted stops on
#  the out-of-order record that follows it
###########################################################
echo -e "    intersect.opt.t171...\c"
printf 'chr1\t10\t20\ta1\nchr1\t5\t8\ta0\n' > a.tmp.bed
printf 'chr1\t15\t25\tb1\n' > b.tmp.bed
$BT intersect -a a.tmp.bed -b b.tmp.bed -sorted -nobuf > obs 2> err && echo "exit 0" >> obs || echo "exit $?" >> obs
cat err >> obs
printf 'chr1\t15\t20\ta1\nexit 1\nError: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t5\t8\ta0\n' > exp
check obs exp
rm obs exp err a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t172: -iobuf 8, the smallest buffer allowed,
#  reads lines longer than the buffer correctly in both files
###########################################################
echo -e "    intersect.opt.t172...\c"
printf 'chr1\t100\t200\tfirst_record\t0\t+\nchr1\t300\t400\tsecond_record\t0\t-\nchr2\t1000\t2000\tthird\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t350\tspanning_b\t9\t+\nchr2\t1999\t3000\tlast_base\t9\t-\n' > b.tmp.bed
printf 'chr1\t100\t200\tfirst_record\t0\t+\tchr1\t150\t350\tspanning_b\t9\t+\t50\n' > exp
printf 'chr1\t300\t400\tsecond_record\t0\t-\tchr1\t150\t350\tspanning_b\t9\t+\t50\n' >> exp
printf 'chr2\t1000\t2000\tthird\t0\t+\tchr2\t1999\t3000\tlast_base\t9\t-\t1\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -iobuf 8 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t173: -iobuf 8 with -sorted (inputs as above)
###########################################################
echo -e "    intersect.opt.t173...\c"
printf 'chr1\t100\t200\tfirst_record\t0\t+\nchr1\t300\t400\tsecond_record\t0\t-\nchr2\t1000\t2000\tthird\t0\t+\n' > a.tmp.bed
printf 'chr1\t150\t350\tspanning_b\t9\t+\nchr2\t1999\t3000\tlast_base\t9\t-\n' > b.tmp.bed
printf 'chr1\t100\t200\tfirst_record\t0\t+\tchr1\t150\t350\tspanning_b\t9\t+\t50\n' > exp
printf 'chr1\t300\t400\tsecond_record\t0\t-\tchr1\t150\t350\tspanning_b\t9\t+\t50\n' >> exp
printf 'chr2\t1000\t2000\tthird\t0\t+\tchr2\t1999\t3000\tlast_base\t9\t-\t1\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -iobuf 8 -sorted > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t174: -iobuf takes K, M and G suffixes and
#  the rows do not change
###########################################################
echo -e "    intersect.opt.t174...\c"
printf 'chr1\t100\t200\ta1\nchr1\t300\t400\ta2\n' > a.tmp.bed
printf 'chr1\t150\t350\tb1\n' > b.tmp.bed
printf 'chr1\t150\t200\ta1\nchr1\t300\t350\ta2\n' > exp
cat exp exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 1K > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 1M >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -iobuf 1G >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t175: output larger than the 16 KB output
#  buffer comes out whole and in order, buffered, with
#  -nobuf, with -iobuf 8 and with -sorted
#  A = B = 2000 disjoint records chr1:100i-(100i+50), named
#  r<i>, written by awk. Each record hits only itself, so
#  -wa -wb prints every line next to itself (about 75 KB).
###########################################################
echo -e "    intersect.opt.t175...\c"
awk 'BEGIN { for (i = 0; i < 2000; i++) printf "chr1\t%d\t%d\tr%d\n", i * 100, i * 100 + 50, i }' > a.tmp.bed
awk 'BEGIN { for (i = 0; i < 2000; i++) printf "chr1\t%d\t%d\tr%d\tchr1\t%d\t%d\tr%d\n", i * 100, i * 100 + 50, i, i * 100, i * 100 + 50, i }' > exp
cat exp exp exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b a.tmp.bed -wa -wb > obs
$BT intersect -a a.tmp.bed -b a.tmp.bed -wa -wb -nobuf >> obs
$BT intersect -a a.tmp.bed -b a.tmp.bed -wa -wb -iobuf 8 >> obs
$BT intersect -a a.tmp.bed -b a.tmp.bed -wa -wb -sorted >> obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.opt.t176: -header with -iobuf 8 and -nobuf: the
#  header lines, longer than the buffer, come out whole and
#  first
###########################################################
echo -e "    intersect.opt.t176...\c"
printf '#first header line\n#second header line\nchr1\t100\t200\ta1\nchr1\t300\t400\ta2\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\n' > b.tmp.bed
printf '#first header line\n#second header line\nchr1\t100\t200\ta1\t1\nchr1\t300\t400\ta2\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header -c -iobuf 8 -nobuf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  INTERVAL GEOMETRY
#
#  BED intervals are half-open: chr1 100 200 covers bases
#  100..199. A zero-length record (start == end) is widened
#  by 1 bp on each side for the overlap test and printed as
#  given (pins the current output: zero-length handling).
###########################################################
###########################################################

###########################################################
#  intersect.opt.t177: intervals that only touch do not
#  overlap, on either side; -wao gives A the NULL row
###########################################################
echo -e "    intersect.opt.t177...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t0\t100\tleft\nchr1\t200\t300\tright\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\t.\t-1\t-1\t.\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wao > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t178: a 1 bp overlap counts, at either end
###########################################################
echo -e "    intersect.opt.t178...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t0\t101\tleft\nchr1\t199\t300\tright\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t0\t101\tleft\t1\nchr1\t100\t200\ta\tchr1\t199\t300\tright\t1\n' > exp
printf 'chr1\t100\t101\ta\nchr1\t199\t200\ta\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t179: nested intervals: default mode prints
#  the inner interval, -wo its length
#  outer 100-500 holds B 200-300; inner 1200-1300 lies in B
#  1000-2000.
###########################################################
echo -e "    intersect.opt.t179...\c"
printf 'chr1\t100\t500\touter\nchr1\t1200\t1300\tinner\n' > a.tmp.bed
printf 'chr1\t200\t300\tb_in\nchr1\t1000\t2000\tb_out\n' > b.tmp.bed
printf 'chr1\t200\t300\touter\nchr1\t1200\t1300\tinner\n' > exp
printf 'chr1\t100\t500\touter\tchr1\t200\t300\tb_in\t100\nchr1\t1200\t1300\tinner\tchr1\t1000\t2000\tb_out\t100\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t180: identical intervals overlap over their
#  whole length and pass -f 1.0 -r
###########################################################
echo -e "    intersect.opt.t180...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t100\t200\tb\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t100\t200\tb\t100\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -f 1.0 -r -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t181: a B record listed twice is two hits:
#  two -wa -wb rows, a count of 2, one -u row
###########################################################
echo -e "    intersect.opt.t181...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t150\t250\tdup\nchr1\t150\t250\tdup\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t150\t250\tdup\nchr1\t100\t200\ta\tchr1\t150\t250\tdup\n' > exp
printf 'chr1\t100\t200\ta\t2\n' >> exp
printf 'chr1\t100\t200\ta\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -c >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -u >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t182: an A record listed twice is reported
#  twice, by -u when it has a hit and by -v when not
###########################################################
echo -e "    intersect.opt.t182...\c"
printf 'chr1\t100\t200\tx\nchr1\t100\t200\tx\nchr1\t500\t600\ty\nchr1\t500\t600\ty\n' > a.tmp.bed
printf 'chr1\t150\t160\tb\n' > b.tmp.bed
printf 'chr1\t100\t200\tx\nchr1\t100\t200\tx\nchr1\t500\t600\ty\nchr1\t500\t600\ty\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -u > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -v >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t183: duplicated records give the same rows
#  with -sorted
###########################################################
echo -e "    intersect.opt.t183...\c"
printf 'chr1\t100\t200\tx\nchr1\t100\t200\tx\n' > a.tmp.bed
printf 'chr1\t150\t250\tdup\nchr1\t150\t250\tdup\n' > b.tmp.bed
printf 'chr1\t100\t200\tx\tchr1\t150\t250\tdup\t50\n' > exp
cat exp exp exp exp exp exp exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t184: coordinates at and above 2^31 and 2^32
#  and at 10^10, unsorted and -sorted
#    a1 2147483600-2147483700: b1 2147483647-2147483648 (1 bp
#       ending at 2^31) and b2 2147483648-2147483650 (2 bp)
#    a2 4294967290-4294967300: b3 4294967295-4294967297 spans
#       2^32 (2 bp)
#    a3 10000000000-10000000100: b4 10000000050-10000000200
#       (50 bp)
###########################################################
echo -e "    intersect.opt.t184...\c"
printf 'chr1\t2147483600\t2147483700\ta1\nchr1\t4294967290\t4294967300\ta2\nchr1\t10000000000\t10000000100\ta3\n' > a.tmp.bed
printf 'chr1\t2147483647\t2147483648\tb1\nchr1\t2147483648\t2147483650\tb2\nchr1\t4294967295\t4294967297\tb3\nchr1\t10000000050\t10000000200\tb4\n' > b.tmp.bed
printf 'chr1\t2147483600\t2147483700\ta1\tchr1\t2147483647\t2147483648\tb1\t1\n' > exp
printf 'chr1\t2147483600\t2147483700\ta1\tchr1\t2147483648\t2147483650\tb2\t2\n' >> exp
printf 'chr1\t4294967290\t4294967300\ta2\tchr1\t4294967295\t4294967297\tb3\t2\n' >> exp
printf 'chr1\t10000000000\t10000000100\ta3\tchr1\t10000000050\t10000000200\tb4\t50\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t185: default mode clips at large positions:
#  the overlap of a1 with b1 is 2147483647-2147483648, and
#  touching at 2^31 is still no overlap
###########################################################
echo -e "    intersect.opt.t185...\c"
printf 'chr1\t2147483600\t2147483648\ta1\nchr1\t2147483648\t2147483700\ta2\n' > a.tmp.bed
printf 'chr1\t2147483647\t2147483648\tb1\n' > b.tmp.bed
printf 'chr1\t2147483647\t2147483648\ta1\n' > exp
printf 'chr1\t2147483648\t2147483700\ta2\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -v >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t186: -wo prints an overlap of 2^31 - 1 bp
###########################################################
echo -e "    intersect.opt.t186...\c"
printf 'chr1\t0\t2147483647\ta\n' > a.tmp.bed
printf 'chr1\t0\t3000000000\tb\n' > b.tmp.bed
printf 'chr1\t0\t2147483647\ta\tchr1\t0\t3000000000\tb\t2147483647\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t187: -c counts an A record far above 2^31
#  hit by B records in several bins
###########################################################
echo -e "    intersect.opt.t187...\c"
printf 'chr1\t5000000000\t5001000000\ta\n' > a.tmp.bed
printf 'chr1\t5000000000\t5000000001\tb1\nchr1\t5000500000\t5000600000\tb2\nchr1\t5000999999\t5001000001\tb3\nchr1\t5001000000\t5001000010\tb4\n' > b.tmp.bed
printf 'chr1\t5000000000\t5001000000\ta\t3\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t188: a zero-length A at 100 overlaps the
#  bases on both sides (99 and 100) and nothing further
###########################################################
echo -e "    intersect.opt.t188...\c"
printf 'chr1\t100\t100\tz\n' > a.tmp.bed
printf 'chr1\t98\t99\tfar_left\nchr1\t99\t100\tleft\nchr1\t100\t101\tright\nchr1\t101\t102\tfar_right\n' > b.tmp.bed
printf 'chr1\t100\t100\tz\tchr1\t99\t100\tleft\nchr1\t100\t100\tz\tchr1\t100\t101\tright\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t189: a zero-length A is printed as given in
#  default mode, and -wo gives it 0 bp of overlap
###########################################################
echo -e "    intersect.opt.t189...\c"
printf 'chr1\t150\t150\tz\n' > a.tmp.bed
printf 'chr1\t100\t200\tb\n' > b.tmp.bed
printf 'chr1\t150\t150\tz\n' > exp
printf 'chr1\t150\t150\tz\tchr1\t100\t200\tb\t0\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t190: zero-length A records with and without
#  hits under -v and -c
###########################################################
echo -e "    intersect.opt.t190...\c"
printf 'chr1\t100\t100\thit\nchr1\t500\t500\tmiss\n' > a.tmp.bed
printf 'chr1\t50\t150\tb\n' > b.tmp.bed
printf 'chr1\t500\t500\tmiss\n' > exp
printf 'chr1\t100\t100\thit\t1\nchr1\t500\t500\tmiss\t0\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -v > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -c >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t191: zero-length B records at 100, 150 and
#  200 overlap A 100-200 (each widened by 1 bp); those at 99
#  and 201 only touch it
###########################################################
echo -e "    intersect.opt.t191...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t99\t99\tz99\nchr1\t100\t100\tz100\nchr1\t150\t150\tz150\nchr1\t200\t200\tz200\nchr1\t201\t201\tz201\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t100\t100\tz100\n' > exp
printf 'chr1\t100\t200\ta\tchr1\t150\t150\tz150\n' >> exp
printf 'chr1\t100\t200\ta\tchr1\t200\t200\tz200\n' >> exp
printf 'chr1\t100\t200\ta\t3\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -c >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t192: a zero-length B inside A: -wo gives 0;
#  -u reports A, -v does not
###########################################################
echo -e "    intersect.opt.t192...\c"
printf 'chr1\t100\t200\ta\n' > a.tmp.bed
printf 'chr1\t150\t150\tz\n' > b.tmp.bed
printf 'chr1\t100\t200\ta\tchr1\t150\t150\tz\t0\n' > exp
printf 'chr1\t100\t200\ta\n' >> exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -u >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -v >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t193: zero-length A and B at the same
#  position overlap
###########################################################
echo -e "    intersect.opt.t193...\c"
printf 'chr1\t100\t100\tza\n' > a.tmp.bed
printf 'chr1\t100\t100\tzb\n' > b.tmp.bed
printf 'chr1\t100\t100\tza\tchr1\t100\t100\tzb\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.opt.t194: zero-length records give the same rows
#  with -sorted
###########################################################
echo -e "    intersect.opt.t194...\c"
printf 'chr1\t100\t100\tz1\nchr1\t300\t300\tz2\n' > a.tmp.bed
printf 'chr1\t99\t100\tleft\nchr1\t100\t101\tright\nchr1\t250\t350\tcover\n' > b.tmp.bed
printf 'chr1\t100\t100\tz1\tchr1\t99\t100\tleft\n' > exp
printf 'chr1\t100\t100\tz1\tchr1\t100\t101\tright\n' >> exp
printf 'chr1\t300\t300\tz2\tchr1\t250\t350\tcover\n' >> exp
cat exp exp > exp2
mv exp2 exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb -sorted >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  ADDITIONS
###########################################################
###########################################################

###########################################################
#  intersect.opt.t195: an A record of known strand never
#  pairs with a B record of unknown strand ('.' or '*'):
#  -S and -s report no pair, and -S -v reports both A
#  records (complements t082 and t083, whose unknown-strand
#  B records only overlap unknown-strand A records)
###########################################################
echo -e "    intersect.opt.t195...\c"
printf 'chr1\t100\t200\tplus\t0\t+\nchr1\t300\t400\tminus\t0\t-\n' > a.tmp.bed
printf 'chr1\t150\t160\tdot1\t0\t.\nchr1\t170\t180\tstar1\t0\t*\nchr1\t350\t360\tdot2\t0\t.\nchr1\t370\t380\tstar2\t0\t*\n' > b.tmp.bed
printf 'chr1\t100\t200\tplus\t0\t+\nchr1\t300\t400\tminus\t0\t-\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -S -wa -wb > obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -s -wa -wb >> obs
$BT intersect -a a.tmp.bed -b b.tmp.bed -S -v >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
