set -e;
echo -e \
"\n###########################################################
#
#  INTERSECT FORMATS
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


htsutil=${htsutil-../../htsutil}

###########################################################
###########################################################
#  BED record types: detection from the first data line,
#  and how each type prints in default and -wb output
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t001: BED4 query, default output
#  A chr1:150-250 (4th column not numeric, so BED4); B is
#  chr1:100-200. Default output is A cut to the overlap,
#  150-200, with the name column kept.
###########################################################
echo -e "    intersect.fmt.t001...\c"
printf 'chr1\t150\t250\tn1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t002: bedGraph query keeps its value text
#  The 4th column -1.50e-3 is numeric, so the file is
#  bedGraph. The value is printed exactly as written (not
#  reformatted as a number); coordinates are the overlap.
###########################################################
echo -e "    intersect.fmt.t002...\c"
printf 'chr1\t150\t250\t-1.50e-3\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\t-1.50e-3\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t003: BED5 query, default output
#  Name and score are kept; coordinates are the overlap.
###########################################################
echo -e "    intersect.fmt.t003...\c"
printf 'chr1\t150\t250\tn1\t7\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\t7\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t004: BED6+ query (7 columns, column 6 a
#  strand), default output
#  Columns 4-7 are printed unchanged after the overlap
#  coordinates 150-200.
###########################################################
echo -e "    intersect.fmt.t004...\c"
printf 'chr1\t150\t250\tn1\t7\t-\tx1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\t7\t-\tx1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t005: BED+ query whose column 6 is not a
#  strand (7 columns), default output
#  Column 6 "q" makes this a BED+ record with 3 fixed
#  columns; columns 4-7 are still printed unchanged.
###########################################################
echo -e "    intersect.fmt.t005...\c"
printf 'chr1\t150\t250\tn1\t7\tq\tx1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\t7\tq\tx1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t006: BED9 query, default output
#  Only start and end are cut to the overlap; thickStart,
#  thickEnd and itemRgb (160, 240, 255,0,0) are unchanged.
#  Pins the current output: columns 7-9 are not adjusted.
###########################################################
echo -e "    intersect.fmt.t006...\c"
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\t7\t-\t160\t240\t255,0,0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t007: BED12 query, default output
#  Start and end become the overlap 150-200; the thick
#  columns and the block columns are printed unchanged.
#  Pins the current output: columns 7-12 are not adjusted.
###########################################################
echo -e "    intersect.fmt.t007...\c"
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t200\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t008: 12 columns with a non-numeric block
#  count are not BED12, so -split does not split them
#  The record has blocks 150-160 and 230-250 if read as
#  BED12, but column 10 is "x", so it is a BED6+ record and
#  is one interval 150-250. B 170-180 lies in the "gap"
#  and still hits it with -split -u.
###########################################################
echo -e "    intersect.fmt.t008...\c"
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\tx\t10,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t170\t180\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\tx\t10,20,\t0,80,\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -split -u > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t009: BED12 with a real block count is
#  split by -split (control for t008)
#  Same record with block count 2: blocks 150-160 and
#  230-250. B 170-180 falls between them: no output.
###########################################################
echo -e "    intersect.fmt.t009...\c"
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t170\t180\n' > b.tmp.bed
: > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -split -u > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t010: bedGraph database, -wb
#  A is cut to the overlap 150-200 and the whole B record
#  follows, its value text unchanged.
###########################################################
echo -e "    intersect.fmt.t010...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\t-1.50e-3\n' > b.tmp.bed
printf 'chr1\t150\t200\tchr1\t150\t250\t-1.50e-3\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t011: BED6+ database, -wb
#  The 7-column B record is printed whole after A's
#  overlap.
###########################################################
echo -e "    intersect.fmt.t011...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\tx1\n' > b.tmp.bed
printf 'chr1\t150\t200\tchr1\t150\t250\tn1\t7\t-\tx1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t012: BED12 database, -wb
#  The 12-column B record is printed whole, with its own
#  start and end (150-250), not the overlap.
###########################################################
echo -e "    intersect.fmt.t012...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > b.tmp.bed
printf 'chr1\t150\t200\tchr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t013: BED12 database, -wo without -split
#  counts the overlap of the whole span
#  A 100-200 and B 150-250 overlap 50 bp; without -split
#  the blocks are ignored, so the count is 50, not the
#  10 bp shared with block 150-160.
###########################################################
echo -e "    intersect.fmt.t013...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\t50\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  The null database record printed by -loj, by type
#  A has one record that hits B (chr1:100-200) and one
#  that does not (chr2:1-2). The miss is followed by the
#  null form of B's record type.
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t014: -loj null record for a BED3 database
#  BED3 null: ". -1 -1".
###########################################################
echo -e "    intersect.fmt.t014...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\nchr2\t1\t2\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t015: -loj null record for a BED4 database
#  BED4 null: ". -1 -1 .".
###########################################################
echo -e "    intersect.fmt.t015...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\nchr2\t1\t2\t.\t-1\t-1\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t016: -loj null record for a bedGraph
#  database
#  bedGraph null: ". -1 -1 ." (the value column is ".").
###########################################################
echo -e "    intersect.fmt.t016...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\t0.5\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\t0.5\nchr2\t1\t2\t.\t-1\t-1\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t017: -loj null record for a BED5 database
#  BED5 (numeric score) null: ". -1 -1 . -1".
###########################################################
echo -e "    intersect.fmt.t017...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\nchr2\t1\t2\t.\t-1\t-1\t.\t-1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t018: -loj null record for 6 columns whose
#  strand column is not a strand
#  Column 6 "q" makes B a BED+ record with 3 fixed columns,
#  so its null is ". -1 -1" and one "." per extra column:
#  ". -1 -1 . . ." (score is "." here, not -1 as in BED6).
#  Pins the current output: BED+ null fields after the
#  coordinates are all ".".
###########################################################
echo -e "    intersect.fmt.t018...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\tq\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\tq\nchr2\t1\t2\t.\t-1\t-1\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t019: -loj null record for a BED6+
#  database (7 columns, column 6 a strand)
#  BED6+ null: ". -1 -1" then "." for name, score, strand
#  and the extra column: 7 fields.
#  Pins the current output: BED6+ null fields after the
#  coordinates are all ".", score included.
###########################################################
echo -e "    intersect.fmt.t019...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\tx1\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\t-\tx1\nchr2\t1\t2\t.\t-1\t-1\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t020: -loj null record for a BED+ database
#  (7 columns, column 6 not a strand)
#  BED+ with 3 fixed columns: ". -1 -1" and four ".".
#  Pins the current output: BED+ null fields after the
#  coordinates are all ".".
###########################################################
echo -e "    intersect.fmt.t020...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\tq\tx1\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\tq\tx1\nchr2\t1\t2\t.\t-1\t-1\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t021: -loj null record for a BED12
#  database
#  BED12 null: the BED6 null ". -1 -1 . -1 ." and six ".".
###########################################################
echo -e "    intersect.fmt.t021...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\t2\t10,20,\t0,80,\nchr2\t1\t2\t.\t-1\t-1\t.\t-1\t.\t.\t.\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t022: -loj null record for 12 columns that
#  fail the BED12 check
#  A non-numeric block count makes B BED6+, whose null is
#  ". -1 -1" and nine "." (score "." rather than -1).
#  Pins the current output: BED6+ null fields after the
#  coordinates are all ".", score included.
###########################################################
echo -e "    intersect.fmt.t022...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\tx\t10,20,\t0,80,\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\tx\t10,20,\t0,80,\nchr2\t1\t2\t.\t-1\t-1\t.\t.\t.\t.\t.\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t023: with several databases the -loj null
#  record has the type of the first database
#  The first database is BED4 and the second BED9. The
#  miss gets a "." file column and the BED4 null
#  ". -1 -1 .", not a 9-column null.
#  Pins the current output: the null record takes the
#  type of the first database.
###########################################################
echo -e "    intersect.fmt.t023...\c"
printf 'chr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\n' > c.tmp.bed
printf 'chr2\t1\t2\t.\t.\t-1\t-1\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
#  intersect.fmt.t024: a BED+ null record is as wide as the
#  widest database
#  The first database is BED6+ with 7 columns, the second
#  BED9. The null is BED6+ padded to 9 columns:
#  ". -1 -1" and six ".".
#  Pins the current output: a BED+ null takes its width
#  from the widest database.
###########################################################
echo -e "    intersect.fmt.t024...\c"
printf 'chr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\tx\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t160\t240\t255,0,0\n' > c.tmp.bed
printf 'chr2\t1\t2\t.\t.\t-1\t-1\t.\t.\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed c.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed c.tmp.bed

###########################################################
###########################################################
#  Coordinate parsing and printing
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t025: leading zeros in BED coordinates are
#  parsed and printed without them
#  0150 and 0250 are 150 and 250; BED start and end are
#  printed from the parsed numbers, so -wa prints 150 250.
#  Pins the current output: BED start and end are
#  reprinted as plain integers.
###########################################################
echo -e "    intersect.fmt.t025...\c"
printf 'chr1\t0150\t0250\tn1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t026: leading zeros in a BED6 database
#  record are dropped in -wb output too
#  Pins the current output: BED start and end are
#  reprinted as plain integers.
###########################################################
echo -e "    intersect.fmt.t026...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t0150\t0250\tn1\t7\t+\n' > b.tmp.bed
printf 'chr1\t150\t200\tchr1\t150\t250\tn1\t7\t+\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t027: BED12 thickStart and thickEnd are
#  printed as written
#  Start and end lose their leading zeros (t025), but
#  thickStart 0160 and thickEnd 0240 are printed as text.
#  Pins the current output: only columns 2 and 3 are
#  reformatted.
###########################################################
echo -e "    intersect.fmt.t027...\c"
printf 'chr1\t0150\t0250\tn1\t7\t-\t0160\t0240\t0\t2\t10,20,\t0,80,\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\t7\t-\t0160\t0240\t0\t2\t10,20,\t0,80,\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t028: scientific notation coordinates
#  15e1 and 25e1 are read as 150 and 250 and printed as
#  integers.
#  Pins the current output: exponent notation is accepted
#  in BED coordinates and printed as an integer.
###########################################################
echo -e "    intersect.fmt.t028...\c"
printf 'chr1\t15e1\t25e1\tn1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t029: a "+" sign before a coordinate is an
#  illegal number
#  Pins the current output: the type checker accepts
#  "+150" as an integer, then the parser rejects it.
###########################################################
echo -e "    intersect.fmt.t029...\c"
printf 'chr1\t+150\t250\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '***** ERROR: illegal number "+150". Exiting...\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t030: non-numeric start after the first
#  line
#  Line 1 sets the type to BED3; line 2's start "abc"
#  cannot be parsed: error naming the bad text, exit 1.
#  Line 1 is on chr2 so that it has no hit to report.
###########################################################
echo -e "    intersect.fmt.t030...\c"
printf 'chr2\t150\t250\nchr1\tabc\t250\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '***** ERROR: illegal number "abc". Exiting...\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t031: decimal end coordinate after the
#  first line is an illegal number
#  Line 1 is on chr2 so that it has no hit to report.
###########################################################
echo -e "    intersect.fmt.t031...\c"
printf 'chr2\t150\t250\nchr1\t150\t250.5\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '***** ERROR: illegal number "250.5". Exiting...\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t032: negative start is an invalid record
#  The error names the file and prints the record.
###########################################################
echo -e "    intersect.fmt.t032...\c"
printf 'chr1\t-5\t250\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.bed. Record is \nchr1\t-5\t250\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t033: start greater than end after the
#  first line of the query
#  Line 2 (300 > 250) is an invalid record: the error
#  names the file and prints the record, exit 1. Line 1 is
#  on chr2 so that it has no hit to report.
###########################################################
echo -e "    intersect.fmt.t033...\c"
printf 'chr2\t150\t250\nchr1\t300\t250\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.bed. Record is \nchr1\t300\t250\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t034: start greater than end in the
#  database
#  B is read before any query record, so the error names
#  b.tmp.bed and nothing is printed.
###########################################################
echo -e "    intersect.fmt.t034...\c"
printf 'chr1\t150\t250\n' > a.tmp.bed
printf 'chr1\t100\t200\nchr1\t220\t210\n' > b.tmp.bed
printf 'Error: Invalid record in file b.tmp.bed. Record is \nchr1\t220\t210\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t035: start greater than end on the first
#  line
#  The type checker rejects BED when the first data line
#  has end < start, and no other type fits: the "unable to
#  determine types" error, exit 1.
###########################################################
echo -e "    intersect.fmt.t035...\c"
printf 'chr1\t250\t150\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: unable to open file or unable to determine types for file a.tmp.bed\n\n- Please ensure that your file is TAB delimited (e.g., cat -t FILE).\n- Also ensure that your file has integer chromosome coordinates in the \n  expected columns (e.g., cols 2 and 3 for BED).\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t036: non-numeric start on the first line
#  Same type-detection error as t035, here for the
#  database file.
###########################################################
echo -e "    intersect.fmt.t036...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tabc\t200\n' > b.tmp.bed
printf 'Error: unable to open file or unable to determine types for file b.tmp.bed\n\n- Please ensure that your file is TAB delimited (e.g., cat -t FILE).\n- Also ensure that your file has integer chromosome coordinates in the \n  expected columns (e.g., cols 2 and 3 for BED).\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t037: two columns are too few for any type
###########################################################
echo -e "    intersect.fmt.t037...\c"
printf 'chr1\t150\nchr1\t160\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: unable to open file or unable to determine types for file a.tmp.bed\n\n- Please ensure that your file is TAB delimited (e.g., cat -t FILE).\n- Also ensure that your file has integer chromosome coordinates in the \n  expected columns (e.g., cols 2 and 3 for BED).\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t038: space-delimited columns are not
#  recognised
#  With no tab the line is one field, so no type fits.
###########################################################
echo -e "    intersect.fmt.t038...\c"
printf 'chr1 150 250\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: unable to open file or unable to determine types for file a.tmp.bed\n\n- Please ensure that your file is TAB delimited (e.g., cat -t FILE).\n- Also ensure that your file has integer chromosome coordinates in the \n  expected columns (e.g., cols 2 and 3 for BED).\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t039: a binary file that is not BAM
#  More than 10% of the bytes are not printable, so the
#  file is binary; it is not BAM, so its type is unknown.
###########################################################
echo -e "    intersect.fmt.t039...\c"
printf '\001\002\003\004\005\006\007\010\016\017\020\021\022\023\024\025\026\027\030\031\n' > a.tmp.bin
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Unable to determine type for file a.tmp.bin\nexit 1\n' > exp
$BT intersect -a a.tmp.bin -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bin b.tmp.bed

###########################################################
#  intersect.fmt.t040: different column counts within the
#  first data lines
#  Lines 1 and 2 have 3 and 4 columns, so the type checker
#  cannot settle on a column count.
###########################################################
echo -e "    intersect.fmt.t040...\c"
printf 'chr1\t150\t250\nchr1\t150\t250\tn2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Type checker found wrong number of fields while tokenizing data line.\nPerhaps you have extra TAB at the end of your line? Check with "cat -t"\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t042: a missing column after the lines the
#  type checker reads
#  Line 5 has only 2 columns.
###########################################################
echo -e "    intersect.fmt.t042...\c"
printf 'chr1\t1\t2\nchr1\t3\t4\nchr1\t5\t6\nchr1\t7\t8\nchr1\t9\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: line number 5 of file a.tmp.bed has 2 fields, but 3 were expected.\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t044: a database file that does not exist
###########################################################
echo -e "    intersect.fmt.t044...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'Error: Unable to open file nosuch.tmp.bed. Exiting.\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b nosuch.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t045: BED12 block lists that disagree with
#  the block count stop -split
#  Block count 2 but one block size: with -split the
#  record cannot be split; the tool exits with status 255
#  (exit(-1)) and an error.
###########################################################
echo -e "    intersect.fmt.t045...\c"
printf 'chr1\t150\t250\tn1\t0\t+\t150\t250\t0\t2\t10,\t0,80,\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: found wrong block counts while splitting entry.\nexit 255\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -split > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t046: without -split the same BED12 block
#  lists are not checked
#  The record is printed as an ordinary interval.
###########################################################
echo -e "    intersect.fmt.t046...\c"
printf 'chr1\t150\t250\tn1\t0\t+\t150\t250\t0\t2\t10,\t0,80,\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t250\tn1\t0\t+\t150\t250\t0\t2\t10,\t0,80,\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  GFF and GTF: 1-based, end-inclusive coordinates
#  GFF start 151, end 250 is the 0-based interval 150-250.
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t047: GFF query, default output prints the
#  overlap in 1-based coordinates
#  GFF 151-250 is 150-250; B 160-170 overlaps 160-170,
#  which GFF writes as 161-170. Every other column
#  (score 12.5, frame 0, attributes with a space) is
#  printed unchanged.
###########################################################
echo -e "    intersect.fmt.t047...\c"
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > a.tmp.gff
printf 'chr1\t160\t170\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t161\t170\t12.5\t+\t0\tID=g1;Name=foo bar\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t048: GFF query, -wa prints the record
#  exactly as read
###########################################################
echo -e "    intersect.fmt.t048...\c"
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > a.tmp.gff
printf 'chr1\t160\t170\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t049: GFF query, -wo counts overlap bases
#  after the 1-based conversion
#  150-250 and 100-200 share 50 bases.
###########################################################
echo -e "    intersect.fmt.t049...\c"
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\tchr1\t100\t200\t50\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t050: GFF database, -wb
#  A 100-200 is cut to the overlap 150-200 (0-based, it is
#  a BED record) and the GFF record follows unchanged.
###########################################################
echo -e "    intersect.fmt.t050...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > b.tmp.gff
printf 'chr1\t150\t200\tchr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t051: GFF database, -wo
#  50 overlapping bases, as in t049 from the other side.
###########################################################
echo -e "    intersect.fmt.t051...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\n' > b.tmp.gff
printf 'chr1\t100\t200\tchr1\tsrc\tgene\t151\t250\t12.5\t+\t0\tID=g1;Name=foo bar\t50\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -wo > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t052: the GFF start base is included
#  GFF start 151 is base 150-151. BED 149-150 ends just
#  before it (no hit); BED 150-151 is that base (hit).
###########################################################
echo -e "    intersect.fmt.t052...\c"
printf 'chr1\t149\t150\nchr1\t150\t151\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t149\t150\t0\nchr1\t150\t151\t1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t053: the GFF end base is included
#  GFF end 250 is base 249-250: BED 249-250 hits, BED
#  250-251 does not.
###########################################################
echo -e "    intersect.fmt.t053...\c"
printf 'chr1\t249\t250\nchr1\t250\t251\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t249\t250\t1\nchr1\t250\t251\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t054: a GFF feature with start == end is
#  one base long
#  GFF 151-151 is base 150-151: only BED 150-151 hits it.
###########################################################
echo -e "    intersect.fmt.t054...\c"
printf 'chr1\t149\t150\nchr1\t150\t151\nchr1\t151\t152\n' > a.tmp.bed
printf 'chr1\tsrc\tSNP\t151\t151\t.\t+\t.\tID=s1\n' > b.tmp.gff
printf 'chr1\t149\t150\t0\nchr1\t150\t151\t1\nchr1\t151\t152\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t055: 8-column GFF (no attributes) query
#  The output has 8 columns and no trailing tab.
###########################################################
echo -e "    intersect.fmt.t055...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t-\t.\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t200\t.\t-\t.\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t056: -loj null record for an 8-column GFF
#  database
#  GFF null: ". . . -1 -1 . . ." (8 columns).
###########################################################
echo -e "    intersect.fmt.t056...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t-\t.\n' > b.tmp.gff
printf 'chr1\t100\t200\tchr1\tsrc\tgene\t151\t250\t.\t-\t.\nchr2\t1\t2\t.\t.\t.\t-1\t-1\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t057: -loj null record for a 9-column GFF
#  database
#  The null gains a "." for the attribute column.
###########################################################
echo -e "    intersect.fmt.t057...\c"
printf 'chr1\t100\t200\nchr2\t1\t2\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t-\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t100\t200\tchr1\tsrc\tgene\t151\t250\t.\t-\t.\tID=g1\nchr2\t1\t2\t.\t.\t.\t-1\t-1\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  intersect.fmt.t058: GFF query with -c
#  The count follows the unchanged GFF record.
###########################################################
echo -e "    intersect.fmt.t058...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > a.tmp.gff
printf 'chr1\t100\t160\nchr1\t240\t300\nchr1\t250\t260\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\t2\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t059: GFF query with -u
#  Only the feature overlapping B (101-200, i.e. 100-200)
#  is printed, once, unchanged.
###########################################################
echo -e "    intersect.fmt.t059...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t301\t400\t.\t+\t.\tID=g2\n' > a.tmp.gff
printf 'chr1\t150\t160\nchr1\t170\t180\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -u > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t060: GFF query with -v
#  Only the feature without a hit (301-400) is printed.
###########################################################
echo -e "    intersect.fmt.t060...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t301\t400\t.\t+\t.\tID=g2\n' > a.tmp.gff
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t301\t400\t.\t+\t.\tID=g2\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -v > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t061: GFF against GFF, -wo
#  151-250 (150-250) and 201-300 (200-300) share 200-250:
#  50 bases.
###########################################################
echo -e "    intersect.fmt.t061...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=a\n' > a.tmp.gff
printf 'chr1\tsrc\texon\t201\t300\t.\t+\t.\tID=b\n' > b.tmp.gff
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=a\tchr1\tsrc\texon\t201\t300\t.\t+\t.\tID=b\t50\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.gff -wo > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.gff

###########################################################
#  intersect.fmt.t062: GFF against GFF, default output
#  The overlap 200-250 is written 1-based as 201-250.
###########################################################
echo -e "    intersect.fmt.t062...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=a\n' > a.tmp.gff
printf 'chr1\tsrc\texon\t201\t300\t.\t+\t.\tID=b\n' > b.tmp.gff
printf 'chr1\tsrc\tgene\t201\t250\t.\t+\t.\tID=a\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.gff > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.gff

###########################################################
#  intersect.fmt.t063: GTF query, default output
#  GTF is read as GFF: the overlap with 100-200 is written
#  151-200 and the attribute column (quotes, "; "
#  separators, trailing ";") is printed unchanged.
###########################################################
echo -e "    intersect.fmt.t063...\c"
printf 'chr1\tHAVANA\texon\t151\t250\t.\t-\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";\n' > a.tmp.gtf
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tHAVANA\texon\t151\t200\t.\t-\t.\tgene_id "G1"; transcript_id "T1"; exon_number "1";\n' > exp
$BT intersect -a a.tmp.gtf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gtf b.tmp.bed

###########################################################
#  intersect.fmt.t064: GFF directives and -header
#  The leading "##gff-version 3" and "#!genome-build"
#  lines are the header and are printed by -header; the
#  "###" line between features is skipped and not printed.
###########################################################
echo -e "    intersect.fmt.t064...\c"
printf '##gff-version 3\n#!genome-build GRCh38\nchr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n###\nchr1\tsrc\tgene\t171\t180\t.\t+\t.\tID=g2\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '##gff-version 3\n#!genome-build GRCh38\nchr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\nchr1\tsrc\tgene\t171\t180\t.\t+\t.\tID=g2\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -u -header > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t065: GFF start greater than end after the
#  first line
#  300 > 200 on line 2: invalid record, printed as read
#  (1-based), exit 1.
###########################################################
echo -e "    intersect.fmt.t065...\c"
printf 'chr2\tsrc\tgene\t151\t250\t.\t+\t.\tID=a\nchr1\tsrc\tgene\t300\t200\t.\t+\t.\tID=b\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.gff. Record is \nchr1\tsrc\tgene\t300\t200\t.\t+\t.\tID=b\nexit 1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t066: GFF start 0 is invalid
#  1-based start 0 becomes -1 after conversion: invalid
#  record, exit 1.
###########################################################
echo -e "    intersect.fmt.t066...\c"
printf 'chr1\tsrc\tgene\t0\t250\t.\t+\t.\tID=a\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.gff. Record is \nchr1\tsrc\tgene\t0\t250\t.\t+\t.\tID=a\nexit 1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t067: non-numeric GFF start after the
#  first line
###########################################################
echo -e "    intersect.fmt.t067...\c"
printf 'chr2\tsrc\tgene\t151\t250\t.\t+\t.\tID=a\nchr1\tsrc\tgene\tx\t200\t.\t+\t.\tID=b\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '***** ERROR: illegal number "x". Exiting...\nexit 1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
###########################################################
#  VCF: POS is 1-based; a record covers
#  [POS-1, POS-1+len(REF)), or for symbolic ALTs the
#  length given by INFO END or SVLEN
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t068: VCF record length comes from REF
#  Against B 0-1000, -wo reports each record's length:
#  SNP A>G 1, deletion ACGT>A 4, insertion A>ACGTT 1 (REF
#  is one base), multi-allelic AC>A,ACC 2 (REF length,
#  not either ALT).
###########################################################
echo -e "    intersect.fmt.t068...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t201\tdel\tACGT\tA\t.\t.\t.\nchr1\t301\tins\tA\tACGTT\t.\t.\t.\nchr1\t401\tma\tAC\tA,ACC\t.\t.\t.\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tchr1\t0\t1000\t1\nchr1\t201\tdel\tACGT\tA\t.\t.\t.\tchr1\t0\t1000\t4\nchr1\t301\tins\tA\tACGTT\t.\t.\t.\tchr1\t0\t1000\t1\nchr1\t401\tma\tAC\tA,ACC\t.\t.\t.\tchr1\t0\t1000\t2\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t069: a SNP at POS 101 is base 100-101
#  BED 99-100 and 101-102 are the neighbouring bases: only
#  100-101 hits.
###########################################################
echo -e "    intersect.fmt.t069...\c"
printf 'chr1\t99\t100\nchr1\t100\t101\nchr1\t101\t102\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\t.\n' >> b.tmp.vcf
printf 'chr1\t99\t100\t0\nchr1\t100\t101\t1\nchr1\t101\t102\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t070: a deletion covers every REF base
#  REF ACGT at POS 201 is 200-204: BED 199-200 and 204-205
#  miss, 200-201 and 203-204 hit.
###########################################################
echo -e "    intersect.fmt.t070...\c"
printf 'chr1\t199\t200\nchr1\t200\t201\nchr1\t203\t204\nchr1\t204\t205\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t201\tdel\tACGT\tA\t.\t.\t.\n' >> b.tmp.vcf
printf 'chr1\t199\t200\t0\nchr1\t200\t201\t1\nchr1\t203\t204\t1\nchr1\t204\t205\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t071: a symbolic ALT ends at INFO END
#  <DEL> at 301 with END=400 is 300-400: the REF base at
#  POS through END, 100 bases.
###########################################################
echo -e "    intersect.fmt.t071...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=400\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;END=400\tchr1\t0\t1000\t100\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t072: a symbolic ALT without END takes
#  its length from SVLEN, ignoring the sign
#  <DEL> at 501 with SVLEN=-50 starts at 500 and reaches
#  about 50 bases; <DUP> at 601 with SVLEN=10 about 10.
#  Probes well inside (520, 608) hit and probes well past
#  the ends (560, 620) do not; the probes stay clear of
#  the last base, whose inclusion is not tested here.
###########################################################
echo -e "    intersect.fmt.t072...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t501\tsv2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50\nchr1\t601\tsv3\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\n' >> a.tmp.vcf
printf 'chr1\t520\t521\tp520\nchr1\t560\t561\tp560\nchr1\t608\t609\tp608\nchr1\t620\t621\tp620\n' > b.tmp.bed
printf 'chr1\t501\tsv2\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=-50\tchr1\t520\t521\tp520\nchr1\t601\tsv3\tN\t<DUP>\t.\tPASS\tSVTYPE=DUP;SVLEN=10\tchr1\t608\t609\tp608\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t074: several SVLEN values: the one with
#  the largest absolute value sets the length
#  SVLEN=-5,20,-30 at 301: with 30 the record reaches
#  past 325, so probe 325-326 hits; -5 or 20 would end
#  before it.
###########################################################
echo -e "    intersect.fmt.t074...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t301\tcnv\tN\t<CNV>,<DEL>,<DUP>\t.\tPASS\tSVLEN=-5,20,-30\n' >> a.tmp.vcf
printf 'chr1\t325\t326\n' > b.tmp.bed
printf 'chr1\t301\tcnv\tN\t<CNV>,<DEL>,<DUP>\t.\tPASS\tSVLEN=-5,20,-30\t1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t075: an <INS> is a zero-length record
#  that hits the base it is attached to
#  <INS> at POS 601 hits BED 600-601 (base 601), and not
#  598-599 or 602-603, two bases away on either side.
###########################################################
echo -e "    intersect.fmt.t075...\c"
printf 'chr1\t598\t599\nchr1\t600\t601\nchr1\t602\t603\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t601\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS\n' >> b.tmp.vcf
printf 'chr1\t598\t599\t0\nchr1\t600\t601\t1\nchr1\t602\t603\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t076: -wo reports 0 bases for an <INS>
#  query
###########################################################
echo -e "    intersect.fmt.t076...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t601\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'chr1\t601\tins\tN\t<INS>\t.\tPASS\tSVTYPE=INS\tchr1\t0\t1000\t0\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t077: any ALT starting "<INS" is an
#  insertion, and its SVLEN is not used
#  <INS:ME:ALU> at 601 with SVLEN=300 stays a point; B
#  610-620 would hit a 300-base record but has no hit.
###########################################################
echo -e "    intersect.fmt.t077...\c"
printf 'chr1\t610\t620\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t601\talu\tN\t<INS:ME:ALU>\t.\tPASS\tSVTYPE=INS;SVLEN=300\n' >> b.tmp.vcf
printf 'chr1\t610\t620\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t078: a symbolic ALT with neither END nor
#  SVLEN is an invalid record
#  The length cannot be determined; the error prints the
#  record as read, exit 1.
###########################################################
echo -e "    intersect.fmt.t078...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.vcf. Record is \nchr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL\nexit 1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t079: a breakend ALT is not symbolic
#  N]chr2:100] does not start with "<", so the length is
#  that of REF: 1 base.
###########################################################
echo -e "    intersect.fmt.t079...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tbnd\tN\tN]chr2:100]\t.\tPASS\tSVTYPE=BND\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'chr1\t101\tbnd\tN\tN]chr2:100]\t.\tPASS\tSVTYPE=BND\tchr1\t0\t1000\t1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wo > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t080: VCF query, default output prints the
#  whole record
#  The deletion covers 300-400 and B 350-360 only part of
#  it, but a VCF record is never cut to the overlap.
###########################################################
echo -e "    intersect.fmt.t080...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\n' >> a.tmp.vcf
printf 'chr1\t350\t360\n' > b.tmp.bed
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t081: VCF query, -wb
#  The whole VCF record, then the whole B record.
###########################################################
echo -e "    intersect.fmt.t081...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\n' >> a.tmp.vcf
printf 'chr1\t350\t360\tb1\n' > b.tmp.bed
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\tchr1\t350\t360\tb1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t082: VCF database, -wb
#  A (BED) is cut to the overlap with 300-400, 350-400,
#  and the VCF record follows unchanged.
###########################################################
echo -e "    intersect.fmt.t082...\c"
printf 'chr1\t350\t500\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\n' >> b.tmp.vcf
printf 'chr1\t350\t400\tchr1\t301\tsv1\tN\t<DEL>\t.\tPASS\tEND=400\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t083: -loj null record for an 8-column VCF
#  database
#  VCF null: ". -1 ." and one "." for each of the other 5
#  columns.
###########################################################
echo -e "    intersect.fmt.t083...\c"
printf 'chr1\t100\t101\nchr2\t1\t2\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' >> b.tmp.vcf
printf 'chr1\t100\t101\tchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr2\t1\t2\t.\t-1\t.\t.\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t084: -loj null record for a VCF database
#  with FORMAT and two sample columns
#  11 columns, so the null is ". -1 ." and eight ".".
###########################################################
echo -e "    intersect.fmt.t084...\c"
printf 'chr1\t100\t101\nchr2\t1\t2\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n' > b.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:7\n' >> b.tmp.vcf
printf 'chr1\t100\t101\tchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:7\nchr2\t1\t2\t.\t-1\t.\t.\t.\t.\t.\t.\t.\t.\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t085: FORMAT and sample columns of a VCF
#  query are printed unchanged
###########################################################
echo -e "    intersect.fmt.t085...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:7\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tGT:DP\t0/1:5\t1/1:7\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t086: -header prints the VCF meta lines
#  and the #CHROM line once, before the records
###########################################################
echo -e "    intersect.fmt.t086...\c"
printf '##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t151\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\tb1\n' > b.tmp.bed
printf '##fileformat=VCFv4.2\n##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tchr1\t0\t1000\tb1\nchr1\t151\tsnp2\tC\tT\t50\tPASS\tDP=3\tchr1\t0\t1000\tb1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wa -wb -header > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t087: -header with -v on a VCF query
#  The header, then only the variant with no hit (151).
###########################################################
echo -e "    intersect.fmt.t087...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t151\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
printf 'chr1\t90\t110\n' > b.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t151\tsnp2\tC\tT\t50\tPASS\tDP=3\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -v -header > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t088: VCF query with -c
#  The deletion 200-204 is hit by 199-201 and 203-210 but
#  not by 204-210: count 2 after the unchanged record.
###########################################################
echo -e "    intersect.fmt.t088...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t201\tdel\tACGT\tA\t.\t.\t.\n' >> a.tmp.vcf
printf 'chr1\t199\t201\nchr1\t203\t210\nchr1\t204\t210\n' > b.tmp.bed
printf 'chr1\t201\tdel\tACGT\tA\t.\t.\t.\t2\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t089: VCF query with -u
#  The variant with two hits is printed once; the one with
#  none is not printed.
###########################################################
echo -e "    intersect.fmt.t089...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t151\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
printf 'chr1\t90\t110\nchr1\t100\t101\n' > b.tmp.bed
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -u > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t090: strand-aware intersect refuses VCF
#  VCF has no strand column; with -s the error names the
#  VCF file (here the database), exit 1.
###########################################################
echo -e "    intersect.fmt.t090...\c"
printf 'chr1\t0\t1000\tn\t0\t+\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' >> b.tmp.vcf
printf '***** ERROR: stranded intersect not supported for VCF file b.tmp.vcf. *****\nexit 1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.vcf -s > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
tail -n 2 obs > obs2
check obs2 exp
rm obs obs2 exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t091: a VCF without the ##fileformat line
#  is not recognised
#  Only "##fileformat=VCF" at the start of the file marks
#  VCF; read as text, column 3 "." is not a coordinate,
#  so no type fits.
###########################################################
echo -e "    intersect.fmt.t091...\c"
printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\t.\tA\tG\t50\tPASS\tDP=10\n' > a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'Error: unable to open file or unable to determine types for file a.tmp.vcf\n\n- Please ensure that your file is TAB delimited (e.g., cat -t FILE).\n- Also ensure that your file has integer chromosome coordinates in the \n  expected columns (e.g., cols 2 and 3 for BED).\nexit 1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t092: VCF POS 0 is an invalid record
###########################################################
echo -e "    intersect.fmt.t092...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t0\t.\tA\tG\t.\t.\t.\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf 'Error: Invalid record in file a.tmp.vcf. Record is \nchr1\t0\t.\tA\tG\t.\t.\t.\nexit 1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t093: a non-numeric POS in the first VCF
#  record is an illegal number
#  The ##fileformat line already fixed the type, so the
#  error comes from parsing POS, not from type detection.
###########################################################
echo -e "    intersect.fmt.t093...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\tx\t.\tA\tG\t.\t.\t.\n' >> a.tmp.vcf
printf 'chr1\t0\t1000\n' > b.tmp.bed
printf '***** ERROR: illegal number "x". Exiting...\nexit 1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t095: VCF query read from stdin as "-",
#  with -header
#  Type detection works on the stream: the header and the
#  hit are printed as from a file.
###########################################################
echo -e "    intersect.fmt.t095...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t501\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
printf 'chr1\t0\t200\n' > b.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' > exp
cat a.tmp.vcf | $BT intersect -a - -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t096: VCF database read from stdin as
#  "stdin"
###########################################################
echo -e "    intersect.fmt.t096...\c"
printf 'chr1\t0\t200\n' > a.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > b.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t501\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> b.tmp.vcf
printf 'chr1\t0\t200\tchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' > exp
cat b.tmp.vcf | $BT intersect -a a.tmp.bed -b stdin -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  intersect.fmt.t097: gzip-compressed VCF query with
#  -header
###########################################################
echo -e "    intersect.fmt.t097...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t501\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
gzip -c a.tmp.vcf > a.tmp.vcf.gz
printf 'chr1\t0\t200\n' > b.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\n' > exp
$BT intersect -a a.tmp.vcf.gz -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.vcf a.tmp.vcf.gz b.tmp.bed

###########################################################
#  intersect.fmt.t098: VCF with CRLF line endings
#  The carriage returns are dropped from the header and
#  the record, and END=170 is still parsed: 150-170
#  overlaps B 100-200 by 20 bases.
###########################################################
echo -e "    intersect.fmt.t098...\c"
printf '##fileformat=VCFv4.2\r\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\r\nchr1\t151\tsv\tN\t<DEL>\t.\t.\tEND=170\r\n' > a.tmp.vcf
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t151\tsv\tN\t<DEL>\t.\t.\tEND=170\tchr1\t100\t200\t20\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -wo -header > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
###########################################################
#  BAM input
#  reads.tmp.bam (from the SAM below) is used by every
#  test in this section and removed at its end.
#    r1   chr1 101 (100-150)  50M          forward, MAPQ 30
#    r2   chr1 201 (200-350)  20M100N30M   reverse, MAPQ 40
#    r3   chr1 501 (500-560)  20M10D30M
#    r4   chr1 701/751        10M          a pair, mates /1 /2
#    r5   unmapped, no position
#    r6   chr1 801 (800-840)  5S20M5I20M5H
#    r7   chr2 101 (100-110)  10M
#    r8   chr1 901 (900-920)  5=1X4=2P10M
#    r9   chr1 1001 (1000-1135) 10M5D10M100N10M
#    r10  unmapped, placed at chr1 1501
#  mapped.tmp.bam holds the same reads without r5 and r10.
###########################################################
###########################################################
printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:10000\n@CO\tformats test reads\n' > reads.tmp.sam
printf 'r1\t0\tchr1\t101\t30\t50M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r2\t16\tchr1\t201\t40\t20M100N30M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r3\t0\tchr1\t501\t50\t20M10D30M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r4\t99\tchr1\t701\t60\t10M\t=\t751\t60\t*\t*\n' >> reads.tmp.sam
printf 'r4\t147\tchr1\t751\t60\t10M\t=\t701\t-60\t*\t*\n' >> reads.tmp.sam
printf 'r5\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r6\t0\tchr1\t801\t20\t5S20M5I20M5H\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r7\t0\tchr2\t101\t0\t10M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r8\t0\tchr1\t901\t10\t5=1X4=2P10M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r9\t0\tchr1\t1001\t10\t10M5D10M100N10M\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
printf 'r10\t4\tchr1\t1501\t0\t*\t*\t0\t0\t*\t*\n' >> reads.tmp.sam
$htsutil samtobam reads.tmp.sam reads.tmp.bam
awk '$1 != "r5" && $1 != "r10"' reads.tmp.sam > mapped.tmp.sam
$htsutil samtobam mapped.tmp.sam mapped.tmp.bam

###########################################################
#  intersect.fmt.t099: -bed prints a BAM alignment as BED12
#  r1 at 1-based 101 with 50M is 100-150. Score is MAPQ
#  30, strand +, thickStart/thickEnd are the alignment
#  ends, itemRgb 0,0,0 and one block of 50.
###########################################################
echo -e "    intersect.fmt.t099...\c"
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t30\t+\t100\t150\t0,0,0\t1\t50,\t0,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t100: -bed, a reverse-strand spliced read
#  Flag 16 gives strand -. 20M100N30M is split at the N
#  into blocks of 20 and 30 bases starting 0 and 120 bases
#  after 200.
###########################################################
echo -e "    intersect.fmt.t100...\c"
printf 'chr1\t205\t210\n' > b.tmp.bed
printf 'chr1\t200\t350\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t101: -bed, a deletion does not split a
#  block
#  20M10D30M covers 500-560 as one block of 60.
###########################################################
echo -e "    intersect.fmt.t101...\c"
printf 'chr1\t505\t510\n' > b.tmp.bed
printf 'chr1\t500\t560\tr3\t50\t+\t500\t560\t0,0,0\t1\t60,\t0,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t102: -bed, mates of a pair are named
#  with /1 and /2
#  Flag 99 is the first mate (+), flag 147 the second
#  mate (-).
###########################################################
echo -e "    intersect.fmt.t102...\c"
printf 'chr1\t700\t760\n' > b.tmp.bed
printf 'chr1\t700\t710\tr4/1\t60\t+\t700\t710\t0,0,0\t1\t10,\t0,\nchr1\t750\t760\tr4/2\t60\t-\t750\t760\t0,0,0\t1\t10,\t0,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t103: -bed, clips and insertions take no
#  reference bases
#  5S20M5I20M5H spans 40 reference bases, 800-840, in one
#  block.
###########################################################
echo -e "    intersect.fmt.t103...\c"
printf 'chr1\t805\t810\n' > b.tmp.bed
printf 'chr1\t800\t840\tr6\t20\t+\t800\t840\t0,0,0\t1\t40,\t0,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t104: -bed, = and X count as aligned
#  bases and P does not
#  5=1X4=2P10M spans 5+1+4+10 = 20 bases, 900-920.
###########################################################
echo -e "    intersect.fmt.t104...\c"
printf 'chr1\t905\t910\n' > b.tmp.bed
printf 'chr1\t900\t920\tr8\t10\t+\t900\t920\t0,0,0\t1\t20,\t0,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t105: -bed, D joins and N splits blocks
#  10M5D10M100N10M: the first block is 10+5+10 = 25 bases
#  (1000-1025), then 100 skipped, then 10 (1125-1135).
###########################################################
echo -e "    intersect.fmt.t105...\c"
printf 'chr1\t1010\t1011\n' > b.tmp.bed
printf 'chr1\t1000\t1135\tr9\t10\t+\t1000\t1135\t0,0,0\t2\t25,10,\t0,125,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wa > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t106: -bed default output cuts only start
#  and end to the overlap
#  B 210-330 overlaps r2 at 210-330. thickStart and
#  thickEnd stay 200 and 350 and the blocks still describe
#  the whole alignment.
#  Pins the current output: columns 7-12 are not adjusted.
###########################################################
echo -e "    intersect.fmt.t106...\c"
printf 'chr1\t210\t330\n' > b.tmp.bed
printf 'chr1\t210\t330\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t107: -bed -wb
#  The cut BED12 form of r2 followed by the B record.
###########################################################
echo -e "    intersect.fmt.t107...\c"
printf 'chr1\t210\t330\tb1\n' > b.tmp.bed
printf 'chr1\t210\t330\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\tchr1\t210\t330\tb1\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wb > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t108: -bed -wo without -split counts the
#  whole alignment span
#  B 210-330 lies within r2's span 200-350: 120 bases,
#  although only 10+10 of them are aligned.
###########################################################
echo -e "    intersect.fmt.t108...\c"
printf 'chr1\t210\t330\tb1\n' > b.tmp.bed
printf 'chr1\t200\t350\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\tchr1\t210\t330\tb1\t120\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -wo > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t109: -split -wo counts bases shared with
#  the blocks, with D inside a block
#  r9 blocks 1000-1025 and 1125-1135; B 1010-1130 shares
#  15 + 5 = 20 bases (the 5D at 1010-1015 counts, the N
#  does not).
###########################################################
echo -e "    intersect.fmt.t109...\c"
printf 'chr1\t1010\t1130\n' > b.tmp.bed
printf 'chr1\t1000\t1135\tr9\t10\t+\t1000\t1135\t0,0,0\t2\t25,10,\t0,125,\tchr1\t1010\t1130\t20\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -split -wo > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t110: -split does not split at a deletion
#  B 522-528 lies inside r3's 10D (520-530). With -split
#  it is still a hit of 6 bases, because D stays inside
#  the block.
###########################################################
echo -e "    intersect.fmt.t110...\c"
printf 'chr1\t522\t528\n' > b.tmp.bed
printf 'chr1\t500\t560\tr3\t50\t+\t500\t560\t0,0,0\t1\t60,\t0,\tchr1\t522\t528\t6\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -split -wo > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t111: -bed -c prints every read with its
#  count
#  Only r2 overlaps B 205-210; every other mapped read,
#  the chr2 read included, gets 0.
###########################################################
echo -e "    intersect.fmt.t111...\c"
printf 'chr1\t205\t210\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t30\t+\t100\t150\t0,0,0\t1\t50,\t0,\t0\n' > exp
printf 'chr1\t200\t350\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\t1\n' >> exp
printf 'chr1\t500\t560\tr3\t50\t+\t500\t560\t0,0,0\t1\t60,\t0,\t0\n' >> exp
printf 'chr1\t700\t710\tr4/1\t60\t+\t700\t710\t0,0,0\t1\t10,\t0,\t0\n' >> exp
printf 'chr1\t750\t760\tr4/2\t60\t-\t750\t760\t0,0,0\t1\t10,\t0,\t0\n' >> exp
printf 'chr1\t800\t840\tr6\t20\t+\t800\t840\t0,0,0\t1\t40,\t0,\t0\n' >> exp
printf 'chr2\t100\t110\tr7\t0\t+\t100\t110\t0,0,0\t1\t10,\t0,\t0\n' >> exp
printf 'chr1\t900\t920\tr8\t10\t+\t900\t920\t0,0,0\t1\t20,\t0,\t0\n' >> exp
printf 'chr1\t1000\t1135\tr9\t10\t+\t1000\t1135\t0,0,0\t2\t25,10,\t0,125,\t0\n' >> exp
$BT intersect -a mapped.tmp.bam -b b.tmp.bed -bed -c > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t112: -bed -u prints the overlapping read
#  once as BED12
#  B has two records over r2; r2 is printed once.
###########################################################
echo -e "    intersect.fmt.t112...\c"
printf 'chr1\t205\t210\nchr1\t330\t340\n' > b.tmp.bed
printf 'chr1\t200\t350\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -bed -u > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t113: -bed -loj
#  r2 is printed with its B hit; every other read is
#  followed by the BED4 null ". -1 -1 .".
###########################################################
echo -e "    intersect.fmt.t113...\c"
printf 'chr1\t205\t210\tb1\n' > b.tmp.bed
printf 'chr1\t100\t150\tr1\t30\t+\t100\t150\t0,0,0\t1\t50,\t0,\t.\t-1\t-1\t.\n' > exp
printf 'chr1\t200\t350\tr2\t40\t-\t200\t350\t0,0,0\t2\t20,30,\t0,120,\tchr1\t205\t210\tb1\n' >> exp
printf 'chr1\t500\t560\tr3\t50\t+\t500\t560\t0,0,0\t1\t60,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr1\t700\t710\tr4/1\t60\t+\t700\t710\t0,0,0\t1\t10,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr1\t750\t760\tr4/2\t60\t-\t750\t760\t0,0,0\t1\t10,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr1\t800\t840\tr6\t20\t+\t800\t840\t0,0,0\t1\t40,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr2\t100\t110\tr7\t0\t+\t100\t110\t0,0,0\t1\t10,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr1\t900\t920\tr8\t10\t+\t900\t920\t0,0,0\t1\t20,\t0,\t.\t-1\t-1\t.\n' >> exp
printf 'chr1\t1000\t1135\tr9\t10\t+\t1000\t1135\t0,0,0\t2\t25,10,\t0,125,\t.\t-1\t-1\t.\n' >> exp
$BT intersect -a mapped.tmp.bam -b b.tmp.bed -bed -loj > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t114: -wb without -bed is ignored for BAM
#  A warning goes to stderr and the output is BAM holding
#  the overlapping read only. obs is the decoded records
#  followed by stderr (also in t116).
###########################################################
echo -e "    intersect.fmt.t114...\c"
printf 'chr1\t205\t210\tb1\n' > b.tmp.bed
printf 'r2\t16\tchr1\t201\t40\t20M100N30M\t*\t0\t0\t*\t*\n' > exp
printf '\n*****\n*****WARNING: -wb and -loj are ignored with bam input, unless bed output is specified with -bed option.\n*****\n' > experr
$BT intersect -a reads.tmp.bam -b b.tmp.bed -wb 2> obserr | $htsutil viewbamrecords > obs
cat obserr >> obs
cat experr >> exp
check obs exp
rm obs exp obserr experr b.tmp.bed

###########################################################
#  intersect.fmt.t116: -header is ignored for BAM input
#  A warning, and the same BAM records as without -header.
###########################################################
echo -e "    intersect.fmt.t116...\c"
printf 'chr1\t205\t210\n' > b.tmp.bed
printf 'r2\t16\tchr1\t201\t40\t20M100N30M\t*\t0\t0\t*\t*\n' > exp
printf '\n*****\n*****WARNING: -header option is not valid for BAM input.\n*****\n' > experr
$BT intersect -a reads.tmp.bam -b b.tmp.bed -header 2> obserr | $htsutil viewbamrecords > obs
cat obserr >> obs
cat experr >> exp
check obs exp
rm obs exp obserr experr b.tmp.bed

###########################################################
#  intersect.fmt.t117: -c with a BAM query needs -bed
#  A count cannot be added to a BAM record: error, exit 1.
###########################################################
echo -e "    intersect.fmt.t117...\c"
printf 'chr1\t205\t210\n' > b.tmp.bed
printf '***** ERROR: writeCount option is not valid with BAM query input, unless bed output is specified with -bed option. *****\nexit 1\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -c > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
tail -n 2 obs > obs2
check obs2 exp
rm obs obs2 exp b.tmp.bed

###########################################################
#  intersect.fmt.t118: -ubam output is BAM with the same
#  records
###########################################################
echo -e "    intersect.fmt.t118...\c"
printf 'chr1\t0\t600\n' > b.tmp.bed
printf 'r1\t0\tchr1\t101\t30\t50M\t*\t0\t0\t*\t*\nr2\t16\tchr1\t201\t40\t20M100N30M\t*\t0\t0\t*\t*\nr3\t0\tchr1\t501\t50\t20M10D30M\t*\t0\t0\t*\t*\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -ubam | $htsutil viewbamrecords > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t119: BAM output carries the input header
#  @HD, both @SQ lines and the @CO line, then the hit.
###########################################################
echo -e "    intersect.fmt.t119...\c"
printf 'chr1\t120\t130\n' > b.tmp.bed
printf '@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000\n@SQ\tSN:chr2\tLN:10000\n@CO\tformats test reads\nr1\t0\tchr1\t101\t30\t50M\t*\t0\t0\t*\t*\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed | $htsutil viewbam > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t120: BAM query read from stdin given as
#  "stdin"
###########################################################
echo -e "    intersect.fmt.t120...\c"
printf 'chr1\t120\t130\n' > b.tmp.bed
printf 'r1\t0\tchr1\t101\t30\t50M\t*\t0\t0\t*\t*\n' > exp
cat reads.tmp.bam | $BT intersect -a stdin -b b.tmp.bed | $htsutil viewbamrecords > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t121: -v on a BAM query writes BAM
#  B covers all of chr1. The reads left are the unplaced
#  unmapped r5, the chr2 read r7 and r10, which is placed
#  on chr1 but unmapped and so never overlaps.
###########################################################
echo -e "    intersect.fmt.t121...\c"
printf 'chr1\t0\t10000\n' > b.tmp.bed
printf 'r5\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\nr7\t0\tchr2\t101\t0\t10M\t*\t0\t0\t*\t*\nr10\t4\tchr1\t1501\t0\t*\t*\t0\t0\t*\t*\n' > exp
$BT intersect -a reads.tmp.bam -b b.tmp.bed -v | $htsutil viewbamrecords > obs
check obs exp
rm obs exp b.tmp.bed

###########################################################
#  intersect.fmt.t122: BAM database records print as BED6
#  A 205-210 hits r2; A 705-755 hits both mates of r4.
#  Each read is chrom, start, end, name (with /1 /2),
#  MAPQ, strand.
#  Pins the current output: the two mates are reported in
#  file order.
###########################################################
echo -e "    intersect.fmt.t122...\c"
printf 'chr1\t205\t210\nchr1\t705\t755\n' > a.tmp.bed
printf 'chr1\t205\t210\tchr1\t200\t350\tr2\t40\t-\nchr1\t705\t755\tchr1\t700\t710\tr4/1\t60\t+\nchr1\t705\t755\tchr1\t750\t760\tr4/2\t60\t-\n' > exp
$BT intersect -a a.tmp.bed -b reads.tmp.bam -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t123: BAM database, -split -wo counts
#  bases shared with the read's blocks
#  A 210-330 and r2's blocks 200-220, 320-350: 10 + 10.
###########################################################
echo -e "    intersect.fmt.t123...\c"
printf 'chr1\t210\t330\n' > a.tmp.bed
printf 'chr1\t210\t330\tchr1\t200\t350\tr2\t40\t-\t20\n' > exp
$BT intersect -a a.tmp.bed -b reads.tmp.bam -split -wo > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t124: BAM database, -wo without -split
#  counts the span
#  A 210-330 lies inside r2's span: 120.
###########################################################
echo -e "    intersect.fmt.t124...\c"
printf 'chr1\t210\t330\n' > a.tmp.bed
printf 'chr1\t210\t330\tchr1\t200\t350\tr2\t40\t-\t120\n' > exp
$BT intersect -a a.tmp.bed -b reads.tmp.bam -wo > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t125: unmapped database reads are never
#  hits, even when placed
#  r10 is placed at chr1 1501 but unmapped: A 1500-1510
#  counts 0.
###########################################################
echo -e "    intersect.fmt.t125...\c"
printf 'chr1\t1500\t1510\n' > a.tmp.bed
printf 'chr1\t1500\t1510\t0\n' > exp
$BT intersect -a a.tmp.bed -b reads.tmp.bam -c > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t126: BAM database read from stdin as "-"
###########################################################
echo -e "    intersect.fmt.t126...\c"
printf 'chr1\t120\t130\n' > a.tmp.bed
printf 'chr1\t120\t130\tchr1\t100\t150\tr1\t30\t+\n' > exp
cat reads.tmp.bam | $BT intersect -a a.tmp.bed -b - -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed

###########################################################
#  intersect.fmt.t127: BAM database, -c counts mapped reads
#  chr1 0-10000 holds r1, r2, r3, both r4 mates, r6, r8
#  and r9: 8. The chr2 read and the unmapped reads do not
#  count.
###########################################################
echo -e "    intersect.fmt.t127...\c"
printf 'chr1\t0\t10000\n' > a.tmp.bed
printf 'chr1\t0\t10000\t8\n' > exp
$BT intersect -a a.tmp.bed -b reads.tmp.bam -c > obs
check obs exp
rm obs exp a.tmp.bed

rm reads.tmp.sam reads.tmp.bam mapped.tmp.sam mapped.tmp.bam

###########################################################
###########################################################
#  Compressed input and stdin
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t128: gzip-compressed database
###########################################################
echo -e "    intersect.fmt.t128...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t300\t400\tb2\n' > b.tmp.bed
gzip -c b.tmp.bed > b.tmp.bed.gz
printf 'chr1\t150\t160\tchr1\t150\t160\tb1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed.gz -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t129: bgzip-compressed database
###########################################################
echo -e "    intersect.fmt.t129...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t300\t400\tb2\n' > b.tmp.bed
$htsutil bgzfcompress b.tmp.bed b.tmp.bed.gz
printf 'chr1\t150\t160\tchr1\t150\t160\tb1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed.gz -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t130: gzip-compressed database on stdin
###########################################################
echo -e "    intersect.fmt.t130...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t300\t400\tb2\n' > b.tmp.bed
gzip -c b.tmp.bed > b.tmp.bed.gz
printf 'chr1\t150\t160\tchr1\t150\t160\tb1\n' > exp
cat b.tmp.bed.gz | $BT intersect -a a.tmp.bed -b - -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t131: a gzip file of two members (two gzip
#  files concatenated) is read to the end
#  Each member holds one record; both are reported.
###########################################################
echo -e "    intersect.fmt.t131...\c"
printf 'chr1\t150\t160\ta1\n' > a1.tmp.bed
printf 'chr1\t170\t180\ta2\n' > a2.tmp.bed
gzip -c a1.tmp.bed > a.tmp.bed.gz
gzip -c a2.tmp.bed >> a.tmp.bed.gz
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed.gz -b b.tmp.bed > obs
check obs exp
rm obs exp a1.tmp.bed a2.tmp.bed a.tmp.bed.gz b.tmp.bed

###########################################################
#  intersect.fmt.t132: a two-member gzip database is read
#  to the end
#  A overlaps one record from each member: count 2.
###########################################################
echo -e "    intersect.fmt.t132...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\n' > b1.tmp.bed
printf 'chr1\t170\t180\tb2\n' > b2.tmp.bed
gzip -c b1.tmp.bed > b.tmp.bed.gz
gzip -c b2.tmp.bed >> b.tmp.bed.gz
printf 'chr1\t100\t200\t2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed.gz -c > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t133: two concatenated bgzip files are
#  read to the end
###########################################################
echo -e "    intersect.fmt.t133...\c"
printf 'chr1\t150\t160\ta1\n' > a1.tmp.bed
printf 'chr1\t170\t180\ta2\n' > a2.tmp.bed
$htsutil bgzfcompress a1.tmp.bed a1.tmp.bed.gz
$htsutil bgzfcompress a2.tmp.bed a2.tmp.bed.gz
cat a1.tmp.bed.gz a2.tmp.bed.gz > a.tmp.bed.gz
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed.gz -b b.tmp.bed > obs
check obs exp
rm obs exp a1.tmp.bed a2.tmp.bed a1.tmp.bed.gz a2.tmp.bed.gz a.tmp.bed.gz b.tmp.bed

###########################################################
#  intersect.fmt.t134: gzip-compressed GFF query
#  Decompression does not change GFF handling: the
#  overlap with 100-200 is written 151-200.
###########################################################
echo -e "    intersect.fmt.t134...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > a.tmp.gff
gzip -c a.tmp.gff > a.tmp.gff.gz
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t200\t.\t+\t.\tID=g1\n' > exp
$BT intersect -a a.tmp.gff.gz -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff a.tmp.gff.gz b.tmp.bed

###########################################################
#  intersect.fmt.t135: gzip-compressed BED with CRLF line
#  endings
#  The carriage returns are dropped after decompression.
###########################################################
echo -e "    intersect.fmt.t135...\c"
printf 'chr1\t150\t160\ta1\r\nchr1\t170\t180\ta2\r\n' > a.tmp.bed
gzip -c a.tmp.bed > a.tmp.bed.gz
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed.gz -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed a.tmp.bed.gz b.tmp.bed

###########################################################
#  intersect.fmt.t136: a gzip query much larger than the
#  read buffers is read completely
#  3000 records (about 80 kB uncompressed) all lie inside
#  B, so -u prints the uncompressed file unchanged.
###########################################################
echo -e "    intersect.fmt.t136...\c"
awk 'BEGIN { for (i = 0; i < 3000; i++) printf "chr1\t%d\t%d\tread_%d\n", i * 10, i * 10 + 5, i }' > exp
gzip -c exp > a.tmp.bed.gz
printf 'chr1\t0\t100000\n' > b.tmp.bed
$BT intersect -a a.tmp.bed.gz -b b.tmp.bed -u > obs
check obs exp
rm obs exp a.tmp.bed.gz b.tmp.bed

###########################################################
#  intersect.fmt.t137: a gzip database much larger than the
#  read buffers is loaded completely
#  The same 3000 records as t136 as B: A counts 3000.
###########################################################
echo -e "    intersect.fmt.t137...\c"
awk 'BEGIN { for (i = 0; i < 3000; i++) printf "chr1\t%d\t%d\tread_%d\n", i * 10, i * 10 + 5, i }' | gzip -c > b.tmp.bed.gz
printf 'chr1\t0\t100000\n' > a.tmp.bed
printf 'chr1\t0\t100000\t3000\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed.gz -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t138: an empty gzip database
#  It has no records, so -loj prints A with the BED3 null.
#  Pins the current output: a database without records
#  gets the BED3 null.
###########################################################
echo -e "    intersect.fmt.t138...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
: | gzip -c > b.tmp.bed.gz
printf 'chr1\t100\t200\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed.gz -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed.gz

###########################################################
#  intersect.fmt.t139: GFF query on stdin
###########################################################
echo -e "    intersect.fmt.t139...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t200\t.\t+\t.\tID=g1\n' > exp
cat a.tmp.gff | $BT intersect -a stdin -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t140: compression is detected from the
#  content, not the file name
#  A plain-text file named .gz is read as text, and gzip
#  data in a file named .bed is decompressed.
###########################################################
echo -e "    intersect.fmt.t140...\c"
printf 'chr1\t150\t160\ta1\n' > a.tmp.bed.gz
printf 'chr1\t100\t200\n' > b.tmp.txt
gzip -c b.tmp.txt > b.tmp.bed
printf 'chr1\t150\t160\ta1\n' > exp
$BT intersect -a a.tmp.bed.gz -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed.gz b.tmp.txt b.tmp.bed

###########################################################
###########################################################
#  Header and comment lines
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t141: browser, track and # lines before
#  the first record form the header printed by -header
#  All three header lines are printed, in order, before
#  the results.
###########################################################
echo -e "    intersect.fmt.t141...\c"
printf 'browser position chr1:1-1000\ntrack name=a description="x y"\n#chrom\tstart\tend\tname\nchr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'browser position chr1:1-1000\ntrack name=a description="x y"\n#chrom\tstart\tend\tname\nchr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t142: without -header the header lines are
#  read and dropped
###########################################################
echo -e "    intersect.fmt.t142...\c"
printf 'browser position chr1:1-1000\ntrack name=a description="x y"\n#chrom\tstart\tend\tname\nchr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t143: comment and track lines after the
#  first record are skipped and are not part of the header
#  Only "#h1" is printed by -header; "# mid" and the late
#  track line are dropped.
###########################################################
echo -e "    intersect.fmt.t143...\c"
printf '#h1\nchr1\t150\t160\ta1\n# mid\ntrack name=late\nchr1\t170\t180\ta2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#h1\nchr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t144: track, browser and chrom header
#  words are recognised in any case
###########################################################
echo -e "    intersect.fmt.t144...\c"
printf 'Track name=x\nBROWSER hide all\nCHROM\tstart\tend\nchr1\t150\t160\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'Track name=x\nBROWSER hide all\nCHROM\tstart\tend\nchr1\t150\t160\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t145: lines starting with "!" or ">" and
#  "visibility" lines are header lines too
###########################################################
echo -e "    intersect.fmt.t145...\c"
printf '!note\n>seq\nvisibility=2\nchr1\t150\t160\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '!note\n>seq\nvisibility=2\nchr1\t150\t160\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t146: header lines of the database are
#  skipped and never printed
#  -header prints the query's header, and A has none.
#  Pins the current output: the two hits of a1 are
#  reported in file order.
###########################################################
echo -e "    intersect.fmt.t146...\c"
printf 'chr1\t150\t160\ta1\n' > a.tmp.bed
printf '#bhdr\ntrack name=b\nchr1\t100\t200\tb1\n# mid\nchr1\t155\t156\tb2\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\tchr1\t100\t200\tb1\nchr1\t150\t160\ta1\tchr1\t155\t156\tb2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t147: -header prints the header even when
#  no record is reported
###########################################################
echo -e "    intersect.fmt.t147...\c"
printf '#h1\nchr1\t500\t600\ta1\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#h1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t148: -header with -c
###########################################################
echo -e "    intersect.fmt.t148...\c"
printf '#h1\nchr1\t150\t160\ta1\nchr1\t500\t600\ta2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#h1\nchr1\t150\t160\ta1\t1\nchr1\t500\t600\ta2\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t149: -header with -u
###########################################################
echo -e "    intersect.fmt.t149...\c"
printf '#h1\nchr1\t150\t160\ta1\nchr1\t500\t600\ta2\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#h1\nchr1\t150\t160\ta1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -u -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t150: empty lines between records are
#  skipped
###########################################################
echo -e "    intersect.fmt.t150...\c"
printf 'chr1\t10\t20\n\nchr1\t30\t40\n\n\nchr1\t50\t60\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t10\t20\nchr1\t30\t40\nchr1\t50\t60\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -v > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t151: lines holding only spaces or a tab
#  are skipped
#  They come after the first four records, past the lines
#  the type checker reads.
###########################################################
echo -e "    intersect.fmt.t151...\c"
printf 'chr1\t10\t20\nchr1\t30\t40\nchr1\t50\t60\nchr1\t70\t80\n   \nchr1\t90\t95\n\t\nchr1\t96\t99\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t10\t20\nchr1\t30\t40\nchr1\t50\t60\nchr1\t70\t80\nchr1\t90\t95\nchr1\t96\t99\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -v > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  Line endings
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t152: CRLF query with -header
#  The carriage return is removed from the header line and
#  from each record.
###########################################################
echo -e "    intersect.fmt.t152...\c"
printf '#hdr\r\nchr1\t150\t160\ta1\t0\t+\r\nchr1\t170\t180\ta2\t0\t-\r\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf '#hdr\nchr1\t150\t160\ta1\t0\t+\nchr1\t170\t180\ta2\t0\t-\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t153: CRLF database, -wb
#  No carriage return reaches the printed B record.
#  Pins the current output: the two hits of A are
#  reported in file order.
###########################################################
echo -e "    intersect.fmt.t153...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t0\t+\r\nchr1\t170\t180\tb2\t0\t-\r\n' > b.tmp.bed
printf 'chr1\t150\t160\tchr1\t150\t160\tb1\t0\t+\nchr1\t170\t180\tchr1\t170\t180\tb2\t0\t-\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t154: a CRLF database is typed as if the
#  carriage return were absent
#  "+\r" would not be a strand; with it removed B is BED6,
#  so the -loj null is the BED6 ". -1 -1 . -1 .".
###########################################################
echo -e "    intersect.fmt.t154...\c"
printf 'chr2\t1\t2\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t0\t+\r\n' > b.tmp.bed
printf 'chr2\t1\t2\t.\t-1\t-1\t.\t-1\t.\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t155: LF and CRLF lines mixed in one file
###########################################################
echo -e "    intersect.fmt.t155...\c"
printf 'chr1\t150\t160\r\nchr1\t170\t180\nchr1\t185\t190\r\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\nchr1\t170\t180\nchr1\t185\t190\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t156: a query without a final newline
#  The last record is still read and printed with a
#  newline.
###########################################################
echo -e "    intersect.fmt.t156...\c"
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t157: a database without a final newline
#  Both records are loaded: count 2.
###########################################################
echo -e "    intersect.fmt.t157...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\nchr1\t170\t180\tb2' > b.tmp.bed
printf 'chr1\t100\t200\t2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t158: a last line ending in a carriage
#  return without a newline
###########################################################
echo -e "    intersect.fmt.t158...\c"
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\r' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\ta1\nchr1\t170\t180\ta2\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t159: CRLF GFF query
#  The attribute column is printed without the carriage
#  return.
###########################################################
echo -e "    intersect.fmt.t159...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\r\n' > a.tmp.gff
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
###########################################################
#  Empty and header-only files
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t160: empty database with -loj
#  An empty file has no record type; the null printed is
#  the BED3 ". -1 -1".
#  Pins the current output: a database without records
#  gets the BED3 null.
###########################################################
echo -e "    intersect.fmt.t160...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
: > b.tmp.bed
printf 'chr1\t100\t200\ta1\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t161: empty database with -c
###########################################################
echo -e "    intersect.fmt.t161...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
: > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t162: header-only database with -loj
#  Treated like an empty file: BED3 null.
#  Pins the current output: a database without records
#  gets the BED3 null.
###########################################################
echo -e "    intersect.fmt.t162...\c"
printf 'chr1\t100\t200\ta1\n' > a.tmp.bed
printf 'track name=b\n#only a header\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t163: a query of only blank lines, with
#  -header
#  No header and no records: empty output, exit 0.
###########################################################
echo -e "    intersect.fmt.t163...\c"
printf '\n\n\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'exit 0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -header > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t164: an empty query with -v
#  Nothing to report: empty output, exit 0.
###########################################################
echo -e "    intersect.fmt.t164...\c"
: > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'exit 0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -v > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  Long names, long and wide lines, non-ASCII text
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t165: a 303-character chromosome name
#  Both files use the same long name, so the records are
#  on the same chromosome and overlap at 15-20.
###########################################################
echo -e "    intersect.fmt.t165...\c"
awk 'BEGIN { s = "chr"; for (i = 0; i < 300; i++) s = s "c"; printf "%s\t10\t20\n", s }' > a.tmp.bed
awk 'BEGIN { s = "chr"; for (i = 0; i < 300; i++) s = s "c"; printf "%s\t15\t25\n", s }' > b.tmp.bed
awk 'BEGIN { s = "chr"; for (i = 0; i < 300; i++) s = s "c"; printf "%s\t15\t20\n", s }' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t166: a BED6+ record with 500 extra
#  columns is printed unchanged by -wa -wb
###########################################################
echo -e "    intersect.fmt.t166...\c"
awk 'BEGIN { printf "chr1\t10\t20\tn\t0\t+"; for (i = 1; i <= 500; i++) printf "\tf%d", i; printf "\n" }' > a.tmp.bed
printf 'chr1\t0\t1000\n' > b.tmp.bed
awk 'BEGIN { printf "chr1\t10\t20\tn\t0\t+"; for (i = 1; i <= 500; i++) printf "\tf%d", i; printf "\tchr1\t0\t1000\n" }' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t167: a 5000-character name after the
#  first lines of the query
#  The long line starts at byte 88 and spans five of the
#  1023-byte chunks the file is read in. Every record lies
#  inside B, so -wa prints the file unchanged.
###########################################################
echo -e "    intersect.fmt.t167...\c"
awk 'BEGIN { for (i = 0; i < 5; i++) printf "chr1\t%d\t%d\tshort%d\n", i * 10, i * 10 + 5, i; s = ""; for (i = 0; i < 5000; i++) s = s "n"; printf "chr1\t100\t105\t%s\nchr1\t110\t115\tafter\n", s }' > a.tmp.bed
printf 'chr1\t0\t1000\n' > b.tmp.bed
cp a.tmp.bed exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wa > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t168: a 5000-character name in the
#  database is printed whole by -wb
#  Only the long record (100-105) overlaps A 100-200.
###########################################################
echo -e "    intersect.fmt.t168...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
awk 'BEGIN { for (i = 0; i < 5; i++) printf "chr1\t%d\t%d\tshort%d\n", i * 10, i * 10 + 5, i; s = ""; for (i = 0; i < 5000; i++) s = s "n"; printf "chr1\t100\t105\t%s\nchr1\t300\t305\tafter\n", s }' > b.tmp.bed
awk 'BEGIN { s = ""; for (i = 0; i < 5000; i++) s = s "n"; printf "chr1\t100\t105\tchr1\t100\t105\t%s\n", s }' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t169: non-ASCII bytes in a name are passed
#  through unchanged
#  The name holds the UTF-8 bytes of "a-umlaut"; the rest
#  of the file is plain text.
###########################################################
echo -e "    intersect.fmt.t169...\c"
printf 'chr1\t150\t160\tgene_n\303\244me_with_a_long_enough_name\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t160\tgene_n\303\244me_with_a_long_enough_name\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  Chromosome naming conventions (unsorted input)
#  The database is read first, so its first record sets
#  the convention. A record that disagrees triggers a
#  warning, printed when found and again at exit; the
#  exit status stays 0.
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t170: query uses "chr1", database "1"
#  The names differ, so there is no overlap; the query
#  record is reported in the warning.
#  Pins the current output: the warning appears twice.
###########################################################
echo -e "    intersect.fmt.t170...\c"
printf 'chr1\t15\t25\n' > a.tmp.bed
printf '1\t10\t20\n' > b.tmp.bed
printf '***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t15\t25\n\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t15\t25\n\nexit 0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t171: a query mixing "chr1" and "1"
#  B uses "chr1". The chr1 record overlaps and is printed
#  (stdout first below); the "1" record is named in the
#  warning (stderr after it).
###########################################################
echo -e "    intersect.fmt.t171...\c"
printf 'chr1\t10\t20\n1\t10\t20\n' > a.tmp.bed
printf 'chr1\t15\t25\n' > b.tmp.bed
printf 'chr1\t15\t20\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\n1\t10\t20\n\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\n1\t10\t20\n\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2> obserr
cat obserr >> obs
check obs exp
rm obs obserr exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t172: -nonamecheck silences the warning
#  Same input as t171: the overlap only, no stderr.
###########################################################
echo -e "    intersect.fmt.t172...\c"
printf 'chr1\t10\t20\n1\t10\t20\n' > a.tmp.bed
printf 'chr1\t15\t25\n' > b.tmp.bed
printf 'chr1\t15\t20\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -nonamecheck > obs 2>&1
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t173: the "chr" test ignores case, but
#  names must match exactly
#  "Chr1" counts as having "chr", so there is no warning;
#  it is a different name from "chr1", so no overlap.
###########################################################
echo -e "    intersect.fmt.t173...\c"
printf 'Chr1\t10\t20\n' > a.tmp.bed
printf 'chr1\t15\t25\n' > b.tmp.bed
printf 'Chr1\t10\t20\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs 2>&1
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t174: chromosome names are compared as
#  whole strings
#  chr1_random and chr10 share a prefix with chr1 but do
#  not match it.
###########################################################
echo -e "    intersect.fmt.t174...\c"
printf 'chr1\t10\t20\nchr1_random\t10\t20\nchr10\t10\t20\n' > a.tmp.bed
printf 'chr1\t15\t25\n' > b.tmp.bed
printf 'chr1\t10\t20\t1\nchr1_random\t10\t20\t0\nchr10\t10\t20\t0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed -c > obs 2>&1
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t175: the warning is given for the first
#  disagreeing record only
#  Both query records use "chr" against a "1" database;
#  only chr1:10-20 is named.
#  Pins the current output: one warning per run.
###########################################################
echo -e "    intersect.fmt.t175...\c"
printf 'chr1\t10\t20\nchr2\t10\t20\n' > a.tmp.bed
printf '1\t15\t25\n' > b.tmp.bed
printf '***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t10\t20\n\n***** WARNING: File a.tmp.bed has inconsistent naming convention for record:\nchr1\t10\t20\n\nexit 0\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  intersect.fmt.t176: a database mixing "chr1" and "1"
#  The second database record is named in the warning,
#  and the query still overlaps the chr1 record.
###########################################################
echo -e "    intersect.fmt.t176...\c"
printf 'chr1\t10\t20\n' > a.tmp.bed
printf 'chr1\t15\t25\n1\t15\t25\n' > b.tmp.bed
printf 'chr1\t15\t20\n***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t15\t25\n\n***** WARNING: File b.tmp.bed has inconsistent naming convention for record:\n1\t15\t25\n\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.bed > obs 2> obserr
cat obserr >> obs
check obs exp
rm obs obserr exp a.tmp.bed b.tmp.bed

###########################################################
###########################################################
#  -wb and -loj with GFF and VCF queries, and databases of
#  different formats in one run
###########################################################
###########################################################

###########################################################
#  intersect.fmt.t177: GFF query with -loj
#  -loj prints the whole query record: the hit row is the
#  GFF record and its B record, the miss row the GFF
#  record and the BED3 null.
###########################################################
echo -e "    intersect.fmt.t177...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\nchr2\tsrc\tgene\t151\t250\t.\t+\t.\tID=g2\n' > a.tmp.gff
printf 'chr1\t160\t170\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\tchr1\t160\t170\nchr2\tsrc\tgene\t151\t250\t.\t+\t.\tID=g2\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t178: VCF query with -loj
#  The whole VCF record on both rows; the miss gets the
#  BED3 null.
###########################################################
echo -e "    intersect.fmt.t178...\c"
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > a.tmp.vcf
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t501\tsnp2\tC\tT\t50\tPASS\tDP=3\n' >> a.tmp.vcf
printf 'chr1\t90\t110\n' > b.tmp.bed
printf 'chr1\t101\tsnp\tA\tG\t50\tPASS\tDP=10\tchr1\t90\t110\nchr1\t501\tsnp2\tC\tT\t50\tPASS\tDP=3\t.\t-1\t-1\n' > exp
$BT intersect -a a.tmp.vcf -b b.tmp.bed -loj > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  intersect.fmt.t179: GFF query with -wb
#  The GFF record is cut to the overlap with 160-170,
#  written 161-170, and B follows.
###########################################################
echo -e "    intersect.fmt.t179...\c"
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > a.tmp.gff
printf 'chr1\t160\t170\tb1\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t161\t170\t.\t+\t.\tID=g1\tchr1\t160\t170\tb1\n' > exp
$BT intersect -a a.tmp.gff -b b.tmp.bed -wb > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  intersect.fmt.t180: databases of three formats in one
#  run, -wa -wb
#  Each hit is printed in its own format after the number
#  of its file (1 GFF, 2 VCF, 3 BED6).
#  Pins the current output: hits are reported in database
#  order.
###########################################################
echo -e "    intersect.fmt.t180...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > c.tmp.vcf
printf 'chr1\t121\tsnp\tA\tG\t50\tPASS\tDP=10\n' >> c.tmp.vcf
printf 'chr1\t190\t300\tb6\t0\t-\n' > d.tmp.bed
printf 'chr1\t100\t200\t1\tchr1\tsrc\tgene\t151\t250\t.\t+\t.\tID=g1\nchr1\t100\t200\t2\tchr1\t121\tsnp\tA\tG\t50\tPASS\tDP=10\nchr1\t100\t200\t3\tchr1\t190\t300\tb6\t0\t-\n' > exp
$BT intersect -a a.tmp.bed -b b.tmp.gff c.tmp.vcf d.tmp.bed -wa -wb > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff c.tmp.vcf d.tmp.bed

[[ $FAILURES -eq 0 ]] || exit 1;
