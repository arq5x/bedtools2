set -e;
BT=${BT-../../bin/bedtools}

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
#  Test default
############################################################
echo -e "    sort.t01...\c"
echo \
"chr3	10	220	f	12	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	210	525	d	21	+
chr7	240	560	x	86	-
chr7	2100	2310	e	32	+
chr9	110	120	a	81	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+" > exp
$BT sort -i a.bed > obs
check obs exp
rm obs exp

###########################################################
#  Test sizeA
############################################################
echo -e "    sort.t02...\c"
echo \
"chr9	110	120	a	81	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+
chr3	10	220	f	12	+
chr7	2100	2310	e	32	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	210	525	d	21	+
chr7	240	560	x	86	-" > exp
$BT sort -i a.bed -sizeA > obs
check obs exp
rm obs exp

###########################################################
#  Test sizeD
############################################################
echo -e "    sort.t03...\c"
echo \
"chr7	240	560	x	86	-
chr7	210	525	d	21	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr3	10	220	f	12	+
chr7	2100	2310	e	32	+
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+
chr9	110	120	a	81	-" > exp
$BT sort -i a.bed -sizeD > obs
check obs exp
rm obs exp


###########################################################
#  Test chrThenSizeA
############################################################
echo -e "    sort.t04...\c"
echo \
"chr3	10	220	f	12	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	2100	2310	e	32	+
chr7	210	525	d	21	+
chr7	240	560	x	86	-
chr9	110	120	a	81	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+" > exp
$BT sort -i a.bed -chrThenSizeA > obs
check obs exp
rm obs exp

###########################################################
#  Test chrThenSizeD
############################################################
echo -e "    sort.t05...\c"
echo \
"chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr3	10	220	f	12	+
chr7	240	560	x	86	-
chr7	210	525	d	21	+
chr7	2100	2310	e	32	+
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+
chr9	110	120	a	81	-" > exp
$BT sort -i a.bed -chrThenSizeD > obs
check obs exp
rm obs exp


###########################################################
#  Test chrThenScoreA
############################################################
echo -e "    sort.t06...\c"
echo \
"chr3	10	220	f	12	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	210	525	d	21	+
chr7	2100	2310	e	32	+
chr7	240	560	x	86	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+
chr9	110	120	a	81	-" > exp
$BT sort -i a.bed -chrThenScoreA > obs
check obs exp
rm obs exp

###########################################################
#  Test chrThenScoreD
############################################################
echo -e "    sort.t07...\c"
echo \
"chr3	100	320	g	96	-
chr3	40	260	p	41	-
chr3	10	220	f	12	+
chr7	240	560	x	86	-
chr7	2100	2310	e	32	+
chr7	210	525	d	21	+
chr9	110	120	a	81	-
chr9	1100	1120	b	12	+
chr9	140	160	z	05	+" > exp
$BT sort -i a.bed -chrThenScoreD > obs
check obs exp
rm obs exp

###########################################################
#  Test faidx
############################################################
echo -e "    sort.t08...\c"
echo \
"chr9	110	120	a	81	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+
chr3	10	220	f	12	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	210	525	d	21	+
chr7	240	560	x	86	-
chr7	2100	2310	e	32	+" > exp
$BT sort -i a.bed -faidx names.txt > obs
check obs exp
rm obs exp

###########################################################
#  Test header
############################################################
echo -e "    sort.t09...\c"
echo \
"#Header line for a.bed
chr3	10	220	f	12	+
chr3	40	260	p	41	-
chr3	100	320	g	96	-
chr7	210	525	d	21	+
chr7	240	560	x	86	-
chr7	2100	2310	e	32	+
chr9	110	120	a	81	-
chr9	140	160	z	05	+
chr9	1100	1120	b	12	+" > exp
$BT sort -i a.bed -header > obs
check obs exp
rm obs exp

###########################################################
#  Test zero-length interval
############################################################
echo -e "    sort.t10...\c"
echo \
"chr1	10	20	RegionA	0	+
chr1	24	30	RegionB	0	+
chr1	25	25	RegionC	0	+
chr1	30	40	RegionD	0	+" > exp
$BT sort -i b.bed > obs
check obs exp
rm obs exp

[[ $FAILURES -eq 0 ]] || exit 1;
