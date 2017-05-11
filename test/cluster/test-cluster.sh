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

# cat in.bed
# chr1	72017	884436	a	1	+
# chr1	72017	844113	b	2	+
# chr1	939517	1011278	c	3	+
# chr1	1142976	1203168	d	4	+
# chr1	1153667	1298845	e	5	-
# chr1	1153667	1219633	f	6	+
# chr1	1155173	1200334	g	7	-
# chr1	1229798	1500664	h	8	-
# chr1	1297735	1357056	i	9	+
# chr1	1844181	1931789	j	10	-

###########################################################
#  basic cluster.
###########################################################
echo -e "    cluster.t1...\c"
echo \
"chr1	72017	884436	a	1	+	1
chr1	72017	844113	b	2	+	1
chr1	939517	1011278	c	3	+	2
chr1	1142976	1203168	d	4	+	3
chr1	1153667	1298845	e	5	-	3
chr1	1153667	1219633	f	6	+	3
chr1	1155173	1200334	g	7	-	3
chr1	1229798	1500664	h	8	-	3
chr1	1297735	1357056	i	9	+	3
chr1	1844181	1931789	j	10	-	4" > exp
$BT cluster -i in.bed > obs
check obs exp
rm obs exp

###########################################################
#  stranded cluster.
###########################################################
echo -e "    cluster.t2...\c"
echo \
"chr1	72017	884436	a	1	+	1
chr1	72017	844113	b	2	+	1
chr1	939517	1011278	c	3	+	2
chr1	1142976	1203168	d	4	+	3
chr1	1153667	1219633	f	6	+	3
chr1	1297735	1357056	i	9	+	4
chr1	1153667	1298845	e	5	-	5
chr1	1155173	1200334	g	7	-	5
chr1	1229798	1500664	h	8	-	5
chr1	1844181	1931789	j	10	-	6" > exp
$BT cluster -i in.bed -s > obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;
