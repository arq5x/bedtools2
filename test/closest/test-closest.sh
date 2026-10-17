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

# cat a.bed
# chr1	20	21

# cat b.bed
# chr1	20	21

# cat b-one-bp-closer.bed
# chr1	19	21

###########################################################
# test 1bp apart; checking for off-by-one errors
###########################################################
echo -e "    closest.t1...\c"
echo \
"chr1	10	20	chr1	20	21	1" > exp
$BT closest -a a.bed -b b.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t1
###########################################################
echo -e "    closest.t2...\c"
echo \
"chr1	20	21	chr1	10	20	1" > exp
$BT closest -a b.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test 0bp apart; checking for off-by-one errors
###########################################################
echo -e "    closest.t3...\c"
echo \
"chr1	10	20	chr1	19	21	0" > exp
$BT closest -a a.bed -b b-one-bp-closer.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test reciprocal of t3
###########################################################
echo -e "    closest.t4...\c"
echo \
"chr1	19	21	chr1	10	20	0" > exp
$BT closest -a b-one-bp-closer.bed -b a.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest without forcing different names ( -N )
###########################################################
echo -e "    closest.t5...\c"
echo \
"chr1	10	20	break1	chr1	40	50	break1	21
chr1	55	58	break2	chr1	60	70	break2	3" > exp
$BT closest -a a.names.bed -b b.names.bed -d > obs
check obs exp
rm obs exp

###########################################################
# test closest with forcing different names ( -N )
###########################################################
echo -e "    closest.t6...\c"
echo \
"chr1	10	20	break1	chr1	60	70	break2	41
chr1	55	58	break2	chr1	40	50	break1	6" > exp
$BT closest -a a.names.bed -b b.names.bed -d -N > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -s yet no matching strands on chrom
###########################################################
echo -e "    closest.t7...\c"
echo \
"chr1	100	200	a	10	+	.	-1	-1	.	-1	." > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -s > obs
check obs exp
rm obs exp

###########################################################
# test closest forcing -S with only an opp strands on chrom
###########################################################
echo -e "    closest.t8...\c"
echo \
"chr1	100	200	a	10	+	chr1	90	120	b	1	-" > exp
$BT closest -a strand-test-a.bed -b strand-test-b.bed -S > obs
check obs exp
rm obs exp


###########################################################
# test reproting of all overlapping features
###########################################################
echo -e "    closest.t9...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed > obs
check obs exp
rm obs exp

###########################################################
# test reproting of first overlapping feature
###########################################################
echo -e "    closest.t10...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	150	201
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t first > obs
check obs exp
rm obs exp

###########################################################
# test reproting of last overlapping feature
###########################################################
echo -e "    closest.t11...\c"
echo \
"chr1	100	101	chr1	100	101
chr1	200	201	chr1	175	375
chr1	300	301	chr1	175	375
chr1	100000	100010	chr1	175	375
chr1	100020	100040	chr1	175	375
chr2	1	10	.	-1	-1
chr2	20	30	.	-1	-1" > exp
$BT closest -a close-a.bed -b close-b.bed -t last > obs
check obs exp
rm obs exp



###########################################################
#
# TEST MULTIPLE DATABASES
#
###########################################################


###########################################################
# test 3 dbs, -each mode, which is the default
###########################################################
echo -e "    closest.t13...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	3	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed > obs
check obs exp
rm obs exp


###########################################################
# test 3 dbs with -names option
###########################################################
echo -e "    closest.t14...\c"
echo \
"chr1	80	100	q1	1	+	a	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	b	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	c	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -names a b c > obs
check obs exp
rm obs exp


###########################################################
# test 3 dbs with -filenames option
###########################################################
echo -e "    closest.t15...\c"
echo \
"chr1	80	100	q1	1	+	mdb1.bed	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	mdb2.bed	chr1	120	170	db2.2	2	-
chr1	80	100	q1	1	+	mdb3.bed	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -filenames > obs
check obs exp
rm obs exp

###########################################################
# test 3 dbs, -all mode
###########################################################
echo -e "    closest.t16...\c"
echo \
"chr1	80	100	q1	1	+	3	chr1	70	90	d3.1	3	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed mdb3.bed -mdb all > obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all
###########################################################
echo -e "    closest.t17...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all > obs
check obs exp
rm obs exp

###########################################################
# test 2 dbs, tie mode = first
###########################################################
echo -e "    closest.t18...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -mdb all -t first > obs
check obs exp
rm obs exp

###########################################################
# test 2 dbs, tie mode = last
###########################################################
echo -e "    closest.t19...\c"
echo \
"chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -mdb all -t last > obs
check obs exp
rm obs exp

###########################################################
# test same strand 
###########################################################
echo -e "    closest.t20...\c"
echo \
"chr1	80	100	q1	1	+	chr1	5	15	d1.1	1	+" > exp
$BT closest -a mq1.bed -b mdb1.bed -s> obs
check obs exp
rm obs exp

###########################################################
# test diff strand 
###########################################################
echo -e "    closest.t21...\c"
echo \
"chr1	80	100	q1	1	+	chr1	20	60	d1.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed -S> obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all, same strand
###########################################################
echo -e "    closest.t22...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	5	15	d1.1	1	+" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all -mdb all -s > obs
check obs exp
rm obs exp


###########################################################
# test 2 dbs, tie mode = all, diff strand
###########################################################
echo -e "    closest.t23...\c"
echo \
"chr1	80	100	q1	1	+	1	chr1	20	60	d1.2	2	-
chr1	80	100	q1	1	+	2	chr1	120	170	db2.2	2	-" > exp
$BT closest -a mq1.bed -b mdb1.bed mdb2.bed -t all -mdb all -S > obs
check obs exp
rm obs exp






###########################################################
#
# TEST -D OPTION
#
###########################################################


###########################################################
# hit on left, forward query, forward hit, ref mode
###########################################################
echo -e "    closest.t24...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D ref > obs
check obs exp
rm obs exp


###########################################################
# hit on left, forward query, forward hit, a mode
###########################################################
echo -e "    closest.t25...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, forward hit, b mode
###########################################################
echo -e "    closest.t26...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q1.bed -b d_d1F.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, reverse hit, ref mode
###########################################################
echo -e "    closest.t27...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D ref > obs
check obs exp
rm obs exp


###########################################################
# hit on left, forward query, reverse hit, a mode
###########################################################
echo -e "    closest.t28...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on left, forward query, reverse hit, b mode
###########################################################
echo -e "    closest.t29...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q1.bed -b d_d1R.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, forward hit, ref mode
###########################################################
echo -e "    closest.t30...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	-21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, forward hit, a mode
###########################################################
echo -e "    closest.t31...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D a > obs
check obs exp
rm obs exp


###########################################################
# hit on left, reverse query, forward hit, b mode
###########################################################
echo -e "    closest.t32...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1F.1	10	+	21" > exp
$BT closest -a d_q2.bed -b d_d1F.bed -D b > obs
check obs exp
rm obs exp




###########################################################
# hit on left, reverse query, reverse hit, ref mode
###########################################################
echo -e "    closest.t33...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on left, reverse query, reverse hit, a mode
###########################################################
echo -e "    closest.t34...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D a > obs
check obs exp
rm obs exp


###########################################################
# hit on left, reverse query, reverse hit, b mode
###########################################################
echo -e "    closest.t35...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	40	60	d1R.1	10	-	-21" > exp
$BT closest -a d_q2.bed -b d_d1R.bed -D b > obs
check obs exp
rm obs exp



###########################################################
# hit on right, forward query, forward hit, ref mode
###########################################################
echo -e "    closest.t36...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, forward hit, a mode
###########################################################
echo -e "    closest.t37...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, forward hit, b mode
###########################################################
echo -e "    closest.t38...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q1.bed -b d_d2F.bed -D b > obs
check obs exp
rm obs exp


###########################################################
# hit on right, forward query, reverse hit, ref mode
###########################################################
echo -e "    closest.t39...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, reverse hit, a mode
###########################################################
echo -e "    closest.t40...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, forward query, reverse hit, b mode
###########################################################
echo -e "    closest.t41...\c"
echo \
"chr1	80	100	d_q1.1	5	+	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q1.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, ref mode
###########################################################
echo -e "    closest.t42...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, a mode
###########################################################
echo -e "    closest.t43...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, forward hit, b mode
###########################################################
echo -e "    closest.t44...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2F.1	10	+	-41" > exp
$BT closest -a d_q2.bed -b d_d2F.bed -D b > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, ref mode
###########################################################
echo -e "    closest.t45...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D ref > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, a mode
###########################################################
echo -e "    closest.t46...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	-41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D a > obs
check obs exp
rm obs exp

###########################################################
# hit on right, reverse query, reverse hit, b mode
###########################################################
echo -e "    closest.t47...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported 
############################################################
echo -e "    closest.t48...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+
chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a a2.bed -b b2.bed > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported, but with -s
############################################################
echo -e "    closest.t49...\c"
echo \
"chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a a2.bed -b b2.bed -s > obs
check obs exp
rm obs exp

############################################################
# Make sure non-overlapping ties are reported, but with -S
############################################################
echo -e "    closest.t50...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+" > exp
$BT closest -a a2.bed -b b2.bed -S > obs
check obs exp
rm obs exp


echo -e "    closest.t51...\c"
echo \
"chr1	80	100	d_q2.1	5	-	chr1	140	160	d2R.1	10	-	41" > exp
$BT closest -a d_q2.bed -b d_d2R.bed -D b > obs
check obs exp
rm obs exp


###########################################################
# see that stranded sweep doesn't prematurely purge 
# records from the cache
###########################################################
echo -e "    closest.t52...\c"
echo \
"chr1	80	100	c1	20	+	chr1	50	60	d2	20	+
chr1	110	130	c2	20	-	chr1	20	40	d1	20	-" > exp
$BT closest -a strand-test-c.bed -b strand-test-d.bed -s > obs
check obs exp
rm obs exp

###########################################################
# check oppostite stranded sweep 
###########################################################
echo -e "    closest.t53...\c"
echo \
"chr1	80	100	c1	20	+	chr1	20	40	d1	20	-
chr1	110	130	c2	20	-	chr1	50	60	d2	20	+" > exp
$BT closest -a strand-test-c.bed -b strand-test-d.bed -S > obs
check exp obs
rm obs exp


###########################################################
# check -iu
###########################################################
echo -e "    closest.t54...\c"
echo \
"chr1	100	120	chr1	200	210	81" > exp
$BT closest -a d.bed -b d_iu.bed -D ref -iu > obs
check exp obs
rm exp obs


##########################################################
# check -id
###########################################################
echo -e "    closest.t55...\c"
echo \
"chr1	100	120	chr1	10	20	-81" > exp
$BT closest -a d.bed -b d_id.bed -D ref -id > obs
check exp obs
rm exp obs

##########################################################
# check ties, single db
###########################################################
echo -e "    closest.t56...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+
chr1	10	20	a1	1	-	chr1	21	22	b2	1	-" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed > obs
check exp obs
rm exp obs


##########################################################
# check ties, single db, -iu
###########################################################
echo -e "    closest.t57...\c"
echo \
"chr1	10	20	a1	1	-	chr1	21	22	b2	1	-	2" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed -D ref -iu > obs
check exp obs
rm exp obs

##########################################################
# check ties, single db, -id
###########################################################
echo -e "    closest.t58...\c"
echo \
"chr1	10	20	a1	1	-	chr1	8	9	b1	1	+	-2" > exp
$BT closest -a bug157_a.bed -b bug157_b.bed -D ref -id > obs 
check exp obs
rm exp obs

##########################################################
# check -header
###########################################################
echo -e "    closest.t59...\c"
echo \
"#Header for file a.bed
chr1	10	20	chr1	20	21" > exp
$BT closest -a a.bed -b b.bed -header > obs
check exp obs
rm exp obs


###########################################################
#
# BUG 281 TESTS: -s and -S (correct cache purging for 
#   opposite strand sweep).
#
###########################################################

###########################################################
# a_med vs b, -s
###########################################################
echo -e "    closest.t60...\c"
echo \
"chr1	249120154	249120155	21176	0	-	chr1	247242115	247242116	551	0	-
chr1	249132529	249132530	6425	0	+	chr1	247495324	247495325	510	0	+" > exp
$BT closest -a bug281_a.medium.bed -b bug281_b.bed -s > obs
check exp obs
rm exp obs


###########################################################
# a_med vs b, -S
###########################################################
echo -e "    closest.t61...\c"
echo \
"chr1	249120154	249120155	21176	0	-	chr1	247495324	247495325	510	0	+
chr1	249132529	249132530	6425	0	+	chr1	247242115	247242116	551	0	-" > exp
$BT closest -a bug281_a.medium.bed -b bug281_b.bed -S > obs
check exp obs
rm exp obs

###########################################################
# a_med.flip vs b, -s
###########################################################
echo -e "    closest.t62...\c"
echo \
"chr1	249120154	249120155	21176	0	+	chr1	247495324	247495325	510	0	+
chr1	249132529	249132530	6425	0	-	chr1	247242115	247242116	551	0	-" > exp
$BT closest -a bug281_a.flip.medium.bed -b bug281_b.bed -s > obs
check exp obs
rm exp obs

###########################################################
# a_med.flip vs b, -S
###########################################################
echo -e "    closest.t63...\c"
echo \
"chr1	249120154	249120155	21176	0	+	chr1	247242115	247242116	551	0	-
chr1	249132529	249132530	6425	0	-	chr1	247495324	247495325	510	0	+" > exp
$BT closest -a bug281_a.flip.medium.bed -b bug281_b.bed -S > obs
check exp obs
rm exp obs

###########################################################
# a_med vs b.flip, -s
###########################################################
echo -e "    closest.t64...\c"
echo \
"chr1	249120154	249120155	21176	0	-	chr1	247495324	247495325	510	0	-
chr1	249132529	249132530	6425	0	+	chr1	247242115	247242116	551	0	+" > exp
$BT closest -a bug281_a.medium.bed -b bug281_b.flip.bed -s > obs
check exp obs
rm exp obs

###########################################################
# a_med vs b.flip, -S
###########################################################
echo -e "    closest.t65...\c"
echo \
"chr1	249120154	249120155	21176	0	-	chr1	247242115	247242116	551	0	+
chr1	249132529	249132530	6425	0	+	chr1	247495324	247495325	510	0	-" > exp
$BT closest -a bug281_a.medium.bed -b bug281_b.flip.bed -S > obs
check exp obs
rm exp obs

###########################################################
# a_med.flip vs b.flip, -s
###########################################################
echo -e "    closest.t66...\c"
echo \
"chr1	249120154	249120155	21176	0	+	chr1	247242115	247242116	551	0	+
chr1	249132529	249132530	6425	0	-	chr1	247495324	247495325	510	0	-" > exp
$BT closest -a bug281_a.flip.medium.bed -b bug281_b.flip.bed -s > obs
check exp obs
rm exp obs

###########################################################
# a_med.flip vs b.flip, -S
###########################################################
echo -e "    closest.t67...\c"
echo \
"chr1	249120154	249120155	21176	0	+	chr1	247495324	247495325	510	0	-
chr1	249132529	249132530	6425	0	-	chr1	247242115	247242116	551	0	+" > exp
$BT closest -a bug281_a.flip.medium.bed -b bug281_b.flip.bed -S > obs
check exp obs
rm exp obs

###########################################################
#  Test intersect -wao with multiple databases and -names
############################################################
echo -e "    closest.t68...\c"
echo \
"1	100	200	a1	ax	b	1	100	200	b1	bx
1	100	200	a1	ax	c	1	500	600	c4	cq
1	300	400	a2	ay	b	1	100	200	b1	bx
1	300	400	a2	ay	c	1	500	600	c4	cq
1	400	500	a3	az	b	1	100	200	b1	bx
1	400	500	a3	az	c	1	500	600	c4	cq
2	500	600	a4	aq	.	.	-1	-1	.	." > exp
$BT closest -a null_a.bed -b null_b.bed null_c.bed -names b c > obs
check exp obs
rm exp obs

###########################################################
#  Test -mdb all
############################################################
echo -e "    closest.t69...\c"
echo \
"chr10	21805001	21805500	shores.bed	chr10	21803197	21805197	0	0
chr10	21805001	21805500	tfbs.bed	chr10	21805031	21805041	V	931	-	0	0
chr10	21805001	21805500	islands.bed	chr10	21805197	21805759	CpG:_22	0	0" > exp
$BT closest -a dmr.bed -b islands.bed tfbs.bed shores.bed -filenames -d -mdb all > obs
check exp obs
rm exp obs

###########################################################
#  Test -mdb all with the islands record having no overlap
############################################################
echo -e "    closest.t70...\c"
echo \
"chr10	21805001	21805500	shores.bed	chr10	21803197	21805197	0	0
chr10	21805001	21805500	tfbs.bed	chr10	21805031	21805041	V	931	-	0	0" > exp
$BT closest -a dmr.bed -b islands.2.bed tfbs.bed shores.bed -filenames -d -mdb all > obs
check exp obs
rm exp obs

###########################################################
#  Test -mdb each with the islands record having no overlap
############################################################
echo -e "    closest.t71...\c"
echo \
"chr10	21805001	21805500	islands.2.bed	chr10	21805597	21805759	CpG:_22	0	98
chr10	21805001	21805500	tfbs.bed	chr10	21805031	21805041	V	931	-	0	0
chr10	21805001	21805500	shores.bed	chr10	21803197	21805197	0	0" > exp
$BT closest -a dmr.bed -b islands.2.bed tfbs.bed shores.bed -filenames -d -mdb each > obs
check exp obs
rm exp obs

###########################################################
#
# NEW TESTS: distances, ties, options, record types, errors
#
###########################################################

###########################################################
#  closest.new.t001: -d for a B to the right of A
#  A chr1:100-200, B chr1:250-300. closest reports
#  B.start - A.end + 1 = 250 - 200 + 1 = 51 for a B on the right
#  (book-ended features are 1 apart, see closest.t1).
###########################################################
echo -e "    closest.new.t001...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t300\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t250\t300\t51\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t002: -d for a B to the left of A
#  A chr1:100-200, B chr1:20-50: A.start - B.end + 1 = 100 - 50 + 1
#  = 51.
###########################################################
echo -e "    closest.new.t002...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t50\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t20\t50\t51\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t003: -d is 0 for an overlapping B, which wins
#  over a book-ended B
#  A chr1:100-200; B chr1:150-160 overlaps it (distance 0) and
#  chr1:200-210 is book-ended (distance 1), so only the overlap
#  is reported.
###########################################################
echo -e "    closest.new.t003...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\nchr1\t200\t210\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t006: the nearer side wins, for each A in turn
#  A1 chr1:100-200: left B chr1:60-90 is 100-90+1 = 11 away, right
#  B chr1:220-250 is 220-200+1 = 21 away -> left B.
#  A2 chr1:1000-1100: left B chr1:950-990 is 11 away, right B
#  chr1:1105-1200 is 1105-1100+1 = 6 away -> right B.
###########################################################
echo -e "    closest.new.t006...\c"
printf 'chr1\t100\t200\nchr1\t1000\t1100\n' > a.tmp.bed
printf 'chr1\t60\t90\nchr1\t220\t250\nchr1\t950\t990\nchr1\t1105\t1200\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t60\t90\t11\nchr1\t1000\t1100\tchr1\t1105\t1200\t6\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t007: -t first breaks a left/right tie with the
#  left B
#  A chr1:100-200; B chr1:50-90 (100-90+1 = 11) and chr1:210-250
#  (210-200+1 = 11) tie. The left one comes first in the B file.
###########################################################
echo -e "    closest.new.t007...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t90\nchr1\t210\t250\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t50\t90\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -t first > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t008: -t last breaks the same tie with the right B
#  Same input as t007; the right B is last in the B file.
###########################################################
echo -e "    closest.new.t008...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t90\nchr1\t210\t250\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t210\t250\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -t last > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t009: ties on the same side are all reported, in
#  B file order
#  A chr1:100-200; B chr1:20-90 and chr1:40-90 both end at 90, so
#  both are 100-90+1 = 11 away; chr1:300-400 is 101 away. pins the
#  current output: tied hits print in B file order.
###########################################################
echo -e "    closest.new.t009...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t90\nchr1\t40\t90\nchr1\t300\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t20\t90\t11\nchr1\t100\t200\tchr1\t40\t90\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t010: -t first on a same-side tie keeps the first
#  B in the file
#  Same input as t009: chr1:20-90 is first.
###########################################################
echo -e "    closest.new.t010...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t90\nchr1\t40\t90\nchr1\t300\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t20\t90\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -t first > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t011: -t last on a same-side tie keeps the last
#  tied B in the file
#  Same input as t009: chr1:40-90 is the last of the two tied B.
###########################################################
echo -e "    closest.new.t011...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t90\nchr1\t40\t90\nchr1\t300\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t40\t90\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -t last > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t012: -io skips an overlapping B
#  A chr1:100-200; B chr1:150-160 overlaps and is ignored; of the
#  rest, chr1:230-240 is 230-200+1 = 31 away and chr1:40-60 is
#  100-60+1 = 41 away.
###########################################################
echo -e "    closest.new.t012...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t40\t60\nchr1\t150\t160\nchr1\t230\t240\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t230\t240\t31\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -io > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t013: -io keeps a book-ended B
#  A chr1:100-200; chr1:150-160 overlaps (ignored), chr1:200-210
#  touches A without sharing a base, so it is kept at distance 1.
#  pins the current output: -io removes only B that overlap A, so a
#  book-ended B is kept; the usage text says "close, yet not
#  touching features only" (bugs.md #19).
###########################################################
echo -e "    closest.new.t013...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\nchr1\t200\t210\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t200\t210\t1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -io > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t014: -io when every B on the chromosome overlaps A
#  Nothing is left, so the BED3 null record and distance -1 are
#  printed. pins the current output: BED3 null record ". -1 -1";
#  the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t014...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t90\t110\nchr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\t.\t-1\t-1\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -io > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t015: null record for a BED4 database
#  A on chr2 has no B on its chromosome; a BED4 null is
#  ". -1 -1 ." (chrom, start, end, name). pins the current output:
#  BED4 null record format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t015...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\tb1\nchr2\t100\t200\t.\t-1\t-1\t.\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t016: null record for a BED5 database
#  BED5 null is ". -1 -1 . -1" (name ".", score -1). pins the
#  current output: BED5 null record format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t016...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t7\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\tb1\t7\nchr2\t100\t200\t.\t-1\t-1\t.\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t017: null record for a BED6 database with -d
#  BED6 null is ". -1 -1 . -1 ." and -d appends distance -1. pins
#  the current output: BED6 null record format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t017...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t7\t-\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\tb1\t7\t-\t0\nchr2\t100\t200\t.\t-1\t-1\t.\t-1\t.\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t018: BED12 database record and its null record
#  The hit prints all 12 fields; the null is the BED6 null plus
#  six "." fields. pins the current output: BED12 null record
#  format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t018...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t180\tb1\t0\t+\t150\t180\t0\t2\t10,10,\t0,20,\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t180\tb1\t0\t+\t150\t180\t0\t2\t10,10,\t0,20,\nchr2\t100\t200\t.\t-1\t-1\t.\t-1\t.\t.\t.\t.\t.\t.\t.\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t019: BED6+2 database record and its null record
#  The hit prints all 8 fields. The null is built as BED3 plus the
#  extra columns: ". -1 -1" followed by five "." fields. pins the
#  current output: BED6+ null record format, whose score column is
#  ".", unlike the -1 of a BED6 null; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t019...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\tb1\t0\t+\tx1\ty1\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\tb1\t0\t+\tx1\ty1\nchr2\t100\t200\t.\t-1\t-1\t.\t.\t.\t.\t.\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t020: bedGraph database record and its null record
#  A 4-column file with a numeric 4th field is bedGraph; its null
#  is ". -1 -1 .". pins the current output: bedGraph null record
#  format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t020...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\t2.5\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\t2.5\nchr2\t100\t200\t.\t-1\t-1\t.\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t021: GFF database record and its null record
#  GFF start 201 is 0-based 200, so chr1:200-210 is book-ended to
#  A chr1:100-200: distance 200 - 200 + 1 = 1 (reading 201 as
#  0-based would give 2). The GFF line is printed as read; the null
#  for a 9-column GFF is ". . . -1 -1 . . . .". pins the current
#  output: GFF null record format; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t021...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t201\t210\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t100\t200\tchr1\tsrc\tgene\t201\t210\t.\t+\t.\tID=g1\t1\nchr2\t100\t200\t.\t.\t.\t-1\t-1\t.\t.\t.\t.\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.gff -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  closest.new.t022: GFF B to the right of A uses the 0-based start
#  GFF start 251 is 0-based 250, so the distance from A
#  chr1:100-200 is 250 - 200 + 1 = 51.
###########################################################
echo -e "    closest.new.t022...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\tsrc\tgene\t251\t300\t.\t+\t.\tID=g1\n' > b.tmp.gff
printf 'chr1\t100\t200\tchr1\tsrc\tgene\t251\t300\t.\t+\t.\tID=g1\t51\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.gff -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.gff

###########################################################
#  closest.new.t023: VCF database record and its null record
#  POS 201 with REF A is 0-based chr1:200-201, book-ended to A
#  chr1:100-200: distance 1 (reading POS as 0-based would give 2).
#  The VCF line is printed as read; the null for an 8-column VCF is
#  ". -1 ." followed by five "." fields. pins the current output:
#  VCF null record format; the usage notes say "none" for the
#  chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t023...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t201\trs1\tA\tG\t50\tPASS\t.\n' > b.tmp.vcf
printf 'chr1\t100\t200\tchr1\t201\trs1\tA\tG\t50\tPASS\t.\t1\nchr2\t100\t200\t.\t-1\t.\t.\t.\t.\t.\t.\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.vcf -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.vcf

###########################################################
#  closest.new.t024: GFF query records print as read and use
#  0-based coordinates for distance
#  GFF chr1 101-200 is 0-based chr1:100-200; B chr1:20-50 on the
#  left is 100 - 50 + 1 = 51 away (reading 101 as 0-based would
#  give 52).
###########################################################
echo -e "    closest.new.t024...\c"
printf 'chr1\tsrc\tgene\t101\t200\t.\t-\t.\tID=a1\n' > a.tmp.gff
printf 'chr1\t20\t50\n' > b.tmp.bed
printf 'chr1\tsrc\tgene\t101\t200\t.\t-\t.\tID=a1\tchr1\t20\t50\t51\n' > exp
$BT closest -a a.tmp.gff -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.gff b.tmp.bed

###########################################################
#  closest.new.t025: VCF query records print as read and use
#  0-based coordinates for distance
#  POS 101 with REF AC is 0-based chr1:100-102; B chr1:150-160 is
#  150 - 102 + 1 = 49 away. The VCF header is not printed.
###########################################################
echo -e "    closest.new.t025...\c"
printf '##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t101\trs1\tAC\tA\t50\tPASS\tDP=3\n' > a.tmp.vcf
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t101\trs1\tAC\tA\t50\tPASS\tDP=3\tchr1\t150\t160\t49\n' > exp
$BT closest -a a.tmp.vcf -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.vcf b.tmp.bed

###########################################################
#  closest.new.t026: BED12 query uses its full span, not its blocks
#  A chr1:100-200 has blocks 100-110 and 190-200; B chr1:150-160
#  lies in the gap but inside the span, so the distance is 0.
#  All 12 fields of A are printed.
###########################################################
echo -e "    closest.new.t026...\c"
printf 'chr1\t100\t200\ta1\t0\t+\t100\t200\t0\t2\t10,10,\t0,90,\n' > a.tmp.bed
printf 'chr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\t100\t200\t0\t2\t10,10,\t0,90,\tchr1\t150\t160\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t027: BED6+3 query and database records print every
#  column
#  A chr1:100-200 and B chr1:210-220 are 210-200+1 = 11 apart.
###########################################################
echo -e "    closest.new.t027...\c"
printf 'chr1\t100\t200\ta1\t5\t+\tp\tq\tr\n' > a.tmp.bed
printf 'chr1\t210\t220\tb1\t6\t-\tu\tv\tw\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t5\t+\tp\tq\tr\tchr1\t210\t220\tb1\t6\t-\tu\tv\tw\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t028: a zero-length query inside B
#  A chr1:150-150 is zero-length; it lies inside B chr1:100-200,
#  so the distance is 0 and A is printed with its original
#  coordinates.
###########################################################
echo -e "    closest.new.t028...\c"
printf 'chr1\t150\t150\n' > a.tmp.bed
printf 'chr1\t100\t200\n' > b.tmp.bed
printf 'chr1\t150\t150\tchr1\t100\t200\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t029: a zero-length database record inside A
#  B chr1:150-150 lies inside A chr1:100-200: distance 0, and B is
#  printed with its original coordinates.
###########################################################
echo -e "    closest.new.t029...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t150\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t150\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t030: A on a chromosome missing from B, between two
#  that B has
#  chr1 A finds chr1:300-400 (300-200+1 = 101), chr2 A gets the
#  null record, and chr3 A finds chr3:50-60 (100-60+1 = 41), not
#  a record from another chromosome. pins the current output: BED3
#  null record ". -1 -1"; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t030...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\nchr3\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t400\nchr3\t50\t60\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t300\t400\t101\nchr2\t100\t200\t.\t-1\t-1\t-1\nchr3\t100\t200\tchr3\t50\t60\t41\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t031: A on a first chromosome that B lacks
#  B starts on chr2, so the chr1 A gets the null record and the
#  chr2 A finds chr2:300-310 (300-200+1 = 101). pins the current
#  output: BED3 null record ". -1 -1"; the usage notes say "none" for
#  the chromosome and -1 for all other fields (bugs.md #9).
###########################################################
echo -e "    closest.new.t031...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr2\t300\t310\n' > b.tmp.bed
printf 'chr1\t100\t200\t.\t-1\t-1\t-1\nchr2\t100\t200\tchr2\t300\t310\t101\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t032: B chromosomes absent from A are skipped
#  B has chr1 and chr3 records before and after A's chr2; the chr2
#  A is reported against chr2:500-600 only (500-200+1 = 301).
###########################################################
echo -e "    closest.new.t032...\c"
printf 'chr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t160\nchr2\t500\t600\nchr3\t150\t160\n' > b.tmp.bed
printf 'chr2\t100\t200\tchr2\t500\t600\t301\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t033: several A records share one B and the cache
#  keeps it
#  A1 chr1:100-200 -> B chr1:300-310 at 300-200+1 = 101.
#  A2 chr1:150-250 -> same B at 300-250+1 = 51.
#  A3 chr1:400-500: left B chr1:300-310 at 400-310+1 = 91 beats
#  right B chr1:600-610 at 600-500+1 = 101.
###########################################################
echo -e "    closest.new.t033...\c"
printf 'chr1\t100\t200\nchr1\t150\t250\nchr1\t400\t500\n' > a.tmp.bed
printf 'chr1\t300\t310\nchr1\t600\t610\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t300\t310\t101\nchr1\t150\t250\tchr1\t300\t310\t51\nchr1\t400\t500\tchr1\t300\t310\t91\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t035: -mdb all keeps the sign of a -D distance
#  db x: chr1:50-60 is left of A chr1:100-200, 100-60+1 = 41,
#  upstream under -D ref, so -41. db y: chr1:300-310 is 101 away.
#  The closest across both databases is the db x record.
###########################################################
echo -e "    closest.new.t035...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\t0\t+\n' > b1.tmp.bed
printf 'chr1\t300\t310\tc1\t0\t+\n' > b2.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tx\tchr1\t50\t60\tb1\t0\t+\t-41\n' > exp
$BT closest -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -names x y -mdb all -D ref > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  closest.new.t036: -mdb all with -k 2 orders hits by distance
#  across databases
#  db x: chr1:250-260 is 51 away. db y: chr1:150-160 overlaps (0)
#  and chr1:400-410 is 201 away. The two closest are y's overlap,
#  then x's record. pins the current output: lines are ordered by
#  distance, not by database.
###########################################################
echo -e "    closest.new.t036...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t260\n' > b1.tmp.bed
printf 'chr1\t150\t160\nchr1\t400\t410\n' > b2.tmp.bed
printf 'chr1\t100\t200\ty\tchr1\t150\t160\t0\nchr1\t100\t200\tx\tchr1\t250\t260\t51\n' > exp
$BT closest -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -names x y -mdb all -k 2 -d > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  closest.new.t037: -f excludes a B whose overlap is too small,
#  and the next nearest B is reported
#  A chr1:100-200 (100 bp); B chr1:150-250 overlaps 50 bp = 0.50 of
#  A, below -f 0.6, so it is neither an overlap nor a neighbour.
#  B chr1:300-400 is 300-200+1 = 101 away.
###########################################################
echo -e "    closest.new.t037...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\nchr1\t300\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t300\t400\t101\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -f 0.6 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t038: -f keeps an overlap of exactly the fraction
#  Same input as t037; 50/100 = 0.50 >= -f 0.5, so chr1:150-250 is
#  an overlap at distance 0.
###########################################################
echo -e "    closest.new.t038...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t250\nchr1\t300\t400\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t250\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -f 0.5 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t039: -F tests the fraction of B
#  A chr1:100-200; B chr1:150-350 (200 bp) overlaps 50 bp = 0.25
#  of B, below -F 0.3, so the next B chr1:400-410 is reported at
#  400-200+1 = 201.
###########################################################
echo -e "    closest.new.t039...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t350\nchr1\t400\t410\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t400\t410\t201\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -F 0.3 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t040: -f with -r requires the fraction of A and B
#  Same input as t039: 50 bp is 0.50 of A but 0.25 of B, so -f 0.5
#  -r rejects chr1:150-350 and reports chr1:400-410 at 201.
###########################################################
echo -e "    closest.new.t040...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t350\nchr1\t400\t410\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t400\t410\t201\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -f 0.5 -r > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t041: -e accepts either fraction
#  Same input as t039: 0.50 of A fails -f 0.9 but 0.25 of B passes
#  -F 0.2, so with -e chr1:150-350 is an overlap at distance 0.
###########################################################
echo -e "    closest.new.t041...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t150\t350\nchr1\t400\t410\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t350\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -f 0.9 -F 0.2 -e > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t042: -s skips B records on the other strand and
#  with an unknown strand
#  A is on +. b0 (-) is 100-20+1 = 81 away, b1 (.) is 210-200+1 = 11
#  away, b2 (+) is 300-200+1 = 101 away. Only b2 has the same
#  known strand.
###########################################################
echo -e "    closest.new.t042...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t10\t20\tb0\t0\t-\nchr1\t210\t220\tb1\t0\t.\nchr1\t300\t310\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tchr1\t300\t310\tb2\t0\t+\t101\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -s > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t043: -S keeps only B records on the opposite known
#  strand
#  Same input as t042: only b0 (-) is opposite to A (+); it is 81
#  away.
###########################################################
echo -e "    closest.new.t043...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t10\t20\tb0\t0\t-\nchr1\t210\t220\tb1\t0\t.\nchr1\t300\t310\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tchr1\t10\t20\tb0\t0\t-\t81\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -S > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t044: -k 3 lists the overlap first, then the
#  neighbours by increasing distance
#  A chr1:100-200; chr1:150-160 overlaps (0), chr1:210-220 is 11
#  away, chr1:20-30 is 100-30+1 = 71 away, chr1:300-310 is 101 away.
#  pins the current output: the lines follow distance, not B file
#  order.
###########################################################
echo -e "    closest.new.t044...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t30\nchr1\t150\t160\nchr1\t210\t220\nchr1\t300\t310\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\t0\nchr1\t100\t200\tchr1\t210\t220\t11\nchr1\t100\t200\tchr1\t20\t30\t71\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 3 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t045: -k larger than the B records on the chromosome
#  reports what there is and warns on stderr
#  A chr1:100-200 has two B: chr1:10-20 (100-20+1 = 81) and
#  chr1:300-310 (101). -k 3 prints both on stdout and the warning
#  on stderr; the test appends stderr after stdout.
###########################################################
echo -e "    closest.new.t045...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t10\t20\nchr1\t300\t310\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t10\t20\t81\nchr1\t100\t200\tchr1\t300\t310\t101\n' > exp
printf 'Warning: Fewer hits (2) found on chr1 than requested (3). It is likely that there are fewer total records on that chromosome than requested.\n' >> exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 3 > obs 2> err.tmp
cat err.tmp >> obs
check obs exp
rm obs exp err.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t046: -k 2 still reports every tie at the last
#  distance taken
#  A chr1:100-200; chr1:210-220 is 11 away; chr1:50-80 (left) and
#  chr1:220-230 (right) are both 21 away. The 2nd closest distance
#  is tied, so both are reported: 3 lines. pins the current output:
#  within a tie the upstream (left) B comes before the downstream
#  one.
###########################################################
echo -e "    closest.new.t046...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t80\nchr1\t210\t220\nchr1\t220\t230\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t210\t220\t11\nchr1\t100\t200\tchr1\t50\t80\t21\nchr1\t100\t200\tchr1\t220\t230\t21\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 2 > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t047: -k 2 -t first takes one record from the tie
#  Same input as t046; the tie at 21 contributes only the first
#  (left) B.
###########################################################
echo -e "    closest.new.t047...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t80\nchr1\t210\t220\nchr1\t220\t230\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t210\t220\t11\nchr1\t100\t200\tchr1\t50\t80\t21\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 2 -t first > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t048: -k 2 -t last takes the last record of the tie
#  Same input as t046; the tie at 21 contributes the right B.
###########################################################
echo -e "    closest.new.t048...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t80\nchr1\t210\t220\nchr1\t220\t230\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t210\t220\t11\nchr1\t100\t200\tchr1\t220\t230\t21\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 2 -t last > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t049: -io with -k 2 skips the overlap and takes the
#  two nearest neighbours
#  Same B as t044: without chr1:150-160, the two nearest are
#  chr1:210-220 (11) and chr1:20-30 (71).
###########################################################
echo -e "    closest.new.t049...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t20\t30\nchr1\t150\t160\nchr1\t210\t220\nchr1\t300\t310\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t210\t220\t11\nchr1\t100\t200\tchr1\t20\t30\t71\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -k 2 -io > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t050: -D ref -fu reports the upstream B even when a
#  downstream B is nearer
#  A chr1:100-200 (+); left B chr1:50-60 is 41 away (upstream, so
#  -41), right B chr1:210-220 is 11 away.
###########################################################
echo -e "    closest.new.t050...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\t0\t+\nchr1\t210\t220\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tchr1\t50\t60\tb1\t0\t+\t-41\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref -fu > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t051: -D ref -fd reports the downstream B even when
#  an upstream B is nearer
#  A chr1:100-200; left B chr1:90-95 is 100-95+1 = 6 away, right B
#  chr1:300-310 is 101 away (downstream, positive).
###########################################################
echo -e "    closest.new.t051...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t90\t95\tb1\t0\t+\nchr1\t300\t310\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tchr1\t300\t310\tb2\t0\t+\t101\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref -fd > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t052: -D a -iu on a minus-strand A ignores B on the
#  right
#  A chr1:100-200 (-). Under -D a, upstream of a minus-strand A is
#  the higher coordinates: right B chr1:210-220 (11 away) is
#  upstream and ignored; left B chr1:50-60 is downstream, +41.
###########################################################
echo -e "    closest.new.t052...\c"
printf 'chr1\t100\t200\ta1\t0\t-\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\t0\t+\nchr1\t210\t220\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t-\tchr1\t50\t60\tb1\t0\t+\t41\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D a -iu > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t053: -D a -id on a minus-strand A ignores B on the
#  left
#  Same input as t052: the left B is downstream and ignored; the
#  right B is upstream, so -11.
###########################################################
echo -e "    closest.new.t053...\c"
printf 'chr1\t100\t200\ta1\t0\t-\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\t0\t+\nchr1\t210\t220\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t-\tchr1\t210\t220\tb2\t0\t+\t-11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D a -id > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t054: -D b -iu ignores a plus-strand B to the right
#  of A
#  Under -D b, A to the left of a + strand B is upstream of it:
#  right B chr1:210-220 (+, 11 away) is upstream and ignored. The
#  left + strand B chr1:50-60 has A downstream of it: +41.
###########################################################
echo -e "    closest.new.t054...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t50\t60\tb1\t0\t+\nchr1\t210\t220\tb2\t0\t+\n' > b.tmp.bed
printf 'chr1\t100\t200\ta1\t0\t+\tchr1\t50\t60\tb1\t0\t+\t41\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D b -iu > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t055: -D ref reports an overlap as 0
#  A chr1:100-200 and B chr1:150-160 overlap; the signed distance is
#  0 even though a left B chr1:50-60 exists.
###########################################################
echo -e "    closest.new.t055...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t50\t60\nchr1\t150\t160\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t150\t160\t0\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t056: -iu with -id is an error
###########################################################
echo -e "    closest.new.t056...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: Request either -iu OR -id, not both.\n*****\n\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref -iu -id > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t057: -iu without -D is an error
###########################################################
echo -e "    closest.new.t057...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: When requesting -iu or -id, you also need to specify -D.\n*****\n\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -iu > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t058: -fd without -D is an error
###########################################################
echo -e "    closest.new.t058...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: When requesting -fu or -fd, you also need to specify -D.\n*****\n\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -fd > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t059: -iu with -fu is an error
###########################################################
echo -e "    closest.new.t059...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf "*****ERROR: Can't both ignore upstream and force upstream.\n*****\n\nexit 1\n" > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref -iu -fu > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t060: -id with -fd is an error
###########################################################
echo -e "    closest.new.t060...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf "*****ERROR: Can't both ignore downstream and force downstream.\n*****\n\nexit 1\n" > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D ref -id -fd > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t061: an unknown -t mode is an error
###########################################################
echo -e "    closest.new.t061...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: Request "all", "first", "last" for Tie Mode (-t)\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -t bogus > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t062: an unknown -D mode is an error
###########################################################
echo -e "    closest.new.t062...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: -D option must be followed with "ref", "a", or "b"\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -D up > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t063: -k 0 is an error
###########################################################
echo -e "    closest.new.t063...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '***** ERROR: -k option must be followed by a positive integer value *****\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -k 0 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t064: a non-numeric -k is an error
###########################################################
echo -e "    closest.new.t064...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '***** ERROR: -k option must be followed by a positive integer value *****\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -k two > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t065: -sortout with -d is an error
###########################################################
echo -e "    closest.new.t065...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '*****ERROR: -sortout (sorted output) is not valid with distance reporting.\n*****\n\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -sortout > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t066: -s with -S is an error
###########################################################
echo -e "    closest.new.t066...\c"
printf 'chr1\t100\t200\ta1\t0\t+\n' > a.tmp.bed
printf 'chr1\t300\t310\tb1\t0\t+\n' > b.tmp.bed
printf '***** ERROR: request -s for sameStrand, or -S for diffStrand, not both. *****\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -s -S > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t067: -f outside (0, 1] is an error
###########################################################
echo -e "    closest.new.t067...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf '***** ERROR: -f must be in the range (0.0, 1.0]. *****\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -f 1.5 > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t068: -names with fewer names than databases is an
#  error
###########################################################
echo -e "    closest.new.t068...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b1.tmp.bed
printf 'chr1\t400\t410\n' > b2.tmp.bed
printf '***** ERROR: Number of database name tags given does not match number of databases. *****\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -names only > out.tmp 2>&1 && echo "exit 0" >> out.tmp || echo "exit $?" >> out.tmp
sed -n '/ERROR/,$p' out.tmp > obs
check obs exp
rm obs exp out.tmp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  closest.new.t069: -header prints the query's header lines, both
#  "#" and "track" lines, before the results
###########################################################
echo -e "    closest.new.t069...\c"
printf '# comment\ntrack name=a\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t300\n' > b.tmp.bed
printf '# comment\ntrack name=a\nchr1\t100\t200\tchr1\t250\t300\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -header > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t070: header lines are skipped without -header
#  Same input as t069.
###########################################################
echo -e "    closest.new.t070...\c"
printf '# comment\ntrack name=a\nchr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t300\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t250\t300\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t071: an unsorted query is an error
#  chr1:50-60 starts before chr1:100-200 on the same chromosome.
###########################################################
echo -e "    closest.new.t071...\c"
printf 'chr1\t100\t200\nchr1\t50\t60\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file a.tmp.bed has the following out of order record\nchr1\t50\t60\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t072: an unsorted database is an error
###########################################################
echo -e "    closest.new.t072...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\nchr1\t20\t30\n' > b.tmp.bed
printf 'Error: Sorted input specified, but the file b.tmp.bed has the following out of order record\nchr1\t20\t30\nexit 1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed > obs 2>&1 && echo "exit 0" >> obs || echo "exit $?" >> obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t073: the query can be read from stdin
#  Same records as t001.
###########################################################
echo -e "    closest.new.t073...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t300\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t250\t300\t51\n' > exp
cat a.tmp.bed | $BT closest -a stdin -b b.tmp.bed -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t074: a gzipped database is read like plain text
#  Same records as t006.
###########################################################
echo -e "    closest.new.t074...\c"
printf 'chr1\t100\t200\nchr1\t1000\t1100\n' > a.tmp.bed
printf 'chr1\t60\t90\nchr1\t220\t250\nchr1\t950\t990\nchr1\t1105\t1200\n' | gzip -c > b.tmp.bed.gz
printf 'chr1\t100\t200\tchr1\t60\t90\t11\nchr1\t1000\t1100\tchr1\t1105\t1200\t6\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed.gz -d > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed.gz

###########################################################
#  closest.new.t075: -N skips an overlapping B with the same name
#  A chr1:100-200 is named x. B x (chr1:150-160) overlaps but has
#  the same name; of the others, y (chr1:210-220) is 11 away and
#  z (chr1:20-30) is 71 away.
###########################################################
echo -e "    closest.new.t075...\c"
printf 'chr1\t100\t200\tx\n' > a.tmp.bed
printf 'chr1\t20\t30\tz\nchr1\t150\t160\tx\nchr1\t210\t220\ty\n' > b.tmp.bed
printf 'chr1\t100\t200\tx\tchr1\t210\t220\ty\t11\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -N > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

###########################################################
#  closest.new.t076: -sortout orders the hits of several databases
#  by position
#  db 1 has chr1:300-310, db 2 has chr1:150-160. By default the
#  lines follow database order (1, 2); with -sortout the db 2 hit,
#  which starts first, comes first.
###########################################################
echo -e "    closest.new.t076...\c"
printf 'chr1\t100\t200\n' > a.tmp.bed
printf 'chr1\t300\t310\n' > b1.tmp.bed
printf 'chr1\t150\t160\n' > b2.tmp.bed
printf 'chr1\t100\t200\t2\tchr1\t150\t160\nchr1\t100\t200\t1\tchr1\t300\t310\n' > exp
$BT closest -a a.tmp.bed -b b1.tmp.bed b2.tmp.bed -sortout > obs
check obs exp
rm obs exp a.tmp.bed b1.tmp.bed b2.tmp.bed

###########################################################
#  closest.new.t077: -nobuf writes the same lines, including a null
#  record
#  chr1 A finds chr1:250-300 at 250-200+1 = 51; chr2 A has no B.
###########################################################
echo -e "    closest.new.t077...\c"
printf 'chr1\t100\t200\nchr2\t100\t200\n' > a.tmp.bed
printf 'chr1\t250\t300\n' > b.tmp.bed
printf 'chr1\t100\t200\tchr1\t250\t300\t51\nchr2\t100\t200\t.\t-1\t-1\t-1\n' > exp
$BT closest -a a.tmp.bed -b b.tmp.bed -d -nobuf > obs
check obs exp
rm obs exp a.tmp.bed b.tmp.bed

STARTWD=$(pwd);
for ADDITIONAL_TEST in \
    sortAndNaming/test-sort-and-naming.sh \
    kclosest/test-kclosest.sh \
; do
    # In case the cd operation fails, combine it with the script execution
    cd $(dirname "${STARTWD}/${ADDITIONAL_TEST}") \
        && bash $(basename "${STARTWD}/${ADDITIONAL_TEST}") \
        || FAILURES=$(expr $FAILURES + 1);
done

[[ $FAILURES -eq 0 ]] || exit 1;
