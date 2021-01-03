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
# test basic shuffle
###########################################################
echo -e "    shuffle.t1...\c"
echo \
"chr3	192943497	192943965	trf	789
chr4	13668420	13668593	trf	346
chr1	114076345	114076585	trf	434
chr13	99270316	99270538	trf	273
chr7	64692734	64692911	trf	187
chr3	160795229	160795394	trf	199
chr11	112718526	112718664	trf	242
chr10	84388736	84388771	trf	70
chr15	74373271	74373368	trf	79
chr7	100517230	100517271	trf	73" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head > obs
check obs exp
rm obs exp


###########################################################
# test basic shuffle with -incl (choose intervals randomly)
###########################################################
echo -e "    shuffle.t2...\c"
echo \
"chr3	494824	495292	trf	789
chr3	155662	155835	trf	346
chr5	978428	978668	trf	434
chr2	566144	566366	trf	273
chr1	2524257	2524434	trf	187
chr1	974662	974827	trf	199
chr3	511406	511544	trf	242
chr4	372392	372427	trf	70
chr3	252210	252307	trf	79
chr1	429351	429392	trf	73
chr3	637074	637106	trf	64
chr1	3632329	3632434	trf	149
chr1	1405460	1405498	trf	58
chr1	4587372	4587842	trf	278
chr3	813140	813610	trf	339
chr3	831383	831811	trf	202
chr3	177788	177831	trf	59
chr1	140167	140207	trf	62
chr3	642846	642881	trf	52
chr1	2627907	2628084	trf	302" > exp
$BT shuffle -incl incl.bed -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head -20 > obs
check obs exp
rm obs exp

##############################################################
# test basic shuffle with -incl (choose chroms randomly first)
##############################################################
echo -e "    shuffle.t3...\c"
echo \
"chr3	494824	495292	trf	789
chr3	155662	155835	trf	346
chr5	978428	978668	trf	434
chr2	566144	566366	trf	273
chr1	2524257	2524434	trf	187
chr1	974662	974827	trf	199
chr3	511406	511544	trf	242
chr4	372392	372427	trf	70
chr3	252210	252307	trf	79
chr1	429351	429392	trf	73
chr3	637074	637106	trf	64
chr1	3632329	3632434	trf	149
chr1	1405460	1405498	trf	58
chr1	4587372	4587842	trf	278
chr3	813140	813610	trf	339
chr3	831383	831811	trf	202
chr3	177788	177831	trf	59
chr1	140167	140207	trf	62
chr3	642846	642881	trf	52
chr1	2627907	2628084	trf	302" > exp
$BT shuffle -incl incl.bed -chromFirst -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head -20 > obs
check obs exp
rm obs exp


##############################################################
# test basic shuffle with -excl
##############################################################
echo -e "    shuffle.t4...\c"
echo -n "" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome \
            -excl excl.bed \
| $BT intersect -a - -b excl.bed > obs
check obs exp
rm obs exp

##############################################################
# test basic shuffle with 
##############################################################
echo -e "    shuffle.t5...\c"
echo \
"chr4	13668420	13668593	trf	346
chr1	114076345	114076585	trf	434
chr1	114076345	114076585	trf	434
chr5	17088394	17088864	trf	339
chr3	53794735	53794769	trf	68
chr2	73265723	73265766	trf	86
chr2	4749579	4749649	trf	68
chr1	15263027	15263097	trf	104
chr5	57165089	57165114	trf	50
chr4	33917224	33917392	trf	150" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome \
| $BT intersect -a - -b excl.bed | head > obs
check obs exp
rm obs exp

###############################################################
# test an interval that is bigger than the max chrom length
###############################################################
echo -e "    shuffle.t6...\c"
echo "Error, line 1: tried 1000 potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on." > exp
$BT shuffle -i <(echo -e "chr1\t0\t110") -g <(echo -e "chr1\t100") &> obs
check obs exp
rm obs exp
[[ $FAILURES -eq 0 ]] || exit 1;
