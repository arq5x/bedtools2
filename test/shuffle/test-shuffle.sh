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
"chr1	49390118	49390586	trf	789
chr13	60419956	60420129	trf	346
chr9	88507701	88507941	trf	434
chr1	176220424	176220646	trf	273
chr7	69225121	69225298	trf	187
chr2	25938602	25938767	trf	199
chr7	133111177	133111315	trf	242
chr3	5894626	5894661	trf	70
chr7	144738233	144738330	trf	79
chr9	28526037	28526078	trf	73" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head > obs
check obs exp
rm obs exp


###########################################################
# test basic shuffle with -incl (choose intervals randomly)
###########################################################
echo -e "    shuffle.t2...\c"
echo \
"chr1	3126067	3126535	trf	789
chr5	846414	846587	trf	346
chr3	747335	747575	trf	434
chr2	451924	452146	trf	273
chr1	1837113	1837290	trf	187
chr1	1389014	1389179	trf	199
chr1	458954	459092	trf	242
chr4	267572	267607	trf	70
chr2	608295	608392	trf	79
chr3	544706	544747	trf	73
chr1	285876	285908	trf	64
chr5	961111	961216	trf	149
chr4	53213	53251	trf	58
chr1	1344769	1345239	trf	278
chr1	1516219	1516689	trf	339
chr1	1062880	1063308	trf	202
chr1	2030329	2030372	trf	59
chr1	4349564	4349604	trf	62
chr1	3052376	3052411	trf	52
chr1	2890010	2890187	trf	302" > exp
$BT shuffle -incl incl.bed -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head -20 > obs
check obs exp
rm obs exp

##############################################################
# test basic shuffle with -incl (choose chroms randomly first)
##############################################################
echo -e "    shuffle.t3...\c"
echo \
"chr1	3126067	3126535	trf	789
chr5	846414	846587	trf	346
chr3	747335	747575	trf	434
chr2	451924	452146	trf	273
chr1	1837113	1837290	trf	187
chr1	1389014	1389179	trf	199
chr1	458954	459092	trf	242
chr4	267572	267607	trf	70
chr2	608295	608392	trf	79
chr3	544706	544747	trf	73
chr1	285876	285908	trf	64
chr5	961111	961216	trf	149
chr4	53213	53251	trf	58
chr1	1344769	1345239	trf	278
chr1	1516219	1516689	trf	339
chr1	1062880	1063308	trf	202
chr1	2030329	2030372	trf	59
chr1	4349564	4349604	trf	62
chr1	3052376	3052411	trf	52
chr1	2890010	2890187	trf	302" > exp
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
"chr1	49390118	49390586	trf	789
chr1	49390118	49390586	trf	789
chr1	176220424	176220646	trf	273
chr1	176220424	176220646	trf	273
chr2	25938602	25938767	trf	199
chr3	5894626	5894661	trf	70
chr2	57803306	57803344	trf	58
chr4	13407024	13407058	trf	68
chr3	18198257	18198308	trf	66
chr1	71587801	71587968	trf	159" > exp
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
