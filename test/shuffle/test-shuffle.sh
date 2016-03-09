BT=${BT-../../bin/bedtools}
check()
{
	if diff $1 $2; then
    	echo ok
		return 1
	else
    	echo fail
		return 0
	fi
}


###########################################################
# test basic shuffle
###########################################################
echo "    shuffle.t1...\c"
echo \
"chr9	108600879	108601347	trf	789
chr12	9186177	9186350	trf	346
chr8	89726287	89726527	trf	434
chr8	40323278	40323500	trf	273
chr8	69904335	69904512	trf	187
chr5	138240459	138240624	trf	199
chr11	96382483	96382621	trf	242
chr8	105834146	105834181	trf	70
chrX	105921488	105921585	trf	79
chrX	125331456	125331497	trf	73" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head > obs
check obs exp
rm obs exp


###########################################################
# test basic shuffle with -incl (choose intervals randomly)
###########################################################
echo "    shuffle.t2...\c"
echo \
"chr3	542223	542691	trf	789
chr5	444343	444516	trf	346
chr1	2520601	2520841	trf	434
chr5	194760	194982	trf	273
chr1	2121545	2121722	trf	187
chr1	2246343	2246508	trf	199
chr1	2724117	2724255	trf	242
chr4	304892	304927	trf	70
chr2	332618	332715	trf	79
chr5	822410	822451	trf	73
chr1	1450982	1451014	trf	64
chr1	3218361	3218466	trf	149
chr4	338952	338990	trf	58
chr3	713207	713677	trf	278
chr1	4378307	4378777	trf	339
chr1	4451988	4452416	trf	202
chr1	1545567	1545610	trf	59
chr1	573175	573215	trf	62
chr4	931201	931236	trf	52
chr1	4215777	4215954	trf	302" > exp
$BT shuffle -incl incl.bed -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head -20 > obs
check obs exp
rm obs exp

##############################################################
# test basic shuffle with -incl (choose chroms randomly first)
##############################################################
echo "    shuffle.t3...\c"
echo \
"chr5	310009	310477	trf	789
chr4	520601	520774	trf	346
chr2	130650	130890	trf	434
chr1	3246343	3246565	trf	273
chr2	968160	968337	trf	187
chr3	332618	332783	trf	199
chr4	638727	638865	trf	242
chr3	218361	218396	trf	70
chr1	2259217	2259314	trf	79
chr3	378307	378348	trf	73
chr4	447387	447419	trf	64
chr3	573175	573280	trf	149
chr2	106791	106829	trf	58
chr3	618697	619167	trf	278
chr2	211901	212371	trf	339
chr5	656883	657311	trf	202
chr2	993338	993381	trf	59
chr2	713531	713571	trf	62
chr2	428268	428303	trf	52
chr4	590632	590809	trf	302" > exp
$BT shuffle -incl incl.bed -chromFirst -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome | head -20 > obs
check obs exp
rm obs exp


##############################################################
# test basic shuffle with -excl
##############################################################
echo "    shuffle.t4...\c"
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
echo "    shuffle.t5...\c"
echo \
"chr1	150415830	150415862	trf	64
chr1	150415830	150415862	trf	64
chr5	78078743	78079213	trf	339
chr4	84711820	84712248	trf	202
chr4	61777751	61777794	trf	59
chr3	28583223	28583400	trf	302
chr1	55933709	55934092	trf	712
chr1	55933709	55934092	trf	712
chr1	39686691	39686725	trf	68
chr2	2555287	2555330	trf	86" > exp
$BT shuffle -seed 42 -i simrep.bed  \
            -g ../../genomes/human.hg19.genome \
| $BT intersect -a - -b excl.bed | head > obs
check obs exp
rm obs exp

###############################################################
# test an interval that is bigger than the max chrom length
###############################################################
echo "    shuffle.t6...\c"
echo "Error, line 1: tried 1000 potential loci for entry, but could not avoid excluded regions.  Ignoring entry and moving on." > exp
$BT shuffle -i <(echo -e "chr1\t0\t110") -g <(echo -e "chr1\t100") &> obs
check obs exp
rm obs exp