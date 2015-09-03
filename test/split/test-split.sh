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

cat /dev/zero |tr "\0" "\n" | head -n 10000 |\
awk '{S=int(rand()*1000.0);E=S+int(rand()*100000); printf("chrX\t%d\t%d\n",S,E);}' > tmp.bed

echo "    split.01.size...\c"
rm -f _tmp.*.bed
echo "_tmp.00001.bed	9943540	200
_tmp.00002.bed	9943482	201
_tmp.00003.bed	9943541	200
_tmp.00004.bed	9943561	200
_tmp.00005.bed	9943471	200
_tmp.00006.bed	9943475	200
_tmp.00007.bed	9943468	200
_tmp.00008.bed	9943487	200
_tmp.00009.bed	9943539	200
_tmp.00010.bed	9943531	200" > exp
${BT} split -i tmp.bed -p _tmp -n 50 -a size | head > _tmp.size.tsv
check exp _tmp.size.tsv


echo "    split.02.simple...\c"
rm -f _tmp.*.bed
echo "_tmp.00001.bed	9952674	200
_tmp.00002.bed	9751661	200
_tmp.00003.bed	9649058	200
_tmp.00004.bed	9929508	200
_tmp.00005.bed	9556713	200
_tmp.00006.bed	10298876	200
_tmp.00007.bed	10043102	200
_tmp.00008.bed	9781861	200
_tmp.00009.bed	9502188	200
_tmp.00010.bed	9991229	200" > exp
${BT} split -i tmp.bed -p _tmp -n 50 -a simple | head > _tmp.simple.tsv
check exp _tmp.simple.tsv

echo "    spliit.03.simple...\c"
rm -f _tmp.*.bed
echo "_tmp.00001.bed	414150	10
_tmp.00002.bed	586843	10
_tmp.00003.bed	503604	10
_tmp.00004.bed	410044	10
_tmp.00005.bed	499400	10
_tmp.00006.bed	537341	10
_tmp.00007.bed	698581	10
_tmp.00008.bed	555258	10
_tmp.00009.bed	474511	10
_tmp.00010.bed	633012	10" > exp
${BT} split -i tmp.bed -p _tmp -n 1000 -a simple | head > _tmp.simple.tsv
check exp _tmp.simple.tsv



rm -f _tmp.*.bed jeter.bed tmp.bed  _tmp.simple.tsv  _tmp.size.tsv
