ORG=hg19
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -D $ORG -P 3306   -e "select chrom,chromStart,chromEnd,K.transcript,X.geneSymbol from knownCanonical as K,kgXref as X where  X.kgId=K.transcript;" | sort -k1,1 -k2,2n > $ORG.knownCanonical.bed
