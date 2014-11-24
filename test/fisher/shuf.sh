obs=$(bedtools intersect -wo -a taa.bed -b tbb.bed | wc -l)

seq 1000 | xargs -P 7 -n 1 bash -c "bedtools intersect -wo -a taa.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i tbb.bed -g tgg.genome) | wc -l " | awk -vobs=$obs '{s += ($1 > obs)}END{print s / NR}'

