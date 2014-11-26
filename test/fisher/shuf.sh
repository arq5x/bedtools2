set -eo pipefail

obs=$(bedtools intersect -wo -a taa.bed -b tbb.bed | wc -l)

p=$(seq 100 | xargs -P 7 -n 1 bash -c "bedtools intersect -wo -a taa.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i tbb.bed -g tgg.genome) | wc -l " | awk -vobs=$obs '{s += ($1 > obs)}END{print (1 + s) / (1 + NR)}')

if [ '1' -eq $(echo $p'< 0.1' | bc -l) ] || [ '1' -eq $(echo $p'> 0.9' | bc -l) ]; then
    p=$(seq 1000 | xargs -P 7 -n 1 bash -c "bedtools intersect -wo -a taa.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i tbb.bed -g tgg.genome) | wc -l " | awk -vobs=$obs '{s += ($1 > obs)}END{print (1 + s) / (1 + NR)}')
fi
if [ '1' -eq $(echo $p'< 0.01' | bc -l) ] || [ '1' -eq $(echo $p'> 0.99' | bc -l) ]; then
    p=$(seq 5000 | xargs -P 7 -n 1 bash -c "bedtools intersect -wo -a taa.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i tbb.bed -g tgg.genome) | wc -l " | awk -vobs=$obs '{s += ($1 > obs)}END{print (1 + s) / (1 + NR)}')
fi

echo $p
