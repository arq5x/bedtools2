set -eo pipefail

echo "fisher,shuffled"

for i in $(seq 1000); do
    fisher=$(python ./sim.py | tail -1 | cut -f 2)
    shuffle=$(bash shuf.sh "bedtools intersect -a taa.bed -b tbb.bed" "wc -l" tgg.genome)
    echo "$fisher,$shuffle"
done
