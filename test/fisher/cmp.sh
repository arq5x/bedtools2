set -eo pipefail

echo "fisher,shuffled"

for i in $(seq 1000); do
    fisher=$(python ./sim.py | tail -1 | cut -f 2)
    shuffle=$(bash shuf.sh)
    echo "$fisher,$shuffle"
done
