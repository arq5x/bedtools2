Fisher Testing
==============

Fisher is now based on the count of interval overlaps, subject to `-f`.
We can compare the output of fisher on simulated data by running `python sim.py`
which will show the output from `bedtools fisher` and then running `bash shuf.sh`
which will repeatedly run

```Shell

bedtools intersect -wo -a taa.bed -b <(bedtools shuffle -allowBeyondChromEnd   -i tbb.bed -g tgg.genome) | wc -l 
```

and then report the proportion of times that number is >= the observed intersection.

In `sim.py` changing the lenght of the intervals (via `maxA`, `maxB`) has the greatest effect on the
correspondence of the simulated p-value from `shuf.sh` and the one from `fisher`. The right-tailed p-value
from `fisher` should correspond well to the value from the simulation.


