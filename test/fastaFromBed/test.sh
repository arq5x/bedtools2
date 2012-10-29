echo $'chr1\t1\t10' | ./bin/bedtools getfasta -fi test/fastaFromBed/t.fa -bed stdin -fo t.txt

./bin/bedtools getfasta -exons -fi test/fastaFromBed/t.fa -bed test/fastaFromBed/blocks.bed -fo stdout
