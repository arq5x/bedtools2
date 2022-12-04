print("""\
@HD\tVN:1.0\tSO:coordinate
@SQ\tSN:c1\tAS:genome.txt\tLN:100""")

#y1      0       1       16      255     5M      *       0       0       *       *
for i in range(1000000):
    print(f'r{i}\t0\tc1\t1\t100\t100M\t*\t0\t0\t*\t*')
