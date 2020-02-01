#################################################################
# 
#  PERFORMANCE TESTS FOR VARIOUS INPUT METHODS.
#
#  We'll test uncompressed data, gzip, bgzip, and Bam files.
#  Each will be tried as a file, and as the 3 types of stdin:
#  redirects, pipes, and fifos.
#  
#  We'll run each one twice in a row to check for improvement 
#  due to a cache warm up.
#
#################################################################


BT=${BT-../../bin/bedtools}


#################################################################
#  Start by generating data, if desired
#  
#################################################################
if true; then
  echo "generating data..."
  mkdir perfData
  cd perfData
  ../$BT random -l 1000 -n 10000000 -g ../human.hg19.genome | sort -k1,1 -k2,2n > a10M.bed 
  ../$BT random -l 1000 -n 10000000 -g ../human.hg19.genome | sort -k1,1 -k2,2n > b10M.bed 
  cp a10M.bed a10M_gzipped.bed
  gzip a10M_gzipped.bed
  ../../htsutil bgzfcompress a10M.bed a10M_bgzipped.bed.gz
  ../$BT bedtobam -i a10M.bed -g ../human.hg19.genome > a10M.bam
  cd ..
fi

###################################################################
#  Begin Tests
###################################################################
echo -e "Test 1 of 16...."
echo "Test 1: Intersect a10M, ten million records, uncompressed from file with b10M" > runLog.txt
runit $BT intersect -a perfData/a10M.bed -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 1 for possible cache speed up." >> runLog.txt
runit $BT intersect -a perfData/a10M.bed -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 2 of 16...."
echo -e "\n\nTest 2: Intersect a10M, ten million records, gzipped from file with b10M" >> runLog.txt
runit $BT intersect -a perfData/a10M_gzipped.bed.gz -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 2 for possible cache speed up." >> runLog.txt
runit $BT intersect -a perfData/a10M_gzipped.bed.gz -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 3 of 16...."
echo -e "\n\nTest 3: Intersect a10M, ten million records, bgzipped from file with b10M" >> runLog.txt
runit $BT intersect -a perfData/a10M_bgzipped.bed.gz -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 3 for possible cache speed up." >> runLog.txt
runit $BT intersect -a perfData/a10M_bgzipped.bed.gz -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 4 of 16...."
echo -e "\n\nTest 4: Intersect a10M, ten million records, bam from file with b10M" >> runLog.txt
runit $BT intersect -a perfData/a10M.bam -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 4 for possible cache speed up." >> runLog.txt
runit $BT intersect -a perfData/a10M.bam -b perfData/b10M.bed -sorted 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt


echo -e "Test 5 of 16...."
echo -e "\n\nTest 5: Intersect a10M, ten million records, uncompressed from redirect with b10M" >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M.bed 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 5 for possible cache speed up." >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M.bed 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 6 of 16...."
echo -e "\n\nTest 6: Intersect a10M, ten million records, gzipped from redirect with b10M" >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M_gzipped.bed.gz 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 6 for possible cache speed up." >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M_gzipped.bed.gz 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 7 of 16...."
echo -e "\n\nTest 7: Intersect a10M, ten million records, bgzipped from redirect with b10M" >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M_bgzipped.bed.gz 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 7 for possible cache speed up." >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M_bgzipped.bed.gz 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 8 of 16...."
echo -e "\n\nTest 8: Intersect a10M, ten million records, bam from redirect with b10M" >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M.bam 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 8 for possible cache speed up." >> runLog.txt
runit $BT intersect -a - -b perfData/b10M.bed -sorted < perfData/a10M.bam 2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 9 of 16...."
echo -e "\n\nTest 9: Intersect a10M, ten million records, uncompressed from pipe with b10M" >> runLog.txt
cat perfData/a10M.bed | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 9 for possible cache speed up." >> runLog.txt
cat perfData/a10M.bed | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 10 of 16...."
echo -e "\n\nTest 10: Intersect a10M, ten million records, gzip from pipe with b10M" >> runLog.txt
cat perfData/a10M_gzipped.bed.gz | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 10 for possible cache speed up." >> runLog.txt
cat perfData/a10M_gzipped.bed.gz | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 11 of 16...."
echo -e "\n\nTest 11: Intersect a10M, ten million records, bgzip from pipe with b10M" >> runLog.txt
cat perfData/a10M_bgzipped.bed.gz | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 11 for possible cache speed up." >> runLog.txt
cat perfData/a10M_bgzipped.bed.gz | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 12 of 16...."
echo -e "\n\nTest 12: Intersect a10M, ten million records, bam from pipe with b10M" >> runLog.txt
cat perfData/a10M.bam | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 12 for possible cache speed up." >> runLog.txt
cat perfData/a10M.bam | runit $BT intersect -a - -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 13 of 16...."
echo -e "\n\nTest 13: Intersect a10M, ten million records, uncompressed from fifo with b10M" >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M.bed) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 13 for possible cache speed up." >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M.bed) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 14 of 16...."
echo -e "\n\nTest 14: Intersect a10M, ten million records, gzipped from fifo with b10M" >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M_gzipped.bed.gz) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 14 for possible cache speed up." >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M_gzipped.bed.gz) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 15 of 16...."
echo -e "\n\nTest 15: Intersect a10M, ten million records, bgzipped from fifo with b10M" >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M_bgzipped.bed.gz) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 15 for possible cache speed up." >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M_bgzipped.bed.gz) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt

echo -e "Test 16 of 16...."
echo -e "\n\nTest 16: Intersect a10M, ten million records, bam from fifo with b10M" >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M.bam) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt
echo -e "\nRe-do test 16 for possible cache speed up." >> runLog.txt
runit $BT intersect -a <(cat perfData/a10M.bam) -b perfData/b10M.bed -sorted  2>&1 >/dev/null | grep -e "user" -e maxrss >> runLog.txt


echo "Tests completed."
rm -rf perfData
