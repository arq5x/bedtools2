fail=0
pass=0
failed=""
passed=""

for dir in $(ls); do
    [ -d $dir ] || continue
    echo "Testing bedtools $dir:"
    cd $dir
    bash test-$dir.sh
    if [ $? -eq 0 ]; then
        pass=$((pass + 1))
        passed="$passed $dir"
    else
        fail=$((fail + 1))
        failed="$failed $dir"
    fi
    cd -
done

echo "--"
echo "Results"
echo "pass: $pass"
echo "passed: $passed"
echo "fail: $fail"
echo "failed: $failed"
