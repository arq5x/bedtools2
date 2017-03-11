set -e;

STARTWD=$(pwd);
FAILURES=0;
TEST_PASSES=0;
TOOL_PASSES=""
TOOL_FAILURES=""

for dir in $(ls); do
    echo $dir
    [ -d $dir ] || continue
    echo "Testing bedtools $dir:"
    (cd $dir bash && bash test-$dir.sh) || FAILURES=$(expr $FAILURES + 1);
    if [ $? -eq 0 ]; then
        TEST_PASSES=$((TEST_PASSES + 1))
        TOOL_PASSES="$TOOL_PASSES $dir"
    else
        failed="$failed $dir"
    fi
done

echo
echo
echo "--------------------------"
echo " Test Results             "
echo "--------------------------"
echo "Tools passing: $TOOL_PASSES"
echo "Tools with failures: $TOOL_FAILURES"
echo "Tests failed: $FAILURES"

if [[ $FAILURES -gt 0 ]]; then
    exit 1;
fi
