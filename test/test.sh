set -e;

STARTWD=$(pwd);
FAILURES=0;
TOOL_PASSES=""
TOOL_FAILURES=""

for tool in $(ls); do
    [ -d $tool ] || continue
    echo "Testing bedtools $tool:"
    (cd $tool bash && bash test-$tool.sh) || FAILURES=$(expr $FAILURES + 1);
    if [ $FAILURES -eq 0 ]; then
        TOOL_PASSES="$TOOL_PASSES $tool"
    else
        TOOL_FAILURES="$TOOL_FAILURES $tool"
    fi
done

echo
echo
echo "--------------------------"
echo " Test Results             "
echo "--------------------------"
echo "Tools passing: $TOOL_PASSES"
echo "Tools failing: $TOOL_FAILURES"

if [[ $FAILURES -gt 0 ]]; then
    exit 1;
fi
