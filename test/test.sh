set -e; # Alert user to any uncaught error

ulimit -c unlimited

STARTWD=$(pwd);
TOOL_PASSES="";
TOOL_FAILURES="";

for tool in $(ls); do
    [ -d "${STARTWD}/${tool}" ] || continue;
    echo "Testing bedtools $tool:";
    cd "${STARTWD}/${tool}";
    bash "test-${tool}.sh" \
        && TOOL_PASSES="$TOOL_PASSES $tool" \
        || TOOL_FAILURES="$TOOL_FAILURES $tool";
done

echo
echo
echo "--------------------------"
echo " Test Results             "
echo "--------------------------"
echo "Tools passing: $TOOL_PASSES"
echo "Tools failing: $TOOL_FAILURES"
echo "NB: the 'negativecontrol' test is supposed to fail. If it wasn't caught, "
echo "something went wrong with this test script."
[ "$TOOL_FAILURES" = " negativecontrol" ] || exit 1;
