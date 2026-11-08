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
    # Further test files for the tool: test-<tool>-<topic>.sh
    for extra in test-${tool}-*.sh; do
        [ -f "$extra" ] || continue;
        name="${tool}/${extra%.sh}";
        name="${name#${tool}/test-}";
        echo "Testing bedtools $name:";
        bash "$extra" \
            && TOOL_PASSES="$TOOL_PASSES $name" \
            || TOOL_FAILURES="$TOOL_FAILURES $name";
    done
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
