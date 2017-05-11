set -e; # Alert user to any uncaught error

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

[ "$TOOL_FAILURES" = "" ] || exit 1;
