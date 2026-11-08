# Shared helpers for the test-intersect-<topic>.sh files.
# Source from test/intersect with:  . ./helpers-intersect.sh <topic>
# Portable to bash 3.2 (macOS) and BSD tools: no GNU-only flags.

export LC_ALL=C

BT=${BT-../../bin/bedtools}
# The intersectBed alias lives next to the bedtools binary.
IBT=${IBT-$(dirname "$BT")/intersectBed}
htsutil=${htsutil-../htsutil}

TOPIC=${1-misc}
FAILURES=0
TESTS=0
# Every scratch file goes here; removed on exit.
TMPD=$(mktemp -d "${TMPDIR:-/tmp}/bt-intersect-$TOPIC.XXXXXX")
trap 'rm -rf "$TMPD"' EXIT

# begin <id>: announce a test case.
begin()
{
	TESTS=$((TESTS + 1))
	printf '    intersect.%s.%s...' "$TOPIC" "$1"
}

# check <observed file> <expected file>: byte comparison.
check()
{
	if diff "$1" "$2" > "$TMPD/diff.out" 2>&1; then
		echo ok
	else
		FAILURES=$((FAILURES + 1))
		echo fail
		cat "$TMPD/diff.out"
	fi
}

# rows: copy stdin to stdout with each single space turned into a tab.
# For records whose fields contain no spaces.
rows()
{
	tr ' ' '\t'
}

# expect_error <substring> <command> [args...]: the command must exit
# non-zero, write nothing to stdout, and print <substring> on stderr.
# Call after begin.
expect_error()
{
	local want=$1 status
	shift
	set +e
	"$@" > "$TMPD/err.stdout" 2> "$TMPD/err.stderr"
	status=$?
	set -e
	if [ $status -ne 0 ] && [ ! -s "$TMPD/err.stdout" ] \
		&& grep -F -e "$want" "$TMPD/err.stderr" > /dev/null; then
		echo ok
	else
		FAILURES=$((FAILURES + 1))
		echo fail
		echo "      exit $status; stderr:"; sed 's/^/        /' "$TMPD/err.stderr"
	fi
}

# expect_status <code> <command> [args...]: exact exit status; output to
# $TMPD/st.stdout and $TMPD/st.stderr for further checks. Prints nothing.
expect_status()
{
	local want=$1 status
	shift
	set +e
	"$@" > "$TMPD/st.stdout" 2> "$TMPD/st.stderr"
	status=$?
	set -e
	[ $status -eq "$want" ]
}

# finish: summary line and exit status for the test file.
finish()
{
	echo "    intersect.$TOPIC: $TESTS tests, $FAILURES failures"
	[ $FAILURES -eq 0 ]
}
