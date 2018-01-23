#!/bin/sh
#
#    Copyright (C) 2017 Genome Research Ltd.
#
#    Author: Robert Davies <rmd@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

# Executes a single test and compares against the expected output
run_test() {
    # Expected result: pass (P) / fail (F) / nonzero exit (N)
    p="$1"; shift
    # File with expected output (empty or '.' if none)
    e="$1"; shift
    # Test result
    r="P"
    # Why the test failed
    y=""
    if [ "x$test_iter" = "x" ]
    then
	test_iter=1
    else
	test_iter=`expr $test_iter + 1`
    fi
    result=`eval ${@+"$@"} 2>_err.tmp > _out.tmp`
    if [ $? != 0 ]
    then
	if [ "$p" != "N" ]
	then
	    # Expected zero exit code, got non-zero
	    r="F"
	    y="exit_code"
	else
	    # Expected non-zero exit code and got it
	    r="P"
	fi
    elif [ "$p" = "N" ]
    then
	# Expected non-zero exit code, but got zero
	r="F"
	y="exit_code"
    elif [ "x$e" != "x" -a "$e" != "." ]
    then
    sed -n 's/.*/&/p' _out.tmp > _out.tmp2
	if cmp -s _out.tmp2 "$e"
	then
	    # Output was as expected
	    r="P"
	    rm -f _out.tmp _out.tmp2 _err.tmp
	else
	    # Output differed
	    r="F"
	    y="output"
	fi
    else
	# Expected zero exit code and got it.
	r="P"
	rm -f _out.tmp _out.tmp2 _err.tmp
    fi

    if [ "$r" = "F" ]
    then
	# Test failed
	case "$p" in
	    [PN])
		echo "FAIL : $@"
		if [ "x$e" != "x" -a "$e" != "." ]
		then
		    keep_output="FAIL-$e.${test_iter}"
		else
		    keep_output="FAIL.${test_iter}"
		fi
		mv _out.tmp "${keep_output}.out"
		mv _err.tmp "${keep_output}.err"
		nufail=`expr $nufail + 1`
		if [ "$y" = "exit_code" ]
		then
		    if [ "$p" != "N" ]
		    then
			echo "Got non-zero exit code"
		    else
			echo "Got unexpected zero exit code"
		    fi
		    echo "See ${keep_output}.{out,err} for output"
		else
		    echo "Output differed from expected result"
		    echo "Compare $e ${keep_output}.out"
		fi
		;;
	    *)
		echo "XFAIL: $@"
		nefail=`expr $nefail + 1`
		;;
	esac
    else
	# Test passed
	case "$p" in
	    "P")
		echo "PASS : $@"
		nepass=`expr $nepass + 1`
		;;
	    "N")
		echo "PASS : $@ (must exit non-zero)"
		nepass=`expr $nepass + 1`
		;;
	    *)
		echo "XPASS: $@"
		nupass=`expr $nupass + 1`
		;;
	esac
    fi
}

tabix_test() {
    nupass=0; nepass=0
    nufail=0; nefail=0

    exec 9<"$1"
    while read -r line <&9
    do
	set -- $line
	case $1 in
            "#"*) # skip comments
		;;
            "")   # skip blank lines too
		;;

	    "INIT")
		shift
		eval ${@+"$@"} > /dev/null
		if [ $? != 0 ]
		then
		    echo "INIT FAIL: $@"
		    return 1
		fi
		;;

            *)
		p=$1;shift
		o=$1;shift
		run_test "$p" "$o" ${@+"$@"}
		;;
	esac
    done
    exec 9<&-

    echo ""
    echo "Expected   passes:   $nepass"
    echo "Unexpected passes:   $nupass"
    echo "Expected   failures: $nefail"
    echo "Unexpected failures: $nufail"
    if [ "$nupass" -gt 0 -o "$nufail" -gt 0 ]
    then
	return 1
    else
	return 0
    fi
}

echo "Testing tabix..."

bgzip="../../bgzip"
tabix="../../tabix"

tabix_test $@

exit $?
