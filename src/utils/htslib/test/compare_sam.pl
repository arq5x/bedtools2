#!/usr/bin/perl -w
#
#    Copyright (C) 2013 Genome Research Ltd.
#
#    Author: James Bonfield <jkb@sanger.ac.uk>
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

# Compares two SAM files to report differences.
# Optionally can skip header or ignore specific types of diff.

use strict;
use Getopt::Long;

my %opts;
GetOptions(\%opts, 'noqual', 'noaux', 'notemplate', 'unknownrg', 'nomd', 'template-1', 'noflag', 'Baux');

my ($fn1, $fn2) = @ARGV;
open(my $fd1, "<", $fn1) || die $!;
open(my $fd2, "<", $fn2) || die $!;

# Headers
my ($c1,$c2)=(1,1);
my (@hd1, @hd2, $ln1, $ln2);
while (<$fd1>) {
    if (/^@/) {
        push(@hd1, $_);
    } else {
        $ln1 = $_;
        last;
    }
    $c1++;
}

while (<$fd2>) {
    if (/^@/) {
        push(@hd2, $_);
    } else {
        $ln2 = $_;
        last;
    }
    $c2++;
}

# FIXME: to do
#print "@hd1\n";
#print "@hd2\n";

# Compare lines
while ($ln1 && $ln2) {
    $ln1 =~ s/\015?\012/\n/;
    $ln2 =~ s/\015?\012/\n/;
    chomp($ln1);
    chomp($ln2);

    # Java CRAM adds RG:Z:UNKNOWN when the read-group is absent
    if (exists $opts{unknownrg}) {
        $ln1 =~ s/\tRG:Z:UNKNOWN//;
        $ln2 =~ s/\tRG:Z:UNKNOWN//;
    }

    if (exists $opts{nomd}) {
        $ln1 =~ s/\tMD:Z:[A-Z0-9^]*//;
        $ln2 =~ s/\tMD:Z:[A-Z0-9^]*//;
        $ln1 =~ s/\tNM:i:\d+//;
        $ln2 =~ s/\tNM:i:\d+//;
    }

    my @ln1 = split("\t", $ln1);
    my @ln2 = split("\t", $ln2);

    # Fix BWA bug: unmapped data should have no alignments
    if ($ln1[1] & 4) { $ln1[4] = 0; $ln1[5] = "*"; }
    if ($ln2[1] & 4) { $ln2[4] = 0; $ln2[5] = "*"; }

    # Canonicalise floating point numbers
    map {s/^(..):f:(.*)/{"$1:f:".($2+0)}/e} @ln1[11..$#ln1];
    map {s/^(..):f:(.*)/{"$1:f:".($2+0)}/e} @ln2[11..$#ln2];


    if (exists $opts{Baux}) {
        # Turn ??:H:<hex> into ??:B:c,<vals> so we can compare
        # Cramtools.jar vs htslib encodings.  Probably doable with (un)pack
        map {s/^(..):H:(.*)/{join(",", "$1:B:C", map {hex $_} $2=~m:..:g)}/e} @ln1[11..$#ln1];
        map {s/^(..):H:(.*)/{join(",", "$1:B:C", map {hex $_} $2=~m:..:g)}/e} @ln2[11..$#ln2];

        # Canonicalise ??:B:? data series to be unsigned
        map {s/^(..):B:c(,?)(.*)/{"$1:B:C$2".join(",",map {($_+256)&255} split(",",$3))}/e} @ln1[11..$#ln1];
        map {s/^(..):B:c(,?)(.*)/{"$1:B:C$2".join(",",map {($_+256)&255} split(",",$3))}/e} @ln2[11..$#ln2];

        map {s/^(..):B:s(,?)(.*)/{"$1:B:S$2".join(",",map {($_+65536)&65535} split(",",$3))}/e} @ln1[11..$#ln1];
        map {s/^(..):B:s(,?)(.*)/{"$1:B:S$2".join(",",map {($_+65536)&65535} split(",",$3))}/e} @ln2[11..$#ln2];

        map {s/^(..):B:i(,?)(.*)/{"$1:B:I$2".join(",",map {$_<0? ($_+4294967296) : $_} split(",",$3))}/e} @ln1[11..$#ln1];
        map {s/^(..):B:i(,?)(.*)/{"$1:B:I$2".join(",",map {$_<0? ($_+4294967296) : $_} split(",",$3))}/e} @ln2[11..$#ln2];
    }

    # Rationalise order of auxiliary fields
    if (exists $opts{noaux}) {
        @ln1 = @ln1[0..10];
        @ln2 = @ln2[0..10];
    } else {
        #my @a=@ln1[11..$#ln1];print "<<<@a>>>\n";
        @ln1[11..$#ln1] = sort @ln1[11..$#ln1];
        @ln2[11..$#ln2] = sort @ln2[11..$#ln2];
    }

    if (exists $opts{noqual}) {
        $ln1[10] = "*";
        $ln2[10] = "*";
    }

    if (exists $opts{notemplate}) {
        @ln1[6..8] = qw/* 0 0/;
        @ln2[6..8] = qw/* 0 0/;
    }

    if (exists $opts{noflag}) {
        $ln1[1] = 0; $ln2[1] = 0;
    }

    if (exists $opts{'template-1'}) {
        if (abs($ln1[8] - $ln2[8]) == 1) {
            $ln1[8] = $ln2[8];
        }
    }

    # Cram doesn't uppercase the reference
    $ln1[9] = uc($ln1[9]);
    $ln2[9] = uc($ln2[9]);

    # Cram will populate a sequence string that starts as "*"
    $ln2[9] = "*" if ($ln1[9] eq "*");

    # Fix 0<op> cigar fields
    $ln1[5] =~ s/(\D|^)0\D/$1/g;
    $ln1[5] =~ s/^$/*/g;
    $ln2[5] =~ s/(\D|^)0\D/$1/g;
    $ln2[5] =~ s/^$/*/g;

    # Fix 10M10M cigar to 20M
    $ln1[5] =~ s/(\d+)(\D)(\d+)(\2)/$1+$3.$2/e;
    $ln2[5] =~ s/(\d+)(\D)(\d+)(\2)/$1+$3.$2/e;

    if ("@ln1" ne "@ln2") {
        print "Diff at lines $fn1:$c1, $fn2:$c2\n";
        my @s1 = split("","@ln1");
        my @s2 = split("","@ln2");
        my $ptr = "";
        for (my $i=0; $i < $#s1; $i++) {
            if ($s1[$i] eq $s2[$i]) {
                $ptr .= "-";
            } else {
                last;
            }
        }
        print "1\t@ln1\n2\t@ln2\n\t$ptr^\n\n";
        exit(1);
    }

    $ln1 = <$fd1>;
    $ln2 = <$fd2>;

    $c1++; $c2++;
}

if (defined($ln1)) {
    print "EOF on $fn1\n";
    exit(1);
}

if (defined($ln2)) {
    print "EOF on $fn2\n";
    exit(1);
}

close($fd1);
close($fd2);

exit(0);
