#!/usr/bin/env perl
#
#    Copyright (C) 2017 Genome Research Ltd.
#
#    Author: Anders Kaplan
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

use strict;

my $log_message_count = 0;
my $file_count = 0;
my $failure_count = 0;

sub check_log_message
{
  my ($message, $filename, $line_num) = @_;
  $log_message_count++;

  unless ($message =~ /^\"([A-Z]|%s)/)
  {
    print "$filename line $line_num:\n";
    print "Log message should begin with a capital letter: $message.\n";
    $failure_count++;
  }

  if ($message =~ /\\n\"$/)
  {
    print "$filename line $line_num:\n";
    print "Log message should NOT end with a newline: $message.\n";
    $failure_count++;
  }

  if ($message =~ /\.\"$/)
  {
    print "$filename line $line_num:\n";
    print "Log message should NOT end with a full stop: $message.\n";
    $failure_count++;
  }
}

sub check_file
{
  my ($filename) = @_;
  $file_count++;

  open(my $fh, '<', $filename) or die "Could not open $filename.";
  my $line_num = 1;
  my $line = <$fh>;
  while ($line)
  {
    if ($line =~ /hts_log_\w+\s*\(\s*(\"[^\"]*\")/)
    {
      unless ($line =~ /\\n\"\s*$/) # string constant continues on next line
      {
        check_log_message($1, $filename, $line_num);
      }
    }

    $line_num++;
    $line = <$fh>;
  }
}

sub check_dir
{
  my ($path) = @_;
  foreach my $filename (glob("$path/*.c"))
  {
    check_file($filename);
  }
}

check_dir("..");
check_dir("../cram");

print "$file_count files scanned\n";
print "$log_message_count log messages checked\n";
print "$failure_count errors found\n";
exit($failure_count > 0);
