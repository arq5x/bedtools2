#!/usr/bin/env perl

use strict;
use warnings;
use File::Basename;
use File::Spec;
use Getopt::Long qw(GetOptionsFromArray);
use IPC::Run qw( run timeout );
use POSIX qw(ceil);

my $AUTHOR 		= 'Luke Goodsell <luke.goodsell@ogt.com>';
my $VERSION 	= '0.1';
my $SCRIPTNAME 	= basename($0);
my $USAGE 		= <<EOUSAGE;
USAGE:
    $SCRIPTNAME [OPTIONS] [-i INPUT_FILE] [-o OUTPUT_FILE]

DESCRIPTION:
    This script performs an inner-subtraction on a BED file. That is, it
    subtracts from each region specified by a BED file the regions above it in
    the same BED file.
    
    Regions should be defined in the file from highest priority to lowest 
    priority. Each region has from it only the higher priority regions 
    subtracted.
    
    Eg:

        chr1    100     200     regionA
        chr1    10      400     regionB
        chr1    40      60      regionC

    ... becomes ...

        chr1    100     200     regionA
        chr1    10      100     regionB
        chr1    200     400     regionB
    
    The subtraction is performed by Bedtools' subtractBed. This script is used
    to automate and optimise the subtraction process. Subtraction is therefore
    conducted via three intermediate files: two input and one output. These
    three files are saved by default in the current working directory, though
    this can be changed with the --tempdir option. They are repeatedly rewritten
    for each subtraction round and, by default, deleted after successful
    completion.
    
    This script is designed for O(log2(N)) complexity; the number of
    subtractions it must perform is log2(number of regions). This is hugely more
    efficient than iteratively subtracting from each region all the preceding
    regions. For example, a file with 100,000 regions (lines) would require just
    17 subtractions rather than 99,999, as would otherwise be needed.
    
    This optimisation is achieved by iteratively dividing the regions into 2
    groups: the 'subtractor' group and the 'subtractee'. The regions in the
    subtractor group are then subtracted from the regions in the subtractee
    group. The division of the groups is arranged such that a lower priority
    region is only placed in the subtractor group after any higher priority
    regions have been subtracted from it; therefore, lower priority regions will
    have no effect on higher priority regions.
    
    Each region, R, is numbered Ri=1 to imax, where imax is the number of 
    regions in the input BED file. Each subtraction round, S, is likewise 
    numbered Sn=1 to nmax, where nmax = ceil(log2(imax)). The regions are split 
    into alternating subtractor and subtractee groups, Snj=1 to jmax, where 
    jmax = 2^n. Where j is odd, regions are in the subtractor group. Where j is 
    even, regions are in the subtractee group.

ARGUMENTS:
    -i --in FILE
        Default: STDIN
        The input BED file.
    -o --out FILE
        Default: STDOUT
        Where to write the output file.
    -m --messages FILE
        Default: STDERR
        File to write messages (information not part of the output data). This
        is typically status information, metadata and error notifications.
    -f --force --force-overwrite
        Default: no
        Forces the script to overwrite the output file if it exists.
    -t --tempdir DIRECTORY
        Default: .
        Where to write the temporary files used during the inner subtraction.
    -d --[no]delete-temp-files
        Default: yes
        When finished, this option determines whether the temporary files should
        be deleted after execution finishes succesfully.
    --show-matrix-only
        Default: no
        If specified, the script will generate and print the subtraction matrix
        and then exit without performing the subtraction.
    -s --silent
        Default: no
        If specified, the script will display no status/progress information.
    -v --version
        Print the version information.
    -h -? --usage --help
        Print this usage information.

EXAMPLES:
    $SCRIPTNAME -v
    
    sort -rn -k 5,5 somefile.bed | $SCRIPTNAME | sortBed

EOUSAGE

$|++;
select((select(STDERR), $|++)[0]);

use constant {
	BEDFILE_COLUMNS	=> [qw(
		chrom
		chromStart
		chromEnd
		name
		score
		strand
		thickStart
		thickEnd
		itemRgb
		blockCount
		blockSizes
		blockStarts
	)],
};

# Subroutine prototypes/declarations:
sub main(@);
sub get_ab_subtractionmatrix($);
sub print_ab_subtractionmatrix($$$$);
sub parse_args(@);
sub version(;$);
sub usage(;$);
sub numlines_infile($);
sub parse_bed_line($);
sub serialise_bed_line($);

# The program:
exit main(@ARGV);

# Subroutine definitions:

# required    int main(@)
# The primary subroutine that controls program flow
# Returns the exit status (0 for success, 1 for failure)
sub main(@) {
	my $parameters = parse_args(@_);
	
	my $outfile_fh;
	if($parameters->{output_file} ne '-') {
		open($outfile_fh, '>', $parameters->{output_file})
			or die "Couldn't open output file, '" . $parameters->{output_file} . "', for writing: $!";
	} else {
		$outfile_fh = \*STDOUT;
	}
	
	my $message_fh;
	if($parameters->{messages_file} ne '-') {
		open($message_fh, '>', $parameters->{messages_file})
			or die "Couldn't open message file, '" . $parameters->{messages_file} . "', for writing: $!";
	} else {
		if($parameters->{silent}) {
			open($message_fh, '>', File::Spec->devnull())
				or die "Couldn't open null device, '" . File::Spec->devnull() . "', for writing: $!";
		} else {
			$message_fh = \*STDOUT;
		}
	}
	
	# Where to save the temporary input/output files
	$parameters->{file_a} = $parameters->{temp_directory} . '/bed_innersubtract_A.bed';
	$parameters->{file_b} = $parameters->{temp_directory} . '/bed_innersubtract_B.bed';
	$parameters->{file_out} = $parameters->{temp_directory} . '/bed_innersubtract_out.bed';
	
	# Used for interpreting input files
	my $regionname_to_regionnum = {};
	
	# Internal store of region data. After subtraction, there may be multiple
	# sub regions, so this is a hashref of listrefs of region data.
	my $regionnum_to_regioninfos = {};
	
	# Since the input file may have non-unique region names, each input region
	# is assigned its own pseudoregion name to track regions between subtraction
	# operations. This hashref is used for finding the original name after
	# subtractions have all been complete.
	my $pseudoregionname_to_originalregionname = {};
	
	my $num_entry_cols; # Used for maintaining consistent input and output.
	
	print $message_fh "Reading input bed file...\n";
	my $num_input_lines;
	if($parameters->{input_file} ne '-') {
		$num_input_lines = numlines_infile($parameters->{input_file});
	}
	
	{
		my $in_fh;
		if($parameters->{input_file} ne '-') {
			open($in_fh, '<', $parameters->{input_file})
				or die "Couldn't open input file, '" . $parameters->{input_file} . "', for reading: $!";
		} else {
			$in_fh = \*STDIN;
		}
		
		my $last_percent;
		my $current_regionnum = 0;
		while(my $bed_line = parse_bed_line(<$in_fh>)) {
			if(defined($num_input_lines)) {
				my $this_percent = sprintf("%5.1f", ($. / $num_input_lines) * 100); #/
				if(!defined($last_percent) or $last_percent < $this_percent) {
					print $message_fh "$this_percent% complete\r";
				}
			}
			
			$bed_line->{type} eq 'entry' or next;
			defined($num_entry_cols) or $num_entry_cols = $bed_line->{num_entry_cols};
			$bed_line->{num_entry_cols} == $num_entry_cols or die "num entry columns inconsistent on line $.";
			
			$current_regionnum++;
			my $pseudoregionname = sprintf("Region_%012d", $current_regionnum);
			$pseudoregionname_to_originalregionname->{ $pseudoregionname } = $bed_line->{entry_data}->{name};
			$bed_line->{entry_data}->{name} = $pseudoregionname;
			
			$regionnum_to_regioninfos->{ $current_regionnum } = [ $bed_line ];
			$regionname_to_regionnum->{ $bed_line->{entry_data}->{name} } = $current_regionnum;
		}
		close($in_fh);
		if(defined($num_input_lines)) {
			print $message_fh ' ' x 14 . "\r";
		}
		$current_regionnum >= 1 or die "No entries read from input bed file";
	}
	
	print $message_fh "Calculating subtraction matrix...\n";
	
	my $regionnum_to_subtractionnum_to_ismaster = get_ab_subtractionmatrix(scalar(keys(%{ $regionnum_to_regioninfos })));
	
	if($parameters->{show_matrix_only}) {
		print_ab_subtractionmatrix(
			$outfile_fh,
			$regionnum_to_subtractionnum_to_ismaster,
			$regionnum_to_regioninfos,
			$pseudoregionname_to_originalregionname
		);
		return 0;
	}
	
	my $num_subtractions = scalar(keys(%{ $regionnum_to_subtractionnum_to_ismaster->{1} }));
	print $message_fh "Performing $num_subtractions subtractions...\n";
	
	foreach my $subtractionnum (sort({$a <=> $b} keys(%{ $regionnum_to_subtractionnum_to_ismaster->{1} }))) {
		print $message_fh "    subtraction $subtractionnum...\r";
		
		open(my $outfileA_fh, '>', $parameters->{file_a})
			or die "Couldn't open temp output file A, '" . $parameters->{file_a} . "', for writing: $!";
		
		open(my $outfileB_fh, '>', $parameters->{file_b})
			or die "Couldn't open temp output file B, '" . $parameters->{file_b} . "', for writing: $!";
		
		foreach my $regionnum (sort({$a <=> $b} keys(%{ $regionnum_to_regioninfos }))) {
			foreach my $entry_data (@{ $regionnum_to_regioninfos->{ $regionnum } }) {
				if($regionnum_to_subtractionnum_to_ismaster->{ $regionnum }->{ $subtractionnum }) {
					print $outfileB_fh serialise_bed_line($entry_data);
				} else {
					print $outfileA_fh serialise_bed_line($entry_data);
				}
			}
		}
		
		if(!run(
			[
				@{ $parameters->{subtractBed_command} },
				'-a', $parameters->{file_a},
				'-b', $parameters->{file_b}
			],
			'>', $parameters->{file_out}, '2>', \my $stderr, , timeout( 30 )
		)) {
			my $errno = $!;
			my $status = $?;
			my $message;
			if ($status == -1) {
				$message = "failed to execute: $errno";
			} elsif ($status & 127) {
				$message = sprintf(
					"child died with signal %d, %s coredump",
					($status & 127),
					($status & 128) ? 'with' : 'without'
				);
			} else {
				$message = sprintf("child exited with value %d", $status >> 8);
			}
			die "Subtract command failed. Error message: " . $message . "\nSTDERR:\n$stderr";
		} 
		
		open(my $in_fh, '<', $parameters->{file_out})
			or die "Couldn't open output file, '" . $parameters->{file_out} . "', for reading: $!";
		
		my $new_regionnum_to_regioninfos = {};
		while(my $bed_line = parse_bed_line(<$in_fh>)) {
			$bed_line->{type} eq 'entry' or die "Non-entry line in output";
			$bed_line->{num_entry_cols} = $num_entry_cols;
			my $this_regionnnum = $regionname_to_regionnum->{ $bed_line->{entry_data}->{name} };
			defined($new_regionnum_to_regioninfos->{ $this_regionnnum }) or $new_regionnum_to_regioninfos->{ $this_regionnnum } = [];
			push(@{ $new_regionnum_to_regioninfos->{ $this_regionnnum } }, $bed_line);
		}
		
		close($in_fh);
		
		foreach my $regionnum (sort({$a <=> $b} keys(%{ $regionnum_to_regioninfos }))) {
			if(!exists($new_regionnum_to_regioninfos->{ $regionnum })) {
				if(!$regionnum_to_subtractionnum_to_ismaster->{ $regionnum }->{ $subtractionnum }) {
					delete($regionnum_to_regioninfos->{ $regionnum });
				}
			} else {
				$regionnum_to_regioninfos->{ $regionnum } = $new_regionnum_to_regioninfos->{ $regionnum };
			}
		}
	}
	
	print $message_fh ' ' x 25 . "\r";
	
	print $message_fh "Writing Output file...\n";
	
	foreach my $regionnum (sort({$a <=> $b} keys(%{ $regionnum_to_regioninfos }))) {
		foreach my $entry (@{ $regionnum_to_regioninfos->{ $regionnum } }) {
			$entry->{entry_data}->{name} = $pseudoregionname_to_originalregionname->{ $entry->{entry_data}->{name} };
			
			print $outfile_fh serialise_bed_line($entry);
		}
	}
	
	if($parameters->{delete_temporary_files}) {
		foreach my $temp_file_param (qw(file_a file_b file_out)) {
			!-f $parameters->{ $temp_file_param } or unlink($parameters->{ $temp_file_param });
		}
	}
	
	print $message_fh "Done.\n";
	
	return 0;
}

# hashref get_ab_subtractionmatrix($)
# Arguments:
#   required    integer   $imax
# Produce a matrix of region num (i) to subtraction num (n) to an integer, where
# 1 means the region is in the subtractor list for that i,n and 0 means the
# region is in the subtractee list.
sub get_ab_subtractionmatrix($) {
	my($imax) = @_;
	
	my $i_to_n_to_isSubtractor = {};
	
	my $nmax = ceil(log($imax) / log(2));
	
	for(my $n = 1; $n <= $nmax; $n++) {
		my $jmax = 2 ** $n;
		for(my $i=1; $i <= $imax; $i ++) {
			my $j = ceil(($i * $jmax) / $imax);
			$i_to_n_to_isSubtractor->{ $i }->{ $n } = $j % 2;
		}
	}
	
	return $i_to_n_to_isSubtractor;
}

# void print_ab_subtractionmatrix($$$$)
# Arguments:
#   required    scalar    $filehandle
#   required    hashref   $regionnum_to_subtractionnum_to_isSubtractor
#   required    hashref   $regionnum_to_regioninfos
#   required    hashref   $pseudoregionname_to_originalregionname
# Print the provided $regionnum_to_subtractionnum_to_isSubtractor subtraction
# matrix to the specified filehandle, using the provided regionnum/name
# arguments to display the name of the region in the input file.
sub print_ab_subtractionmatrix($$$$) {
	my ($filehandle, $regionnum_to_subtractionnum_to_isSubtractor, $regionnum_to_regioninfos, $pseudoregionname_to_originalregionname) = @_;
	
	my $imax = scalar(keys(%{ $regionnum_to_subtractionnum_to_isSubtractor }));
	my $nmax = scalar(keys(%{ $regionnum_to_subtractionnum_to_isSubtractor->{1} }));
	
	print $filehandle join("\t",
			'name',
			'i \ n',
			(1 .. $nmax)
		) . "\n\n";
	
	for(my $i=1; $i <= $imax; $i ++) {
		my $original_region_name = $pseudoregionname_to_originalregionname->{ $regionnum_to_regioninfos->{ $i }->[0]->{entry_data}->{name} };
		print $filehandle join("\t",
			$original_region_name,
			$i,
			map({ $regionnum_to_subtractionnum_to_isSubtractor->{ $i }->{ $_ } } (1 .. $nmax))
		) . "\n";
	}
	
	return;
}

# hashref parse_args(@)
# Arguments:
#   required    list      @args
# Parse the provided command line parameters into a hashref
sub parse_args(@) {
	my @args = @_;
	
	my $parameters = {};
	$parameters->{subtractBed_command} = [];
	
	{
	    # Capture warnings from GetOptions
	    local $SIG{__WARN__} = sub { my $error = $_[0]; chomp($error); usage($error); };
		GetOptionsFromArray(\@args,
			'i|in=s'						=> \$parameters->{input_file},
			'o|out=s'						=> \$parameters->{output_file},
			'm|messages=s'					=> \$parameters->{messages_file},
			
			't|tempdir=s'					=> \$parameters->{temp_directory},
			'd|delete-temp-files!'			=> \$parameters->{delete_temporary_files},
			
			'show-matrix-only!'				=> \$parameters->{show_matrix_only},
			
			'subtractBed-command=s{,}'		=> $parameters->{subtractBed_command},
			
			'f|force|force-overwrite!'		=> \$parameters->{force_overwrite},
			's|silent!'						=> \$parameters->{silent},
			
			'h|?|usage|help'				=> \&usage,
			'v|version'						=> \&version,
		) or usage();
	}
	
	defined($parameters->{input_file}) or $parameters->{input_file} = '-';
	
	$parameters->{input_file} eq '-' or -f $parameters->{input_file}
		or usage("No such input file: " . $parameters->{input_file});
	
	(defined($parameters->{output_file}) and length($parameters->{output_file}) > 0)
		or $parameters->{output_file} = '-';
	$parameters->{output_file} eq '-' or ! -f $parameters->{output_file} or $parameters->{force_overwrite}
		or usage("Cannot write to output file, '$parameters->{output_file}', as it already exists. This can be overriden with the --force-overwrite parameter.");
	
	(defined($parameters->{messages_file}) and length($parameters->{messages_file}) > 0)
		or $parameters->{messages_file} = '-';
	$parameters->{messages_file} eq '-' or ! -f $parameters->{messages_file} or $parameters->{force_overwrite}
		or usage("Cannot write to messages file, '$parameters->{messages_file}', as it already exists. This can be overriden with the --force-overwrite parameter.");
	
	$parameters->{input_file} ne $parameters->{output_file} or $parameters->{input_file} eq '-'
		or usage("In-place manipulation is not supported, as it is too error-prone. Please specify a different output file.");
	
	defined($parameters->{temp_directory}) or $parameters->{temp_directory} = '.';
	-d $parameters->{temp_directory}
		or usage("Specified temporary directory, '" . $parameters->{temp_directory} . "' is not a directory");
	
	defined($parameters->{delete_temporary_files}) or $parameters->{delete_temporary_files} = 1;
	
	scalar(@{ $parameters->{subtractBed_command} }) > 0
		or push(@{ $parameters->{subtractBed_command} }, 'bedtools', 'subtract');
	
	return $parameters;
}

# void version(;$)
# Arguments:
#   optional    scalar   $exit_after
# Prints a formatted version string. If EXIT_AFTER is supplied, exits cleanly
# afterwards.
sub version(;$) {
	print STDERR <<EOVERSION;
$SCRIPTNAME version $VERSION
$AUTHOR

EOVERSION
	
	!(@_) or exit;
	return;
}

# void usage(;$)
# Arguments:
#   optional    scalar   $error_message
# Prints the version and usage information and exits. If ERROR_MESSAGE is 
# supplied, it prints the error message and exits with an error state, otherwise
# it exits cleanly. 
sub usage(;$) {
	my $error;
	if(@_) {
		$error = shift;
		if(defined($error) and ($error eq 'usage' or $error eq 'h' or $error eq 'help' or $error eq '?')) {
			$error = undef;
		}
	}
	
	version();
	
	if($error) {
		print STDERR "ERROR: $error\n\n";
		my $usage_num_lines = () = $USAGE =~ /\n/g;
		if($usage_num_lines > 30) {
			print STDERR "More information can be gained by running\n    $SCRIPTNAME --usage\n";
			exit 1;
			return; # In case exit() has been redefined;
		}
	}
	
	print STDERR $USAGE;
	
	if($error) {
		exit 1;
	} else {
		exit 0;
	}
}

# int numlines_infile(;$)
# Arguments:
#   required    scalar   $filepath
# Returns the number of lines in a file.
# Empty files have a line count of one.
# Files with no line-feeds, but with a non-empty first line have a line count of
# one.
# If the final line does not have a trailing line-feed, it is still counted as a
# separate line.
sub numlines_infile($) {
	my ($filepath) = @_;
	
	defined($filepath)
		or die "No filepath supplied to numlines_infile";
	
	open(my $fh, $filepath)
		or die "Couldn't open file '$filepath' for reading: $!";
	
	# Number of lines in file.
	my $num_lines = 0;
	
	# Last non-empty buffer read. Only empty if nothing read from file at all.
	my $last_read_buffer = '';
	
	while(sysread($fh, my $buffer, 4096)) {
		$last_read_buffer = $buffer;
		$num_lines += ($buffer =~ tr/\n/\n/);
	}
	
	# Some boundary condition checks:
	if(length($last_read_buffer) > 0) {
		if($num_lines == 0) {
			# non-empty file with no line-feeds
			$num_lines ++;
		} elsif(substr($last_read_buffer, -1, 1) ne "\n") {
			# non-empty file whose last line wasn't terminated with a line-feed
			$num_lines ++;
		}
	}
	
	close $fh;
	return $num_lines;
}

# hashref parse_bed_line($)
# Arguments:
#   required    scalar   $line
# Parses the supplied BED file line into a hashref of data
sub parse_bed_line($) {
	my ($line) = @_;
	
	defined($line) or return;
	$line =~ s/\015?\012?$//; # Support windows/unix/mac line endings
	length($line) > 0 or return;
	
	my $line_data = {
		type			=> undef,
		entry_data		=> undef,
		browser_data	=> undef,
		track_data		=> undef,
		num_entry_cols	=> undef,
	};
		
	if($line =~ /^chr/i) {
		# entry line
		$line_data->{type} = 'entry';
		$line_data->{entry_data} = {};
		
		my @values = split("\t", $line, -1);
		scalar(@values) >= 3 or die "Bad bed file line: $line";
		$line_data->{num_entry_cols} = scalar(@values);
		
		for(my $col_idx = 0; $col_idx <= $#values; $col_idx++) {
			if($col_idx < 12) {
				my $column_name = BEDFILE_COLUMNS->[ $col_idx ];
				$line_data->{entry_data}->{ $column_name } = $values[$col_idx];
			} else {
				defined($line_data->{entry_data}->{rest})
					or $line_data->{entry_data}->{rest} = [];
				push(@{ $line_data->{entry_data}->{rest} }, $values[$col_idx]);
			}
		}
	} elsif($line =~ /^browser /) {
		# browser line
		$line_data->{type} = 'browser';
		$line_data->{browser_data} = {};
		
		my @elements = quotewords('\s+', 0, $line);
		scalar(@elements) > 2
			or die "Bad browser line: $line";
		my ($browser, $key, @values) = @elements;
		$line_data->{browser_data}->{ $key } = \@values;
		
	} elsif($line =~ /^track/) {
		# track line
		$line_data->{type} = 'track';
		$line_data->{track_data} = {};
		
		my @mappings = quotewords('\s+', 0, $line);
		shift(@mappings);
		foreach my $mapping (@mappings) {
			$mapping =~ /^([^=]+)=(.+)$/ or die "Bad track variable mapping: $mapping";
			my ($key, $value) = ($1, $2);
			!defined($line_data->{track_data}->{ $key })
				or die "Multiple assignments of '$key' property in track line: $line";
			
			$line_data->{track_data}->{ $key } = $value;
		}
		
	} else {
		die "Bad bedfile line: $line";
	}
	
	return $line_data;
}

# hashref serialise_bed_line($)
# Arguments:
#   required    hashref  $line
# Serialises the supplied BED line data into a string
sub serialise_bed_line($) {
	my ($line) = @_;
	
	scalar(keys(%{ $line })) > 0
		or die "No data provided to write";
	
	defined($line->{type})
		or die "No line type provided";
	
	my $line_string;
	
	if($line->{type} eq 'browser') {
		my @browser_vars = keys(%{ $line->{browser_data} });
		scalar(@browser_vars) == 1
			or die "Invalid number of browser vars for a single line: " . scalar(@browser_vars);
		my $browser_var = shift(@browser_vars);
		
		$line_string = "browser $browser_var " . join(" ", @{ $line->{browser_data}->{ $browser_var } });
		
	} elsif($line->{type} eq 'track') {
		$line_string = "track";
		foreach my $track_var (keys(%{ $line->{track_data} })) {
			my $this_var_value = $line->{track_data}->{ $track_var };
		
			$this_var_value =~ s/\\/\\\\/g;
			$this_var_value =~ s/"/\\"/g;
			$this_var_value = '"' . $this_var_value . '"';
			
			$line_string .= " ${track_var}=${this_var_value}";
		}
		
	} elsif($line->{type} eq 'entry') {
		my $num_columns = $line->{num_entry_cols};
		$num_columns ||= 12;
		
		my @values_to_print = ();
		
		for(my $col_idx = 0; $col_idx < $num_columns; $col_idx++) {
			my $column_name = BEDFILE_COLUMNS->[ $col_idx ];
			push(@values_to_print, (defined($line->{entry_data}->{ $column_name }) ? $line->{entry_data}->{ $column_name } : ""));
		}
		
		if(defined($line->{entry_data}->{rest})) {
			push(@values_to_print, @{ $line->{entry_data}->{rest} });
		}
		
		$line_string = join("\t", @values_to_print);
		
	} else {
		die "Unexpected line type: '$line->{type}'";
	}
	
	return $line_string  . "\n";
}


__END__







