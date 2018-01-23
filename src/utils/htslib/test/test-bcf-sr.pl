#!/usr/bin/env perl
#
# Author: petr.danecek@sanger
#
# Test bcf synced reader's allele pairing
#

use strict;
use warnings;
use Carp;
use Data::Dumper;
use List::Util 'shuffle';
use File::Temp qw/ tempfile tempdir /;
use FindBin;
use lib "$FindBin::Bin";

my $opts = parse_params();
run_test($opts);

exit;

#--------------------------------

sub error
{
    my (@msg) = @_;
    if ( scalar @msg ) { confess @msg; }
    print 
        "Usage: test-bcf-sr.pl [OPTIONS]\n",
        "Options:\n",
        "   -s, --seed <int>        Random seed\n",
        "   -t, --temp-dir <dir>    When given, temporary files will not be removed\n",
        "   -v, --verbose           \n",
        "   -h, -?, --help          This help message\n",
        "\n";
    exit -1;
}
sub parse_params
{
    my $opts = {};
    while (defined(my $arg=shift(@ARGV)))
    {
        if ( $arg eq '-t' || $arg eq '--temp-dir' ) { $$opts{keep_files}=shift(@ARGV); next }
        if ( $arg eq '-v' || $arg eq '--verbose' ) { $$opts{verbose}=1; next }
        if ( $arg eq '-s' || $arg eq '--seed' ) { $$opts{seed}=shift(@ARGV); next }
        if ( $arg eq '-?' || $arg eq '-h' || $arg eq '--help' ) { error(); }
        error("Unknown parameter \"$arg\". Run -h for help.\n");
    }
    $$opts{tmp} = exists($$opts{keep_files}) ? $$opts{keep_files} : tempdir(CLEANUP=>1);
    if ( $$opts{keep_files} ) { cmd("mkdir -p $$opts{keep_files}"); }
    if ( !exists($$opts{seed}) )
    {
        $$opts{seed} = time();
        print STDERR "Random seed is $$opts{seed}\n";
    }
    srand($$opts{seed});
    return $opts;
}

sub _cmd
{
    my ($cmd) = @_;
    my $kid_io;
    my @out;
    my $pid = open($kid_io, "-|");
    if ( !defined $pid ) { error("Cannot fork: $!"); }
    if ($pid)
    {
        # parent
        @out = <$kid_io>;
        close($kid_io);
    }
    else
    {
        # child
        exec('bash', '-o','pipefail','-c', $cmd) or error("Cannot execute the command [/bin/sh -o pipefail -c $cmd]: $!");
    }
    return ($? >> 8, join('',@out));
}
sub cmd
{
    my ($cmd) = @_;
    my ($ret,$out) = _cmd($cmd);
    if ( $ret ) { error("The command failed [$ret]: $cmd\n", $out); }
    return $out;
}

sub save_vcf
{
    my ($opts,$vars,$fname) = @_;
    open(my $fh,"| $FindBin::Bin/../bgzip -c > $fname") or error("$FindBin::Bin/../bgzip -c > $fname: !");
    print $fh qq[##fileformat=VCFv4.3\n];
    print $fh qq[##FILTER=<ID=PASS,Description="All filters passed">\n];
    print $fh qq[##contig=<ID=1>\n];
    print $fh qq[##contig=<ID=2>\n];
    print $fh '#'. join("\t", qw(CHROM POS ID  REF ALT QUAL    FILTER  INFO))."\n";
    for my $var (@$vars)
    {
        my @als = split(/,/,$var);
        my @alts = ();
        my $ref;
        for my $al (@als)
        {
            my ($xref,$alt) = split(/>/,$al);
            $ref = $xref;
            push @alts,$alt;
        }
        print $fh join("\t", (1,100,'.',$ref,join(',',@alts),'.','.','.'))."\n";
    }
    for my $var (@$vars)
    {
        my @als = split(/,/,$var);
        my @alts = ();
        my $ref;
        for my $al (@als)
        {
            my ($xref,$alt) = split(/>/,$al);
            $ref = $xref;
            push @alts,$alt;
        }
        print $fh join("\t", (1,300,'.',$ref,join(',',@alts),'.','.','.'))."\n";
    }
    for my $var (@$vars)
    {
        my @als = split(/,/,$var);
        my @alts = ();
        my $ref;
        for my $al (@als)
        {
            my ($xref,$alt) = split(/>/,$al);
            $ref = $xref;
            push @alts,$alt;
        }
        print $fh join("\t", (2,100,'.',$ref,join(',',@alts),'.','.','.'))."\n";
    }
    close($fh) or error("close failed: bgzip -c > $fname");
    cmd("$FindBin::Bin/../tabix -f $fname");
}

sub random_alt
{
    my ($ref,$is_snp) = @_;
    my @acgt = qw(A C G T);
    my $alt = $acgt[rand @acgt];
    if ( $ref eq $alt ) { return '.'; }   # ref
    if ( !$is_snp ) { $alt = $ref.$alt; }
    return $alt;
}

sub check_outputs
{
    my ($fname_bin,$fname_perl) = @_;
    my %out = ();
    open(my $fh,'<',$fname_bin) or error("$fname_bin: $!");
    while (my $line=<$fh>)
    {
        my ($pos,@vals) = split(/\t/,$line);
        chomp($vals[-1]);
        $vals[-1] =~ s/\r$//;
        push @{$out{$pos}},join("\t",@vals);
    }
    close($fh) or error("close failed: $fname_bin");
    if ( keys %out != 3 ) { error("Expected 3 positions, found ",scalar keys %out,": $fname_bin\n"); }
    my $n;
    for my $pos (keys %out)
    {
        if ( !defined $n ) { $n = scalar @{$out{$pos}}; }
        if ( @{$out{$pos}} != $n ) { error("Expected $n positions, found ",scalar keys %{$out{$pos}},"\n"); }
    }
    my @blines = @{$out{(keys %out)[0]}};

    my @plines = ();
    open($fh,'<',$fname_perl) or error("$fname_perl: $!");
    while (my $line=<$fh>)
    {
        chomp($line);
        $line =~ s/\r$//;
        push @plines,$line;
    }
    close($fh) or error("close failed: $fname_perl");
    if ( @blines != @plines ) { error("Different number of lines: ",scalar @blines," vs ",scalar @plines," in $fname_bin vs $fname_perl\n"); }
    @blines = sort @blines;
    @plines = sort @plines;
    for (my $i=0; $i<@plines; $i++)
    {
        if ( $blines[$i] ne $plines[$i] )
        {
            #error("Different lines in $fname_bin vs $fname_perl:\n\t$blines[$i].\nvs\n\t$plines[$i].\n"); 
            error("Different lines in $fname_bin vs $fname_perl:\n\t".join("\n\t",@blines)."\nvs\n\t".join("\n\t",@plines)."\n"); 
        }
    }
}

sub run_test
{
    my ($opts) = @_;
    my @acgt = qw(A C G T);
    my $ref  = $acgt[rand @acgt];
    my @vcfs = ();
    my $nvcf = 1 + int(rand(10));
    for (my $i=0; $i<$nvcf; $i++)
    {
        my %vars  = ();
        my $nvars = 1 + int(rand(6));
        for (my $j=0; $j<$nvars; $j++)
        {
            my $snp = int(rand(2));
            my $alt = random_alt($ref,$snp);
            my $var = "$ref>$alt";
            if ( $alt ne '.' && !int(rand(5)) )    # create multiallelic site
            {
                my $alt2 = random_alt($ref,$snp);
                if ( $alt2 ne '.' && $alt ne $alt2 )
                {
                    $var .= ",$ref>$alt2"; 
                }
            }
            $vars{$var} = 1;
        }
        my $ndup = 1 + int(rand(4));
        for (my $j=0; $j<$ndup; $j++)
        {
            my @keys = shuffle keys %vars;
            push @vcfs, \@keys;
        }
    }
    @vcfs = shuffle @vcfs;
    open(my $fh,'>',"$$opts{tmp}/list.txt") or error("$$opts{tmp}/list.txt: $!");
    my %groups = ();
    my @group_list = ();
    for (my $i=0; $i<@vcfs; $i++)
    {
        my $vcf = $vcfs[$i];
        my $key = join(';',sort @$vcf);
        if ( !exists($groups{$key}) )
        {
            push @group_list,$key;
            $groups{$key}{vars} = [@$vcf];
            $groups{$key}{key} = $key;
        }
        push @{$groups{$key}{vcfs}},$i;
        save_vcf($opts,$vcf,"$$opts{tmp}/$i.vcf.gz");
        print $fh "$$opts{tmp}/$i.vcf.gz\n";
    }
    close($fh);

    my @groups = ();
    for my $group (@group_list) { push @groups, $groups{$group}; }
    for my $logic (qw(snps indels both snps+ref indels+ref both+ref exact some all))
    #for my $logic (qw(snps))
    {
        print STDERR "$FindBin::Bin/test-bcf-sr $$opts{tmp}/list.txt -p $logic > $$opts{tmp}/rmme.bin.out\n" unless !$$opts{verbose};
        cmd("$FindBin::Bin/test-bcf-sr $$opts{tmp}/list.txt -p $logic > $$opts{tmp}/rmme.bin.out");

        open(my $fh,'>',"$$opts{tmp}/rmme.perl.out") or error("$$opts{tmp}/rmme.perl.out: $!");
        $$opts{fh} = $fh;
        $$opts{logic} = $logic;
        pair_lines($opts,\@groups);
        close($fh) or error("close failed: $$opts{tmp}/rmme.perl.out");

        check_outputs("$$opts{tmp}/rmme.bin.out","$$opts{tmp}/rmme.perl.out");
    }
}

sub pair_lines
{
    my ($opts,$groups) = @_;

    #print 'groups: '.Dumper($groups);

    # get a list of all unique variants and their groups
    my %vars = ();
    my @var_list = ();
    for (my $igrp=0; $igrp<@$groups; $igrp++)
    {
        my $grp = $$groups[$igrp];
        for (my $ivar=0; $ivar<@{$$grp{vars}}; $ivar++)
        {
            my $var = $$grp{vars}[$ivar];
            if ( !exists($vars{$var}) ) { push @var_list,$var; }  # just to keep the order
            push @{$vars{$var}}, { igrp=>$igrp, ivar=>$ivar, cnt=>scalar @{$$grp{vcfs}} };
        }
    }

    # each variant has a list of groups that it is present in
    my @vars = ();
    for my $var (@var_list) { push @vars, $vars{$var}; }

    #print STDERR 'unique variants: '.Dumper(\@var_list);
    # for (my $i=0; $i<@vars; $i++)
    # {
    #     my $igrp = $vars[$i][0]{igrp};
    #     my $jvar = $vars[$i][0]{ivar};
    #     my $var  = $$groups[$igrp]{vars}[$jvar];
    #     print STDERR "$i: $var\n"; 
    # }

    # initialize variant sets - combinations of compatible variants across multiple reader groups
    my @var_sets = ();
    for (my $i=0; $i<@vars; $i++) { push @var_sets,[$i]; }

    my @bitmask = ();
    my @pmatrix = ();
    for (my $iset=0; $iset<@var_sets; $iset++)
    {
        $pmatrix[$iset] = [(0) x (scalar @$groups)];
        $bitmask[$iset] = 0;
    }
    my @max;
    for (my $iset=0; $iset<@var_sets; $iset++)
    {
        my $tmp_max = 0;
        for my $ivar (@{$var_sets[$iset]})
        {
            my $var = $vars[$ivar];
            for my $grp (@$var)
            {
                my $igrp = $$grp{igrp};
                $pmatrix[$iset][$igrp] += $$grp{cnt};
                if ( $bitmask[$iset] & (1<<$igrp) ) { error("Uh!"); }
                $bitmask[$iset] |= 1<<$igrp;
                $tmp_max += $$grp{cnt};
            }
        }
        push @max, $tmp_max;
    }

    # pair the lines
    while ( @var_sets )
    {
        my $imax = 0;
        for (my $iset=1; $iset<@var_sets; $iset++)
        {
            if ( $max[$iset] > $max[$imax] ) { $imax = $iset; }
        }
        # if ( @var_sets == @vars ) { dump_pmatrix($groups,\@vars,\@var_sets,\@pmatrix,\@bitmask); }

        my $ipair = undef;
        my $max_score = 0;
        for (my $iset=0; $iset<@var_sets; $iset++)
        {
            if ( $bitmask[$imax] & $bitmask[$iset] ) { next; }     # cannot merge
            my $score = pairing_score($opts,$groups,\@vars,$var_sets[$imax],$var_sets[$iset]);
            if ( $max_score < $score ) { $max_score = $score; $ipair = $iset; }
        }

        # merge rows thus creating a new variant set
        if ( defined $ipair && $ipair != $imax )
        {
            $imax = merge_rows($groups,\@vars,\@var_sets,\@pmatrix,\@bitmask,\@max,$imax,$ipair);
            next;
        }

        output_row($opts,$groups,\@vars,\@var_sets,\@pmatrix,\@bitmask,\@max,$imax);
        # dump_pmatrix($groups,\@vars,\@var_sets,\@pmatrix,\@bitmask);
    }
}

sub merge_rows
{
    my ($grps,$vars,$var_sets,$pmat,$bitmask,$max,$ivset,$jvset) = @_;
    if ( $ivset > $jvset ) { my $tmp = $ivset; $ivset = $jvset; $jvset = $tmp; }
    push @{$$var_sets[$ivset]}, @{$$var_sets[$jvset]};
    for (my $igrp=0; $igrp<@{$$pmat[$ivset]}; $igrp++)
    {
        $$pmat[$ivset][$igrp] += $$pmat[$jvset][$igrp];
    }
    $$max[$ivset] += $$max[$jvset];
    $$bitmask[$ivset] |= $$bitmask[$jvset];
    splice(@$var_sets,$jvset,1);
    splice(@$pmat,$jvset,1);
    splice(@$bitmask,$jvset,1);
    splice(@$max,$jvset,1);
    return $ivset;
}

sub output_row
{
    my ($opts,$grps,$vars,$var_sets,$pmat,$bitmask,$max,$ivset) = @_;
    my $varset = $$var_sets[$ivset];
    my @tmp = ();
    for my $grp (@$grps)
    {
        for my $vcf (@{$$grp{vcfs}}) { push @tmp, '-'; }
    }
    for my $idx (@$varset)
    {
        for my $var (@{$$vars[$idx]})
        {
            my $igrp = $$var{igrp};
            my $jvar = $$var{ivar};
            my $str  = $$grps[$igrp]{vars}[$jvar];
            $str =~ s/[^>]>//g;
            for my $ivcf (@{$$grps[$igrp]{vcfs}}) { $tmp[$ivcf] = $str; }
        }
    }
    print {$$opts{fh}} join("\t",@tmp)."\n";
    splice(@$var_sets,$ivset,1);
    splice(@$pmat,$ivset,1);
    splice(@$bitmask,$ivset,1);
    splice(@$max,$ivset,1);
}

sub dump_pmatrix
{
    my ($grps,$vars,$var_sets,$pmat,$bitmask) = @_;
    for (my $ivset=0; $ivset<@$var_sets; $ivset++)
    {
        my $varset = $$var_sets[$ivset];
        my @tmp = ();
        for my $ivar (@$varset)
        {
            my $igrp = $$vars[$ivar][0]{igrp};
            my $jvar = $$vars[$ivar][0]{ivar};
            push @tmp, $$grps[$igrp]{vars}[$jvar];
        }
        printf STDERR "%-10s",join(',',@tmp);
        for (my $igrp=0; $igrp<@{$$pmat[0]}; $igrp++)
        {
            print STDERR "\t$$pmat[$ivset][$igrp]";
        }
        print STDERR "\n";
    }
    print STDERR "\n";
}

sub var_type
{
    my ($vars) = @_;
    my %type = ();
    for my $var (split(/,/,$vars))
    {
        my ($ref,$alt) = split(/>/,$var);
        if ( $ref eq $alt or $alt eq '.' ) { $type{ref} = 1; }
        elsif ( length($ref)==length($alt) && length($ref)==1 ) { $type{snp} = 1; }
        else { $type{indel} = 1; }
    }
    return keys %type;
}
sub multi_is_subset
{
    my ($avar,$bvar) = @_;
    my %avars = ();
    my %bvars = ();
    for my $var (split(/,/,$avar)) { $avars{$var} = 1; }
    for my $var (split(/,/,$bvar)) { $bvars{$var} = 1; }
    for my $var (keys %avars)
    {
        if ( exists($bvars{$var}) ) { return 1; }
    }
    for my $var (keys %bvars)
    {
        if ( exists($avars{$var}) ) { return 1; }
    }
    return 0;
}
sub multi_is_exact
{
    my ($avar,$bvar) = @_;
    my %avars = ();
    my %bvars = ();
    for my $var (split(/,/,$avar)) { $avars{$var} = 1; }
    for my $var (split(/,/,$bvar)) { $bvars{$var} = 1; }
    for my $var (keys %avars)
    {
        if ( !exists($bvars{$var}) ) { return 0; }
    }
    for my $var (keys %bvars)
    {
        if ( !exists($avars{$var}) ) { return 0; }
    }
    return 1;
}
sub pairing_score
{
    my ($opts,$grps,$vars,$avset,$bvset) = @_;

    my $score = {};
    if ( $$opts{logic}=~/both/ or $$opts{logic}=~/snps/ or $$opts{logic}=~/all/ ) 
    { 
        $$score{snp}{snp} = 3;
        if ( $$opts{logic}=~/ref/ or $$opts{logic}=~/all/ ) { $$score{snp}{ref} = 2; }
    }
    if ( $$opts{logic}=~/both/ or $$opts{logic}=~/indels/ or $$opts{logic}=~/all/ ) 
    {
        $$score{indel}{indel} = 3; 
        if ( $$opts{logic}=~/ref/ or $$opts{logic}=~/all/ ) { $$score{indel}{ref} = 2; }
    }
    if ( $$opts{logic}=~/all/ )
    {
        $$score{snp}{indel} = 1;
        $$score{indel}{snp} = 1;
    }
    for my $a (keys %$score)
    {
        for my $b (keys %{$$score{$a}})
        {
            $$score{$b}{$a} = $$score{$a}{$b};
        }
    }

    my $max_int = 0xFFFFFFFF;
    my $min = $max_int;
    for my $ia (@$avset)
    {
        for my $ib (@$bvset)
        {
            my $avar = $$grps[ $$vars[$ia][0]{igrp} ]{vars}[ $$vars[$ia][0]{ivar} ];
            my $bvar = $$grps[ $$vars[$ib][0]{igrp} ]{vars}[ $$vars[$ib][0]{ivar} ];

            if ( $avar eq $bvar ) { return $max_int; }
            if ( $$opts{logic} eq 'exact' )
            {
                if ( multi_is_exact($avar,$bvar) ) { return $max_int; }
                next;
            }
            elsif ( multi_is_subset($avar,$bvar) ) { return $max_int; }

            my @atype = var_type($avar);
            my @btype = var_type($bvar);
            my $max = 0;
            for my $a (@atype)
            {
                for my $b (@btype)
                {
                    if ( !exists($$score{$a}{$b}) ) { next; }
                    if ( $max < $$score{$a}{$b} ) { $max = $$score{$a}{$b}; }
                }
            }
            if ( !$max ) { return 0; }      # some of the variants in the two groups are not compatible
            if ( $min > $max ) { $min = $max; }
        }
    }
    if ( $$opts{logic} eq 'exact' ) { return 0; }

    my $cnt = 0;
    for my $ivar (@$avset,@$bvset)
    {
        my $var = $$vars[$ivar];
        for my $grp (@$var)
        {
            $cnt += $$grp{cnt};
        }
    }
    return (1<<(28+$min)) + $cnt;
}


