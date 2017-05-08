#!/usr/bin/perl

use strict;
use warnings;
use lib "";#lib path for chisqure 
use Distribute qw(chisqure);

(my $mode, my $input_file, my $output_file, my $bp_threshold_selection, my $bp_threshold, my $s_avg_ref, my $p_avg_ref, my $s_avg_scaf, my $p_avg_scaf, my $depth_limit) = @ARGV;
my $p_limit=0.05;
die "perl $0 <mode> <file_input> <file_output> <bp_threshold_selection> <bp_threshold> <s_avg_ref> <p_avg_ref>  <s_avg_scaf> <p_avg_scaf> <depth_limit>" unless $mode && $input_file && $output_file && $bp_threshold_selection && $bp_threshold && $s_avg_ref && $p_avg_ref && $s_avg_scaf && $p_avg_scaf && $depth_limit;
die "mode should be \"insertion\" , \"deletion\" or \"inversion\"" unless $mode eq "insertion" || $mode eq "deletion" || $mode eq "inversion";

#validator judgement
my $validator_over;
my $validator_under;
if($mode eq "insertion")
{
	$validator_over = \&valid_insertion_over;
	$validator_under = \&valid_insertion_under;
}
elsif ($mode eq "deletion")
{
	$validator_over = \&valid_deletion_over;
	$validator_under = \&valid_deletion_under;
}
elsif ($mode eq "inversion")
{
	$validator_over = \&valid_inversion_over;
	$validator_under = \&valid_inversion_under;
}

#bp_threshold selector judgement
if($bp_threshold_selection eq "ref")
{
	$bp_threshold_selection = 2;
}
elsif($bp_threshold_selection eq "scaf")
{
	$bp_threshold_selection = 5;
}
else
{
	die "bp_threshold should be \"ref\" or \"scaf\"\n";
}

open my $fh, $input_file or die "Error opening input file $input_file!\n";
open my $O, ">$output_file" or die "Error opening output file $output_file!\n";

my $total=0; my $matched=0;

$/=">"; <$fh>; $/="\n";
while(my $l1=<$fh>)
{
	#data inport
	chomp $l1;
	my @a1 = split /\s+/, $l1;
	
	$/=">";
	my $seg = <$fh>;
	chomp $seg;
	$/="\n";
	
	my @a;
	my @line = split /\n/, $seg;
	for (my $i = 0; $i < @line; $i++){
		push @a, [split /\s+/, $line[$i]];
	}

	++$total;

	#data effectiveness validation
	if($a[1][0] eq "*")
	{
		next;
	}
	for(my $i=0; $i<=5; ++$i)
	{
		if (!defined($a[$i][7]) || $a[$i][7]==0) {$a[$i][7]=1;} 
		if (!defined($a[$i][8]) || $a[$i][8]==0) {$a[$i][8]=1;}
	}	
	
#data processing
	#single and both data smoother
	for(my $i=0; $i<=5; ++$i)
	{
		if(($a[$i][3]+$a[$i][4]+$a[$i][5]+$a[$i][6]) == 0)
		{
			$a[$i][3]=1;
		}
		$a[$i][7] /= ($a[$i][3]+$a[$i][4]+$a[$i][5]+$a[$i][6]);
		$a[$i][8] /= ($a[$i][3]+$a[$i][4]+$a[$i][5]+$a[$i][6]);
	}

	#chi_square
	my @chi;
	my @p;
	for(my $i=0; $i<=5; $i+=2)
	{
		($chi[$i], $p[$i]) = chisqure([ [$a[$i][7], $a[$i][8]], [$s_avg_ref, $p_avg_ref] ]);
	}
	for(my $i=1; $i<=5; $i+=2)
	{
		($chi[$i], $p[$i]) = chisqure([ [$a[$i][7], $a[$i][8]], [$s_avg_scaf, $p_avg_scaf] ]);
	}
	
	#depth
	my @depth;
	for(my $i=0; $i<=5; ++$i)
	{
		($depth[$i]) = $a[$i][7] + $a[$i][8];
	}

	#judgements
	if(&judge_threshold(@a1))
	{
		&$validator_over(\@a1, \@chi, \@p, \@depth);
	}
	else
	{
		&$validator_under(\@a1, \@chi, \@p, \@depth);
	}
}

printf STDERR "$matched / $total";

	sub valid_inversion_over
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$p[1] > $p_limit) && ($$p[5] > $p_limit))
		{
			&output_info("1", @$info);
		}
	}
	sub valid_insertion_over
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$p[0] <= $p_limit) && ($$p[4] <= $p_limit))
		{
			&output_info("1", @$info);
		}
	}
	sub valid_deletion_over
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$p[0] <= $p_limit) && ($$p[4] <= $p_limit))
		{
				&output_info("1", @$info);
		}
	}
	sub valid_inversion_under
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$depth[3] >= $depth_limit))
		{
			&output_info("0", @$info);
		}
	}
	sub valid_insertion_under
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$depth[3] >= $depth_limit))
		{
			&output_info("0", @$info);
		}
	}
	sub valid_deletion_under
	{
		(my $info, my $chi, my $p, my $depth) = @_;
		if(($$depth[2] <= $depth_limit))
		{
				&output_info("0", @$info);
		}
	}
	sub judge_threshold
	{
		my @info = @_;
		$info[$bp_threshold_selection] =~ /(\S+),(\S+),(\S+)/;
		return (($3 - $2 + 1) >= $bp_threshold);
	}
	sub output_info
	{
		my @info = @_;
		++$matched;
		$info[4] =~ s/^>(\S+)/$1/;
		$info[3] =~ /(\S+),(\S+),(\S+)/;
		my $rn=$1; my $rs=$2; my $re=$3;
		$info[6] =~ /(\S+),(\S+),(\S+)/;
		my $sn=$1; my $ss=$2; my $se=$3;
		print $O "$mode$info[0]\t$rn\t$rs\t$re\t$sn\t$ss\t$se\n";
	}
