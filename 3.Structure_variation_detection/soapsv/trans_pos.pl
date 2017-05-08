#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;


my (%info);
my ($info_file, $align_file) = @ARGV;
die "perl $0 scaffold_info_file alignment_file\n" unless $info_file && $align_file;

read_info($info_file, \%info);

open my $align_fh, $align_file or die $!;
while (<$align_fh>){
	chomp;
	my @tmp = split;

	my $target_seq = <$align_fh>;
	my $queue_seq = <$align_fh>;
	my $carritage = <$align_fh>;

	if ($tmp[7] eq "+"){
		print join(" ", @tmp) , "\n";
		print "$target_seq";
		print "$queue_seq";
		print "$carritage";
	}else{
		my $len = $info{$tmp[4]};
		my ($start, $end) = ($len - $tmp[6] + 1, $len - $tmp[5] + 1);
		print join(" ", @tmp[0..4]) . " $start $end $tmp[7] $tmp[8]\n";
		print "$target_seq";
		print "$queue_seq";
		print "$carritage";
	}
}
close $align_fh;



sub read_info{
	my ($file, $ref) = @_;
	open my $fh, "<", $file or die $!;
	while (<$fh>){
		chomp;
		my @tmp = split;
		$ref->{$tmp[0]} = $tmp[1];
	}
	close $fh;
}
