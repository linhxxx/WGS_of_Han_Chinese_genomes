#!/usr/bin/perl
use strict;
use warnings;


my $fa_file=shift;
die "perl $0 fa_file\n" unless $fa_file;

open IN,$fa_file or die $!;
$/=">";$/=<IN>;$/="\n";
while (<IN>){
	chomp;
	(my $id=$_)=~s/\s+.*$//;
	
	$/=">";
	my $seq=<IN>;
	chomp $seq;
	$seq=~s/\s+//g;
	$/="\n";

	my $length = length($seq);
	$seq =~ s/N//g;
	my $seq_len = length($seq);
	print "$id\t$length\t$seq_len\n";
}
close IN;

