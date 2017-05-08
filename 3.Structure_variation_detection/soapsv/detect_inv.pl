#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;


my (@scaffold, %strand);
my ($align_file, $inver_file) = @ARGV;
die "perl $0 alignment_file_bh2_oi inversion_output_file\n" unless   $align_file && $inver_file;



#########################

my $prev_id = "";
my $prev_chr = "";
open my $align_fh, $align_file or die $!;
open my $inver_fh, ">", $inver_file or die $!;
while (<$align_fh>){
	chomp;
	my (undef, $chr, $start_s, $end_s, $scaff_id, $start_q, $end_q, $strand, undef) = split /\s+/;

	if ($chr ne $prev_chr || $scaff_id ne $prev_id){
		if(@scaffold){
			search_inv(\@scaffold, \%strand, $prev_id);
			@scaffold = ();
			%strand = ();
		}
	}
	
	$prev_chr = $chr;
	$prev_id = $scaff_id;
	push @scaffold,  [$chr, $start_s, $end_s, $start_q, $end_q, $strand];
	$strand{$strand}++;
}
close $align_fh;

if (@scaffold){
	search_inv(\@scaffold, \%strand, $prev_id);
	@scaffold = ();
	%strand = ();
}



#########################
#####output inversion
sub search_inv{
	my ($array_p,  $strand_ref, $id) = @_;
	
	@$array_p = sort {$a->[3]<=>$b->[3]} @$array_p;
	my $strand = ($strand_ref->{'+'}) ? ( $strand_ref->{'-'} 
							? (($strand_ref->{'-'} < $strand_ref->{'+'}) ? '+' : '-') : '+') : '-';
	
	for (my $i = 0; $i < @$array_p; $i++){
    	my @scaff = @{$array_p->[$i]};
		print {$inver_fh}  "$id\tInversion\t" . join("\t", @scaff[3,4,0..2]) . "\t$strand\t$scaff[5]\n"  if $scaff[5] ne $strand; ###   inversion		
	}
}
