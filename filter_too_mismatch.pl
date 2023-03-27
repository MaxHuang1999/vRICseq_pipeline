#!/usr/bin/perl
die "perl GM12878_rep1_interaction.unique.sam max_#of_mismatch\n" if(@ARGV != 2);
my $in_sam=shift;
my $cutoff=shift;
my %bad_reads;
open(IN,$in_sam) || die;
while(my $line=<IN>){
	chomp $line;
	my @sub=split/\s+/,$line;
	my @id_info=split/_/,$sub[0];
	my $mismatch;
	foreach my $s (@sub){
		if($s =~ /NM:i:(\d+)/){
			$mismatch=$1;
			last;
		}
	}
	if($mismatch > $cutoff){
		$bad_reads{$id_info[0]."_".$id_info[1]}=1;
	}
}	
close IN;

open(IN,$in_sam) || die;
while(my $line=<IN>){
	my @sub=split/\s+/,$line;
	my @id_info=split/_/,$sub[0];
	if($bad_reads{$id_info[0]."_".$id_info[1]}){
		next;
	}
	else{
		print $line;
	}
}
