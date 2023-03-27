#!/usr/bin/perl
die "perl $0 read1.bed read2.bed covid_annotation.bed" if(@ARGV != 3);
my $read1_bed=shift;
my $read2_bed=shift;
my $covid19_anno=shift;

my %all_trans_introns;
open(CA,$covid19_anno) || die;
while(my $line=<CA>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($sub[4] !~ /NCIntron/){
		next;
	}	
	$all_trans_introns{$sub[3]}=$sub[1]."\t".$sub[2];
}

open(RA,$read1_bed) || die;
open(RB,$read2_bed) || die;
while(my $reada=<RA>){
	my $readb=<RB>;
	chomp $reada;
	chomp $readb;
	my @info_a=split/\s+/,$reada;
	my @info_b=split/\s+/,$readb;
	
	if($info_a[5] ne $info_b[5]){
		die;
	}
	if($info_a[5] eq "+"){
		my $min_fragment_len=$info_a[2]-$info_b[1];
		my $min_trans="VirusPlus";
		foreach my $trans (keys %all_trans_introns){
			my ($intron_start,$intron_end)=split/\s+/,$all_trans_introns{$trans};
			if($intron_start >= $info_b[2] and $intron_end <= $info_a[1]){
				my $dis=$info_a[2]-$info_b[1]-$intron_end+$intron_start;	
				if($dis < $min_fragment_len){
					$min_fragment_len=$dis;
					$min_trans=$trans;
				}
			}
		}
		print $info_a[3],"\t",$info_a[2]-$info_b[1],"\t",$min_fragment_len,"\t$min_trans\n";
	}
	elsif($info_a[5] eq "-"){	#report raw fragment distance
		print $info_a[3],"\t",$info_b[2]-$info_a[1],"\t",$info_b[2]-$info_a[1],"\tVirusMinus\n";
	}
}

