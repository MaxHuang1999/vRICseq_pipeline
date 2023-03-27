#!/usr/bin/perl
die "perl $0 pcp_rep2_interaction.loops \n" if(@ARGV != 1);
my $loop=shift;
my $bin=1;

my $loop_num;
open(LP,$loop) || die;
while(my $line=<LP>){
	#if($loop_num >= 10000){
	#	last;
	#}
	chomp $line;
	my @sub=split/\s+/,$line;
	my $end_one_start=$sub[1];
	my $end_two_start=$sub[4];
	my $loci_one=convert($end_one_start);
	my $loci_two=convert($end_two_start);
	if($loci_one =~ /bad/ or $loci_two =~ /bad/){
		next;
	}
	else{
		$loop_num++;
		print $loci_one,"\t",$loci_two,"\t";
		my $chr_one=(split/\s+/,$loci_one)[0];
		my $chr_two=(split/\s+/,$loci_two)[0];
		my $color=choose_color($chr_one,$chr_two);
		print "color=",$color,"\n";
	}
	
}

sub choose_color{
	my $chr1=shift;
	my $chr2=shift;
	if($chr1 eq "hs1" or $chr2 eq "hs1"){
		return "purple"
	}
	elsif($chr1 eq "hs2" or $chr2 eq "hs2"){
		return "red";
	}
}

sub convert{
	my $loci=shift;
	my $new_start=$loci;
	my $new_end=$new_start+$bin-1;
	return "hs1\t".$new_start."\t".$new_end;
}

