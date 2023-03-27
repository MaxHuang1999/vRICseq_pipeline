#!/usr/bin/perl
die "perl $0 class_of_AlignPair.log\n" if(@ARGV != 1);
my $class_log=shift;

my %hash;
open(CL,$class_log) || die;
while(my $line=<CL>){
	chomp $line;
	my @sub=split/\s+/,$line;
	$hash{$sub[1]."\t".$sub[2]}++;
}

foreach (sort {$hash{$b} <=> $hash{$a}} keys %hash){
	print $_,"\t",$hash{$_},"\n";
}
