#!/usr/bin/perl
die "perl $0 test.sam" if(@ARGV != 1);
my $align_sam_1=shift;


$tmp_count=pair_part($align_sam_1);
print "Part_from_Align_Read1:\t",$tmp_count,"\n";

sub pair_part{
	my $sam=shift;
	my $count=0;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		chomp $line;
		if($line=~/^@/){
			next;
		}
		else{
			my @sub=split/\s+/,$line;
			my $cigar=$sub[5];
			my $start_in_read=0;
			my $start_in_genome=$sub[3];

			my @cigar_blocks=split/N/,$cigar;

			foreach my $i (0..$#cigar_blocks-1){
				my $end_in_read=$start_in_read;
				my $end_in_genome=$start_in_genome-1;
				my $left_part=$cigar_blocks[$i];
				my $gap_part;
				if($left_part=~s/(\d+)$//){
					$gap_part=$1;
				}
				else{
					die;
				}
				my $right_part=$cigar_blocks[$i+1];
				$right_part=~s/\d+$//;
				warn $left_part,"\t",$gap_part,"\t",$right_part,"\tgapps\n";
				#print left part
				while($left_part=~/(\d+)(\w)/g){
					my $tmp_len=$1;
					my $tmp_content=$2;
					if($tmp_content eq "M"){
						$end_in_read+=$tmp_len;
						$end_in_genome+=$tmp_len;
					}
					elsif($tmp_content eq "I"){
						$end_in_read+=$tmp_len;
					}
					elsif($tmp_content eq "D"){
						$end_in_genome+=$tmp_len;
					}
					elsif($tmp_content eq "S" or $tmp_content eq "H"){
						$end_in_read+=$tmp_len;
					}
				}

				warn $start_in_read,"\t",$end_in_read,"\taa\n";
				print $sub[0],"\t",$sub[1],"\t",$sub[2],"\t",$start_in_genome,"\t255\t",$left_part,"\t*\t0\t0\t";
				print substr($sub[9],$start_in_read,$end_in_read-$start_in_read),"\t";
				print substr($sub[10],$start_in_read,$end_in_read-$start_in_read),"\n";

				#warn $end_in_genome,"\tgg\n";

				#deal with gap
				$start_in_read=$end_in_read;
				$start_in_genome=$end_in_genome+$gap_part;

				#warn $start_in_genome,"\tgg\n";
				
				#print right part;
				while($right_part=~/(\d+)(\w)/g){
					my $tmp_len=$1;
					my $tmp_content=$2;
					if($tmp_content eq "M"){
						$end_in_read+=$tmp_len;
						$end_in_genome+=$tmp_len;
					}
					elsif($tmp_content eq "I"){
						$end_in_read+=$tmp_len;
					}
					elsif($tmp_content eq "D"){
						$end_in_genome+=$tmp_len;
					}
					elsif($tmp_content eq "S" or $tmp_content eq "H"){
						$end_in_read+=$tmp_len;
					}
				}

				warn $start_in_read,"\t",$end_in_read,"\taa\n";
			
				print $sub[0],"\t",$sub[1],"\t",$sub[2],"\t",$start_in_genome+1,"\t255\t",$right_part,"\t*\t0\t0\t";
				print substr($sub[9],$start_in_read,$end_in_read-$start_in_read),"\t";
				print substr($sub[10],$start_in_read,$end_in_read-$start_in_read),"\n";

				$start_in_genome=$start_in_genome+1;
				#$start_in_read=$end_in_read;
				#$start_in_genome=$end_in_genome+1;
				$count++;
				
			}
		}
	}
	close SM;
	return $count;
}

sub header{
	my $sam=shift;
	my $out;
	open(SM,$sam) || die;
	while(my $line=<SM>){
		if($line=~/^@/){
			$out.=$line;
		}
	}
	close SM;
	return $out;
}

