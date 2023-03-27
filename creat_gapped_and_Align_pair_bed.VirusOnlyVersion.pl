die "perl $0 interaction.sam\n" if(@ARGV != 1);
my $sam=shift;

my $output_part=$sam;
my $output_alignpair_A_virus=$sam;
my $output_alignpair_B_virus=$sam;

$output_part=~s/sam$/gapped.bed/;
$output_alignpair_A_virus=~s/sam$/Alignpair_read1.Virus.bed/;
$output_alignpair_B_virus=~s/sam$/Alignpair_read2.Virus.bed/;

open(PT,">$output_part") || die;
open(VAPA,">$output_alignpair_A_virus") || die;
open(VAPB,">$output_alignpair_B_virus") || die;

open(SM,$sam) || die;
while(my $frag_a=<SM>){
        if($frag_a=~/^@/){
                next;
        }
        else{
                my $frag_b=<SM>;
                my @sub_a=split/\s+/,$frag_a;
                my @sub_b=split/\s+/,$frag_b;
                my @id_a_info=split/_/,$sub_a[0];
                my @id_b_info=split/_/,$sub_b[0];
                my $strand_a=$id_a_info[2];
                my $strand_b=$id_b_info[2];
                if($id_a_info[1] ne $id_a_info[1]){     #same pair
                        print $frag_a,$frag_b;
                        die "did not belong to the same pair\n";
                }
                else{#same read name
                        my @sub_a=split/\s+/,$frag_a;
                        my @sub_b=split/\s+/,$frag_b;

                        my $chr_a=$sub_a[2];
                        my $loci_a=$sub_a[3];
			my $cigar_a=$sub_a[5];
			

                        my $chr_b=$sub_b[2];
                        my $loci_b=$sub_b[3];
			my $cigar_b=$sub_b[5];

			if($sub_a[0] =~ /^AlignPair/){
                                $cigar_a=~/(\d+)M/;
                                my $match_a=$1;
                                $cigar_b=~/(\d+)M/;
                                my $match_b=$1;

                                if($strand_a eq $strand_b){
                                        if($strand_a eq "Plus"){
                                                my $fragment_len=$loci_a+$match_a-$loci_b;
                                                if($fragment_len <= 300){	#may be normal
                                                }
                                                else{
							if($chr_a =~ /NC_019843.3/ and $chr_b =~ /NC_019843.3/){
								print VAPA $chr_a,"\t",$loci_a-1,"\t",$loci_a+$match_a-1,"\t",$sub_a[0],"\t255\t+\n";    #read 1
								print VAPB $chr_b,"\t",$loci_b-1,"\t",$loci_b+$match_b-1,"\t",$sub_b[0],"\t255\t+\n";    #read 2
							}
							else{
								die;
							}
                                                }
                                        }
                                        elsif($strand_a eq "Minus"){
                                                my $fragment_len=$loci_b+$match_b-$loci_a;
                                                if($fragment_len <= 300){
                                                }
                                                else{
							if($chr_a =~ /NC_019843.3/ and $chr_b =~ /NC_019843.3/){
								print VAPA $chr_a,"\t",$loci_a-1,"\t",$loci_a+$match_a-1,"\t",$sub_a[0],"\t255\t-\n";    #read 1
								print VAPB $chr_b,"\t",$loci_b-1,"\t",$loci_b+$match_b-1,"\t",$sub_b[0],"\t255\t-\n";    #read 2
							}
							else{
								die;
							}
                                                }
                                        }
                                }
			}
			elsif($sub_a[0] =~ /^Part/){	#overlap junction site are classified as normal;
				my $match_a=get_start_and_len($cigar_a);
				my $match_b=get_start_and_len($cigar_b);
	                        print PT $chr_a,"\t",$loci_a+$match_a-1,"\t",$match_a+$loci_a+1,"\t";
	                        print PT $chr_b,"\t",$loci_b-1,"\t",$loci_b+1,"\t";
	                        print PT $sub_a[0],"\t60\t+\t+\n";
			}
			else{	#chimeric reads; all classified as chimeric
			}
		}
	}
}

sub get_start_and_len{
        my $cigar=shift;
	my $only_block_len;
        while($cigar=~/(\d+)(\w)/g){
                my $num=$1;
                my $class=$2;
                if($class eq "S" or $class eq "H"){
                }
                elsif($class eq "M" or $class eq "D"){
                        $only_block_len+=$1;
                }
                elsif($class eq "I"){
			next;
		}
		else{
                        die;
                }
        }
	return $only_block_len;
}
	






