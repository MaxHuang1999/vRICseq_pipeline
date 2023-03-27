#!/usr/bin/perl
die "perl $0 pcp_rep1_interaction.sam want_chr want_start want_end> xxx.bedpe\n" if(@ARGV != 4);
my $contact_sam=shift;
my $want_chr=shift;
my $want_start=shift;
my $want_end=shift;

my %unique;
my %gap_len;
open(CS,$contact_sam) || die;
while(my $frag_a=<CS>){
        if($frag_a=~/^@/){
                next;
        }
        my $frag_b=<CS>;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	my @id_info_a=split/_/,$sub_a[0];
	my @id_info_b=split/_/,$sub_b[0];

     	if($id_info_a[0]."_".$id_info_a[1] ne $id_info_b[0]."_".$id_info_b[1]){
        	next;
      	}
     	else{#same read name
        	my $chr_a=$sub_a[2];
		my ($loci_a,$end_a)=get_map_start_end($frag_a);
		$sub_a[1] = $sub_a[1] > 255 ? $sub_a[1]-256 : $sub_a[1];

  		my $chr_b=$sub_b[2];
		my ($loci_b,$end_b)=get_map_start_end($frag_b);
		$sub_b[1] = $sub_b[1] > 255 ? $sub_b[1]-256 : $sub_b[1];

		if($id_info_a[2] ne $id_info_b[2] or $id_info_a[2] ne "Plus"){	#same reads; same RNA; same strand
			next;
		}
		
		if($chr_a ne $want_chr or $chr_b ne $want_chr){
			next;
		}

                if($sub_a[0] =~ /^AlignPair/){
			next;
                }
                elsif($sub_a[0] =~ /^Part/){
                        if($end_a >= $loci_b){	#impossible situation
                                print $frag_a,$frag_b;
                                die "this Part-pair is incorrect\n";
                        }
			my $tmp_gap_len=$loci_b-$end_a-1;
			$gap_len{$tmp_gap_len}++;
			print $want_chr,"\t",$end_a,"\t",$end_a+1,"\t";
			print $want_chr,"\t",$loci_b-1,"\t",$loci_b,"\n";
		}
                elsif($sub_a[0] =~ /^Chimeric/){
			#next;
			if($sub_a[1] eq "0"){
				if($end_b > $loci_a){	#maybe some Bugs in chimeric read sets; only a little, so filtered.
					#print $frag_a,$frag_b;
					next;
				}
				print $want_chr,"\t",$loci_b-1,"\t",$loci_b,"\t";
				print $want_chr,"\t",$end_a-1,"\t",$end_a,"\n";
			}
			elsif($sub_a[1] eq "16"){
				if($end_a > $loci_b){	#maybe some Bugs in chimeric read sets; only a little, so filtered.
					#print $frag_a,$frag_b;
					next;
				}
				print $want_chr,"\t",$loci_a-1,"\t",$loci_a,"\t";
				print $want_chr,"\t",$end_b-1,"\t",$end_b,"\n";
			}
			else{
				print $frag_a,$frag_b;
				die "wrong strand\n";
			}
		}
		else{
			print $frag_a,$frag_b;
			die "wrong reads class\n";
		}
	}
}

sub get_map_start_end{
	my $frag=shift;
	my @sub=split/\s+/,$frag;
	my $cigar=$sub[5];
	my $chr_start=$sub[3];
	my $chr_end=$sub[3];
	while($cigar=~/(\d+)(\w)/g){
	    	my $tmp_len=$1;
	    	my $tmp_content=$2;
	    	if($tmp_content eq "M"){
			$chr_end+=$tmp_len;
		}
		elsif($tmp_content eq "I"){
		}
		elsif($tmp_content eq "D"){
			$chr_end+=$tmp_len;
		}
		elsif($tmp_content eq "S" or $tmp_content eq "H"){
		}
		elsif($tmp_content eq "N"){
			warn "No gap exists in chimeric segments\n";
			warn $frag;
			die;
		}
	}
	return ($chr_start,$chr_end-1);
}
