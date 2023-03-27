#!/usr/bin/perl
die "perl $0 merge.interaction.sam starChimericRead1.sam starChimericRead2.sam bwaChimericRead1.sam bwaChimericRead2.sam\n" if(@ARGV != 5);
my $input_interaction_sam=shift;
my $star_read1_chimeric_sam=shift;
my $star_read2_chimeric_sam=shift;
my $bwa_read1_chimeric_sam=shift;
my $bwa_read2_chimeric_sam=shift;

#store head for pair;
my %bwa_read1_chimeric_head=read_bwa_chimeric($bwa_read1_chimeric_sam);
my %bwa_read2_chimeric_head=read_bwa_chimeric($bwa_read2_chimeric_sam);
my %star_read1_chimeric_head=read_star_chimeric($star_read1_chimeric_sam);
my %star_read2_chimeric_head=read_star_chimeric($star_read2_chimeric_sam);

open(IISA,$input_interaction_sam) || die;
while(my $frag_a=<IISA>){
	if($frag_a=~/^@/){
		next;
	}
	my $frag_b=<IISA>;
	my @sub_a=split/\s+/,$frag_a;
	my @sub_b=split/\s+/,$frag_b;
	my @info_a=split/_/,$sub_a[0];
	my @info_b=split/_/,$sub_b[0];
	if($info_a[1] ne $info_b[1] or $info_a[-1] ne $info_b[-1]){
		die;
	}
	
	if($sub_a[0]!~/AlignPair/){
		print $frag_a,$frag_b;
		next;
	}

	#warn $info_a[-1],"\t",$info_b[-1],"\n";
	#read 1
	if(exists $bwa_read1_chimeric_head{$info_a[-1]} and exists $star_read1_chimeric_head{$info_a[-1]}){
		warn $info_a[-1],"\taa\n";
		#die;
	}
	elsif(exists $bwa_read1_chimeric_head{$info_a[-1]} and !exists $star_read1_chimeric_head{$info_a[-1]}){
		print "AlignPairSegment_$info_a[1]","_";
		print flag_to_strand($bwa_read1_chimeric_head{$info_a[-1]},1);
		print "_";
		print $bwa_read1_chimeric_head{$info_a[-1]},"\n";
	}
	elsif(!exists $bwa_read1_chimeric_head{$info_a[-1]} and exists $star_read1_chimeric_head{$info_a[-1]}){
		print "AlignPairSegment_$info_a[1]","_";
		print flag_to_strand($star_read1_chimeric_head{$info_a[-1]},1);
		print "_";
		print $star_read1_chimeric_head{$info_a[-1]},"\n";
	}
	else{
		$frag_a=~s/AlignPairWhole/AlignPairSegment/;
		print $frag_a;
	}
	#read2
	if(exists $bwa_read2_chimeric_head{$info_b[-1]} and exists $star_read2_chimeric_head{$info_b[-1]}){
		#die;
	}
	elsif(exists $bwa_read2_chimeric_head{$info_b[-1]} and !exists $star_read2_chimeric_head{$info_b[-1]}){
		print "AlignPairSegment_$info_b[1]","_";
		print flag_to_strand($bwa_read2_chimeric_head{$info_b[-1]},2);
		print "_";
		print $bwa_read2_chimeric_head{$info_b[-1]},"\n";
	}
	elsif(!exists $bwa_read2_chimeric_head{$info_b[-1]} and exists $star_read2_chimeric_head{$info_b[-1]}){
		print "AlignPairSegment_$info_b[1]","_";
		print flag_to_strand($star_read2_chimeric_head{$info_b[-1]},2);
		print "_";
		print $star_read2_chimeric_head{$info_b[-1]},"\n";
	}
	else{
		$frag_b=~s/AlignPairWhole/AlignPairSegment/;
		print $frag_b;
	}
	#last;
}

sub flag_to_strand{
	my $frag=shift;
	my $read1_or_read2=shift;
	my $flag=(split/\s+/,$frag)[1];
	if(($flag == 0 and $read1_or_read2 == 1) or ($flag == 16 and $read1_or_read2 == 2)){
		return "Minus";
	}
	elsif(($flag == 16 and $read1_or_read2 == 1) or ($flag == 0 and $read1_or_read2 == 2)){
		return "Plus";
	}
	else{
		warn $flag,"\t",$read1_or_read2,"\tddd\n";
		die;
	}
}

sub read_star_chimeric{
	my $file=shift;
	my %hash;
	open(IN,$file) || die;
	while(my $frag_a=<IN>){
		if($frag_a =~ /^@/){
			next;
		}
		my $frag_b=<IN>;

                my @sub=split/\s+/,$frag_a;
                my @blocks;
                if($sub[5]=~/[^0-9DIMNS]/){#too complicated
                        #print $sub[5],"\n";
                        next;
                }

		if($sub[0]=~s/Head_//){
		}
		else{
			die;
		}
		$sub[1] = $sub[1] >= 256 ? $sub[1]-256 : $sub[1];
	
                my @content;
                my @lens;
                while($sub[5]=~/(\d+)(\w)/g){
                        push (@lens,$1);
                        push (@content,$2);
                }

                my $already_match=$sub[3];
                my $start_in_read;
                my $end_in_read;
                foreach my $i (0..$#content){
                        if($content[$i] eq "M"){
                                my $tmp_block=$sub[0]."\t".$sub[1]."\t".$sub[2]."\t".$already_match."\t255\t".$lens[$i]."M\t*\t0\t0\t";
                                $already_match+=$lens[$i]-1;
                                my $seq=substr($sub[9],$start_in_read,$lens[$i]);
                                my $qual=substr($sub[10],$start_in_read,$lens[$i]);
                                $start_in_read+=$lens[$i];
                                $tmp_block.=$seq."\t".$qual;
                                push (@blocks,$tmp_block);
                        }
                        elsif($content[$i] eq "N"){
                                $already_match+=$lens[$i]+1;
                        }
                        elsif($content[$i] eq "S"){
                                $start_in_read+=$lens[$i];
                        }
                        elsif($content[$i] eq "I"){
                                $start_in_read+=$lens[$i];
                        }
                        elsif($content[$i] eq "D"){
                                $already_match+=$lens[$i];
                        }
                        else{
                                die;
                        }
                }
                if($sub[1] eq "0"){
                        $hash{$sub[0]}=$blocks[0];
                }
                elsif($sub[1] eq "16"){
                        $hash{$sub[0]}=$blocks[-1];
                }
                else{
                        die;
                }
        }
        close IN;
        return %hash;
}

sub read_bwa_chimeric{
	my $file=shift;
	my %hash;
	open(IN,$file) || die;
	while(my $frag_a=<IN>){
		if($frag_a =~ /^@/){
			next;
		}
		my $frag_b=<IN>;
		if($frag_a !~ /Chimeric/){
			next;
		}
		my @sub_a=split/\s+/,$frag_a;
		my @sub_b=split/\s+/,$frag_b;
		my @info_a=split/_/,$sub_a[0];
		my @info_b=split/_/,$sub_b[0];
		if($info_a[1] ne $info_b[1] or $info_a[-1] ne $info_b[-1]){
			die;
		}
		$sub_a[0]=$info_a[-1];
		$hash{$info_a[-1]}=join"\t",@sub_a;
	}
	close IN;
	return %hash;
}

	
