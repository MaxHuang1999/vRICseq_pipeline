#!/usr/bin/perl
die "perl $0 find_all_pairs/interaction.sam find_all_pairs/num_of_interactions_from_part.list find_all_pairs/num_of_interactions.list bwa/out1.read1.chimeric.sam bwa/out1.read2.chimeric.sam output_dir\n" if(@ARGV != 6);
my $star_pet_sam=shift;
my $num_part_pet=shift;
my $num_all_pet=shift;
my $bwa_read1_pet_sam=shift;
my $bwa_read2_pet_sam=shift;
my $output_dir=shift;
$output_dir=~s/\/$//;

my %gapped_frag_already_by_star;
open(STARPET,$star_pet_sam) || die;
while(my $line=<STARPET>){
	if($line!~/Part/){
		next;
	}
	my @sub=split/\s+/,$line;
	my $frag_id=(split/_/,$sub[0])[-1];
	$gapped_frag_already_by_star{$frag_id}=1;
}
close STARPET;

my $bwa_read1_part_pets;
my $bwa_read1_chimeric_pets;
open(BRPSA,$bwa_read1_pet_sam) || die;
while(my $frag_a=<BRPSA>){
	my $frag_b=<BRPSA>;
	if($frag_a=~/^Part/){
		my @sub_a=split/\s+/,$frag_a;
		my @sub_b=split/\s+/,$frag_b;
		my $id_a=(split/_/,$sub_a[0])[-1];
		my $id_b=(split/_/,$sub_b[0])[-1];
		if($id_a ne $id_b){
			die;
		}
		if(exists $gapped_frag_already_by_star{$id_a}){
			next;
		}
		else{
			$bwa_read1_part_pets.=$frag_a.$frag_b;	#only store bwa specific gapped reads
		}
	}
	elsif($frag_a=~/^Chimeric/){
		$bwa_read1_chimeric_pets.=$frag_a.$frag_b;
	}
}

my $bwa_read2_part_pets;
my $bwa_read2_chimeric_pets;
open(BRPSB,$bwa_read2_pet_sam) || die;
while(my $frag_a=<BRPSB>){
	my $frag_b=<BRPSB>;
	if($frag_a=~/^Part/){
		my @sub_a=split/\s+/,$frag_a;
		my @sub_b=split/\s+/,$frag_b;
		my $id_a=(split/_/,$sub_a[0])[-1];
		my $id_b=(split/_/,$sub_b[0])[-1];
		if($id_a ne $id_b){
			die;
		}
		if(exists $gapped_frag_already_by_star{$id_a}){
			next;
		}
		else{
			$bwa_read2_part_pets.=$frag_a.$frag_b;	#only store bwa specific gapped reads
		}
	}
	elsif($frag_a=~/^Chimeric/){
		$bwa_read2_chimeric_pets.=$frag_a.$frag_b;
	}
}

my $num_paired_pet;
my $num_gapped_pet;
my $num_chimeric_read1_pet;
my $num_chimeric_read2_pet;
open(NAP,$num_all_pet) || die;
while(my $line=<NAP>){
	chomp $line;
	if($line=~s/pair\s+//){
		$num_paired_pet=$line;
	}
	elsif($line=~s/gapped\s+//){
		$num_gapped_pet=$line;
	}
	elsif($line=~s/C1\s+//){
		$num_chimeric_read1_pet=$line;
	}
	elsif($line=~s/C2\s+//){
		$num_chimeric_read2_pet=$line;
	}
	else{
		die;
	}
}

my $num_part_from_alignR1;
my $num_part_from_alignR2;
my $num_part_from_C1;
my $num_part_from_C2;
open(NPP,$num_part_pet) || die;
while(my $line=<NPP>){
	chomp $line;
	if($line=~s/Part_from_Align_Read1:\s+//){
		$num_part_from_alignR1=$line;
	}
	elsif($line=~s/Part_from_Align_Read2:\s+//){
		$num_part_from_alignR2=$line;
	}
	elsif($line=~s/Part_from_Chimeric_Read1:\s+//){
		#$num_part_from_C1=$line;
		$num_part_from_C1=0;
	}
	elsif($line=~s/Part_from_Chimeric_Read2:\s+//){
		#$num_part_from_C2=$line;
		$num_part_from_C1=0;
	}
	else{
		die;
	}
}

my $index_after_part_from_chimeric_read1=$num_paired_pet+$num_part_from_alignR1+$num_part_from_alignR2+$num_part_from_C1;
my $index_after_part_from_chimeric_read2=$num_paired_pet+$num_part_from_alignR1+$num_part_from_alignR2+$num_part_from_C1+$num_part_from_C2;
my $index_after_chimeric_read1=$num_paired_pet+$num_gapped_pet+$num_chimeric_read1_pet;
my $index_after_chimeric_read2=$num_paired_pet+$num_gapped_pet+$num_chimeric_read1_pet+$num_chimeric_read2_pet;

my $num_paired_pet_updated=$num_paired_pet;
my $num_gapped_pet_updated=$num_gapped_pet;
my $num_chimeric_read1_pet_updated=$num_chimeric_read1_pet;
my $num_chimeric_read2_pet_updated=$num_chimeric_read2_pet;

my $num_part_from_alignR1_updated=$num_part_from_alignR1;
my $num_part_from_alignR2_updated=$num_part_from_alignR2;
my $num_part_from_C1_updated=$num_part_from_C1;
my $num_part_from_C2_updated=$num_part_from_C2;

my $pet_index;
my $updated_index;
my $sam_header;
my $print_sam_header=1;
open(SPS,$star_pet_sam) || die;
while(my $frag_a=<SPS>){
	if($frag_a=~/^@/){
		$sam_header.=$frag_a;
	}
	else{
		if($print_sam_header){
			print $sam_header;
			$print_sam_header=0;
		}
		$frag_b=<SPS>;
		$pet_index++;
		$updated_index++;
		
		print_updated_line($frag_a);
		print_updated_line($frag_b);

		if($pet_index == $index_after_part_from_chimeric_read1){
			my @arr_bwa_read1_part_pets=split/\n/,$bwa_read1_part_pets;
			if($#arr_bwa_read1_part_pets < 1){	#none
				next;
			}			
			foreach my $i (0..int($#arr_bwa_read1_part_pets/2)){
				$updated_index++;
				my $arr_index=2*$i;
				#print "New:";
				print_updated_line($arr_bwa_read1_part_pets[$arr_index]);
				#print "New:";
				print_updated_line($arr_bwa_read1_part_pets[$arr_index+1]);
				#last;
			}
			$num_part_from_C1_updated+=int($#arr_bwa_read1_part_pets/2)+1;
			$num_gapped_pet_updated+=int($#arr_bwa_read1_part_pets/2)+1;
		}
		if($pet_index == $index_after_part_from_chimeric_read2){
			my @arr_bwa_read2_part_pets=split/\n/,$bwa_read2_part_pets;
			if($#arr_bwa_read2_part_pets < 1){      #none
				next;
			}
			foreach my $i (0..int($#arr_bwa_read2_part_pets/2)){
				$updated_index++;
				my $arr_index=2*$i;
				#print "New:";
				print_updated_line($arr_bwa_read2_part_pets[$arr_index]);
				#print "New:";
				print_updated_line($arr_bwa_read2_part_pets[$arr_index+1]);
				#last;
			}
			$num_part_from_C2_updated+=int($#arr_bwa_read2_part_pets/2)+1;
			$num_gapped_pet_updated+=int($#arr_bwa_read2_part_pets/2)+1;
		}
		if($pet_index == $index_after_chimeric_read1){
			my @arr_bwa_read1_chimeric_pets=split/\n/,$bwa_read1_chimeric_pets;
			if($#arr_bwa_read1_chimeric_pets < 1){	#none
				next;
			}
			foreach my $i (0..int($#arr_bwa_read1_chimeric_pets/2)){
				$updated_index++;
				my $arr_index=2*$i;
				#print "New:";
				print_updated_line($arr_bwa_read1_chimeric_pets[$arr_index]);
				#print "New:";
				print_updated_line($arr_bwa_read1_chimeric_pets[$arr_index+1]);
				#last;
			}
			$num_chimeric_read1_pet_updated+=int($#arr_bwa_read1_chimeric_pets/2)+1;
		}
		if($pet_index == $index_after_chimeric_read2){
			my @arr_bwa_read2_chimeric_pets=split/\n/,$bwa_read2_chimeric_pets;
			if($#arr_bwa_read2_chimeric_pets < 1){	#none
				next;
			}
			foreach my $i (0..int($#arr_bwa_read2_chimeric_pets/2)){
				$updated_index++;
				my $arr_index=2*$i;
				#print "New:";
				print_updated_line($arr_bwa_read2_chimeric_pets[$arr_index]);
				#print "New:";
				print_updated_line($arr_bwa_read2_chimeric_pets[$arr_index+1]);
				#last;
			}
			$num_chimeric_read2_pet_updated+=int($#arr_bwa_read2_chimeric_pets/2)+1;
		}
	}	
}

#update num_of_interactions.list
open(OUTNI,">$output_dir/mergeBwa.num_of_interactions.list") || die;
print OUTNI "pair\t$num_paired_pet_updated\n";
print OUTNI "gapped\t$num_gapped_pet_updated\n";
print OUTNI "C1\t$num_chimeric_read1_pet_updated\n";
print OUTNI "C2\t$num_chimeric_read2_pet_updated\n";

#update num_of_interactions_from_part.list
open(OUTNIP,">$output_dir/mergeBwa.num_of_interactions_from_part.list") || die;
print OUTNIP "Part_from_Align_Read1:\t$num_part_from_alignR1_updated\n";
print OUTNIP "Part_from_Align_Read2:\t$num_part_from_alignR2_updated\n";
print OUTNIP "Part_from_Chimeric_Read1:\t$num_part_from_C1_updated\n";
print OUTNIP "Part_from_Chimeric_Read2:\t$num_part_from_C2_updated\n";

sub print_updated_line{
	my $sam_line=shift;
	my @sub=split/\s+/,$sam_line;
	my @info=split/_/,$sub[0];
	$info[1]=$updated_index;
	$sub[0]=join"_",@info;
	my $new_line=join"\t",@sub;
	print $new_line,"\n";
}


