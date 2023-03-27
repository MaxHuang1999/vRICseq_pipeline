#!usr/bin/python3q
import sys,os
samples=(
"3.PEDV-1",
"4.PEDV-2"
)
for i in samples:
        name = i.split(".")[1]
        cmd1 = "mkdir %s/find_all_pairs"%(i) # create folders named as samples in 2.trim
        os.system(cmd1)
        cmd2 = "ln -s /public/home/dgwei/LDH-220324T/vRICseq_data_analysis/7.mapping_and_pairs/3.PEDV-1/%s_read1_toGenome_Aligned.out.sam %s/find_all_pairs"%(name,i)
        os.system(cmd2)
        cmd3 = "cp /public/home/dgwei/LDH-220324T/vRICseq_data_analysis/7.mapping_and_pairs/find_all_pairs/*.pl %s/find_all_pairs"%(i)
        os.system(cmd3)
        cmd4 = "ln -s /public/home/dgwei/LDH-220324T/vRICseq_data_analysis/7.mapping_and_pairs/3.PEDV-1/%s_read1_toGenome_Chimeric.out.sam %s/find_all_pairs"%(name,i)
        os.system(cmd4)
        cmd5 = "ln -s /public/home/dgwei/LDH-220324T/vRICseq_data_analysis/7.mapping_and_pairs/3.PEDV-1/%s_read2_toGenome_Aligned.out.sam %s/find_all_pairs"%(name,i)
        os.system(cmd5)
        cmd6 = "ln -s /public/home/dgwei/LDH-220324T/vRICseq_data_analysis/7.mapping_and_pairs/3.PEDV-1/%s_read2_toGenome_Chimeric.out.sam %s/find_all_pairs"%(name,i)
        os.system(cmd6)
        
        outfile = open(r"%s/find_all_pairs/work.sh"%(i),"w")
#process mapped results
        print("samtools view -@ 10 -b -S -o %s_read1_toGenome_Aligned.out.bam %s_read1_toGenome_Aligned.out.sam"%(name,name),file=outfile)
        print("samtools view -@ 10 -b -S -o %s_read2_toGenome_Aligned.out.bam %s_read2_toGenome_Aligned.out.sam"%(name,name),file=outfile)
        print("samtools sort -@ 10 -o %s_read1_toGenome_Aligned.out.sort.bam %s_read1_toGenome_Aligned.out.bam"%(name,name),file=outfile)
        print("samtools sort -@ 10 -o %s_read2_toGenome_Aligned.out.sort.bam %s_read2_toGenome_Aligned.out.bam"%(name,name),file=outfile) 
        print("samtools view -@ 10 -h -q 30 -F 256 -b -o %s_read1_toGenome_Aligned.out.sort.uniq.bam %s_read1_toGenome_Aligned.out.sort.bam"%(name,name),file=outfile)
        print("samtools view -@ 10 -h -q 30 -F 256 -b -o %s_read2_toGenome_Aligned.out.sort.uniq.bam %s_read2_toGenome_Aligned.out.sort.bam"%(name,name),file=outfile)
        print("samtools view -h -o %s_read1_toGenome_Aligned.out.sort.uniq.sam %s_read1_toGenome_Aligned.out.sort.uniq.bam"%(name,name),file=outfile)
        print("samtools view -h -o %s_read2_toGenome_Aligned.out.sort.uniq.sam %s_read2_toGenome_Aligned.out.sort.uniq.bam"%(name,name),file=outfile)
        print("perl precess_Chimeric_sam.pl %s_read1_toGenome_Chimeric.out.sam > %s_read1_toGenome_Chimeric.out.processed.sam"%(name,name),file=outfile)
        print("perl precess_Chimeric_sam.pl %s_read2_toGenome_Chimeric.out.sam > %s_read2_toGenome_Chimeric.out.processed.sam"%(name,name),file=outfile)
#fin#d pairs
        print("perl obtain_pairs_from_pair.pl %s_read1_toGenome_Aligned.out.sort.uniq.sam %s_read2_toGenome_Aligned.out.sort.uniq.sam"%(name,name),file=outfile)
        print("samtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_1.bam interaction_from_pair_mapped_reads_1.sam\nsamtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_2.bam interaction_from_pair_mapped_reads_2.sam\nsamtools sort -n -@ 10 -o interaction_from_pair_mapped_reads_1.sort.bam interaction_from_pair_mapped_reads_1.bam\nsamtools sort -n -@ 10 -o interaction_from_pair_mapped_reads_2.sort.bam interaction_from_pair_mapped_reads_2.bam\nsamtools view -h -o interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_1.sort.bam\nsamtools view -h -o interaction_from_pair_mapped_reads_2.sort.sam interaction_from_pair_mapped_reads_2.sort.bam",file=outfile)
        print("perl obtain_pairs_from_gapped_reads.pl %s_read1_toGenome_Aligned.out.sort.uniq.sam %s_read2_toGenome_Aligned.out.sort.uniq.sam %s_read1_toGenome_Chimeric.out.processed.sam %s_read2_toGenome_Chimeric.out.processed.sam interaction_from_gapped_reads.sam"%(name,name,name,name),file=outfile)
        print("perl merge_interaction.pl num_of_interactions_from_part.list interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam interaction_from_gapped_reads.sam %s_read1_toGenome_Chimeric.out.processed.sam %s_read2_toGenome_Chimeric.out.processed.sam %s_interaction.sam"%(name,name,name),file=outfile)
        print("perl count_link_for_each_kind.pl interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam num_of_interactions_from_part.list %s_read1_toGenome_Chimeric.out.processed.sam %s_read2_toGenome_Chimeric.out.processed.sam > num_of_interactions.list"%(name,name),file=outfile)
