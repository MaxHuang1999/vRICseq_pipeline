# -*- coding: utf-8 -*-
# """
# Created on 2023.2.23
#vRICseq data analysis protocol——RNAseq
# @author: Max
#conda 环境迁移：
    #first step:conda list --explicit > spec-list.txt
    #second step:conda creat --name vRICseq --file spec-list.txt
# """

import os,argparse
import sys
import subprocess
#mapping_and_pairs
def vRICpipeline(refseq,fdata,sdata,title):
    os.system("fastp -i %s -I %s -o %s.read1.clean.fq -O %s.read2.clean.fq -z 4 -q 20 -u 30 -f 10"%(fdata,sdata,title,title))#fastp v0.23.2;质控
    os.system("perl /storx/max/workspace/vRICseq/RICpipe-master/step0.remove_PCR_duplicates/scripts/remove_duplicated_reads_PE.pl \
              %s.read1.clean.fq %s.read2.clean.fq %s.read1.clean.rm.fq %s.read2.clean.rm.fq"%(title,title,title,title))#remove_duplicated_reads_PE.pl;去冗余，注意修改路径
    os.system("STAR --runMode genomeGenerate --runThreadN 80 --genomeSAindexNbases 6 \
              --genomeDir align_to_reference_genome --genomeFastaFiles %s"%(refseq))#STAR v2.7.4a;建索引,
    # The parameter ‘--genomeSAindexNbases’ should be set as the minimum between 14 and log2(Genome length)/2 – 1 according to the user manual of the STAR software.（log2(29903)/2 – 1=6 ）
    os.system("STAR --runMode alignReads --genomeDir align_to_reference_genome --readFilesIn  %s.read1.clean.rm.fq --outFileNamePrefix %s_read1_toGenome  --outReadsUnmapped Fastx --outFilterMultimapNmax 100 \
              --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 \
              --limitIObufferSize 1500000000 --runThreadN 80 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --chimOutType SeparateSAMold"%(title,title))#比对
    os.system("STAR --runMode alignReads --genomeDir align_to_reference_genome --readFilesIn  %s.read2.clean.rm.fq --outFileNamePrefix %s_read2_toGenome  --outReadsUnmapped Fastx --outFilterMultimapNmax 100 \
              --outSAMattributes All --alignIntronMin 1 --scoreGapNoncan -4 --scoreGapATAC -4 --chimSegmentMin 15 --chimJunctionOverhangMin 15 --limitOutSJcollapsed 10000000 \
              --limitIObufferSize 1500000000 --runThreadN 80 --alignSJoverhangMin 15 --alignSJDBoverhangMin 10 --alignSJstitchMismatchNmax 5 -1 5 5 --outFilterMatchNminOverLread 0.3 --outFilterScoreMinOverLread 0.3 --chimOutType SeparateSAMold"%(title,title))
    command1='''awk '{if(/^@/){gsub(/\/1$/,"");print}else{print}}' %s.read1.clean.rm.fq > %s.read1.clean.rm.mod.fq'''%(title,title)
    subprocess.call(command1,shell=True)
    command2='''awk '{if(/^@/){gsub(/\/2$/,"");print}else{print}}' %s.read2.clean.rm.fq > %s.read2.clean.rm.mod.fq'''%(title,title)
    subprocess.call(command2, shell=True)#修改fastq文件第一行id号
    os.system("perl ../../select_not_fully_aligned.pl %s.read1.clean.rm.mod.fq %s_read1_toGenomeChimeric.out.sam \
              %s_read1_toGenomeAligned.out.sam > %s_read1_toGenomeUnmapped_really.fq"%(title,title,title,title))
    os.system("perl ../../select_not_fully_aligned.pl %s.read2.clean.rm.mod.fq %s_read2_toGenomeChimeric.out.sam \
                  %s_read2_toGenomeAligned.out.sam > %s_read2_toGenomeUnmapped_really.fq" % (title, title, title, title))#select_not_fully_aligned.pl,挑出未比对上的reads
    os.system("bwa index %s"%(refseq))
    os.system("mkdir z1.second_round_byBWA")
    os.system("bwa mem -t 16 -k 12 -T 15 -o ./z1.second_round_byBWA/read1_futher_by_bwa.sam %s %s_read1_toGenomeUnmapped_really.fq"%(refseq,title))
    os.system("bwa mem -t 16 -k 12 -T 15 -o ./z1.second_round_byBWA/read2_futher_by_bwa.sam %s %s_read2_toGenomeUnmapped_really.fq" % (refseq, title))
    os.system("perl ../../collect_chimeric_ligation_from_sam.pl ./z1.second_round_byBWA/read1_futher_by_bwa.sam 1 ./z1.second_round_byBWA > ./z1.second_round_byBWA/out1.read1.chimeric.sam")
    os.system("perl ../../collect_chimeric_ligation_from_sam.pl ./z1.second_round_byBWA/read2_futher_by_bwa.sam 2 ./z1.second_round_byBWA > ./z1.second_round_byBWA/out1.read2.chimeric.sam")
    #bwa比对结果不会像STAR一样将chimeric reads 挑选出来，需要我们自己写脚本挑选出来

    # #find_all_pairs
    os.system("samtools view -@ 10 -b -S -o %s_read1_toGenomeAligned.out.bam %s_read1_toGenomeAligned.out.sam;samtools view -@ 10 -b -S -o %s_read2_toGenomeAligned.out.bam \
              %s_read2_toGenomeAligned.out.sam"%(title,title, title, title))
    os.system("samtools sort -@ 10 -o %s_read1_toGenomeAligned.out.sort.bam %s_read1_toGenomeAligned.out.bam;\
              samtools sort -@ 10 -o %s_read2_toGenomeAligned.out.sort.bam %s_read2_toGenomeAligned.out.bam"%(title,title, title, title))
    os.system("samtools view -@ 10 -h -q 30 -F 256 -b -o  %s_read1_toGenomeAligned.out.sort.uniq.bam  %s_read1_toGenomeAligned.out.sort.bam;\
              samtools view -@ 10 -h -q 30 -F 256 -b -o  %s_read2_toGenomeAligned.out.sort.uniq.bam  %s_read2_toGenomeAligned.out.sort.bam"%(title,title, title, title))
    os.system("samtools view -h -o %s_read1_toGenomeAligned.out.sort.uniq.sam %s_read1_toGenomeAligned.out.sort.uniq.bam;samtools view -h -o \
              %s_read2_toGenomeAligned.out.sort.uniq.sam %s_read2_toGenomeAligned.out.sort.uniq.bam"%(title,title, title, title))
    os.system("perl ../../precess_Chimeric_sam.pl %s_read1_toGenomeChimeric.out.sam > %s_read1_toGenomeChimeric.out.processed.sam;\
              perl ../../precess_Chimeric_sam.pl %s_read2_toGenomeChimeric.out.sam > %s_read2_toGenomeChimeric.out.processed.sam"%(title,title, title, title))
    os.system("perl ../../obtain_pairs_from_pair.pl %s_read1_toGenomeAligned.out.sort.uniq.sam %s_read2_toGenomeAligned.out.sort.uniq.sam"%(title,title))
    os.system("samtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_1.bam interaction_from_pair_mapped_reads_1.sam; \
              samtools view -@ 10 -b -S -o interaction_from_pair_mapped_reads_2.bam interaction_from_pair_mapped_reads_2.sam")
    os.system("samtools sort -n -@ 10 -o interaction_from_pair_mapped_reads_1.sort.bam interaction_from_pair_mapped_reads_1.bam;\
              samtools sort -n -@ 10 -o interaction_from_pair_mapped_reads_2.sort.bam interaction_from_pair_mapped_reads_2.bam")
    os.system("samtools view -h -o interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_1.sort.bam;\
              samtools view -h -o interaction_from_pair_mapped_reads_2.sort.sam interaction_from_pair_mapped_reads_2.sort.bam")
    os.system("perl ../../obtain_pairs_from_gapped_reads.pl %s_read1_toGenomeAligned.out.sort.uniq.sam %s_read2_toGenomeAligned.out.sort.uniq.sam \
              %s_read1_toGenomeChimeric.out.processed.sam %s_read2_toGenomeChimeric.out.processed.sam interaction_from_gapped_reads.sam"%(title,title, title, title))
    os.system("perl ../../merge_interaction.pl num_of_interactions_from_part.list interaction_from_pair_mapped_reads_1.sort.sam \
              interaction_from_pair_mapped_reads_2.sort.sam interaction_from_gapped_reads.sam  %s_read1_toGenomeChimeric.out.processed.sam \
              %s_read2_toGenomeChimeric.out.processed.sam %s_interaction.sam"%(title,title,title))
    os.system("perl ../../count_link_for_each_kind.pl interaction_from_pair_mapped_reads_1.sort.sam interaction_from_pair_mapped_reads_2.sort.sam \
              num_of_interactions_from_part.list %s_read1_toGenomeChimeric.out.processed.sam %s_read2_toGenomeChimeric.out.processed.sam > num_of_interactions.list"%(title,title))
    os.system("mkdir z2.update_interaction_sam ; mkdir z3.update_alignPair")
    os.system("perl ../../supp_bwaSam_and_update.pl %s_interaction.sam num_of_interactions_from_part.list num_of_interactions.list ./z1.second_round_byBWA/out1.read1.chimeric.sam \
              ./z1.second_round_byBWA/out1.read2.chimeric.sam ./z2.update_interaction_sam > ./z2.update_interaction_sam/mergeBwa.interaction.sam"%(title))
    os.system("perl ../../replace_AlignPair_by_ChimericFragment.pl ./z2.update_interaction_sam/mergeBwa.interaction.sam \
              %s_read1_toGenomeChimeric.out.processed.sam %s_read2_toGenomeChimeric.out.processed.sam \
              ./z1.second_round_byBWA/out1.read1.chimeric.sam ./z1.second_round_byBWA/out1.read2.chimeric.sam > ./z3.update_alignPair/update.mergeBwa.interaction.sam"%(title,title))

    #intra_gene_interaction_chimeric
    os.system("perl ../../creat_gapped_and_Align_pair_bed.VirusOnlyVersion.pl ./z3.update_alignPair/update.mergeBwa.interaction.sam")
    #注意修改perl脚本中的id号
    # command3='''awk -F '\t' 'BEGIN{OFS="\t"}{if(($5-$2)>20) print $1,$2,$3,$4,$5,$6,"SARS",$8,$9,$10 }' ./z3.update_alignPair/update.mergeBwa.interaction.gapped.bed |sort |uniq -c |awk -F ' ' 'BEGIN{OFS="\t"}{if($1>=10) print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }' |awk -F '\t' 'BEGIN{OFS="\t"}{if($1>=10) print $1,$2,$3,$4,$5,$6,"SARS_"NR,$8,$9,$10 }' > %s-totalRNAseq_interaction.junction.bed'''%(title)
    # subprocess.call(command3,shell=True)
    # os.system("python ../1.create_element.py %s-totalRNAseq_interaction.junction.bed > %s-totalRNAseq_interaction.element.bed"%(title,title))
    os.system("touch %s-totalRNAseq_interaction.junction.bed;touch %s-totalRNAseq_interaction.element.bed"%(title,title))
    os.system("bedtools pairtopair -a ./z3.update_alignPair/update.mergeBwa.interaction.gapped.bed -b %s-totalRNAseq_interaction.junction.bed -is > gapped_reads_overlapped_with_exon_junction.list"%(title))
    os.system("perl ../../distance_on_virus.pl ./z3.update_alignPair/update.mergeBwa.interaction.Alignpair_read1.Virus.bed ./z3.update_alignPair/update.mergeBwa.interaction.Alignpair_read2.Virus.bed \
              %s-totalRNAseq_interaction.element.bed > pairreads_distance_in_virus.list"%(title))
    os.system("perl ../../split_intra_to_Normal_Chimeric.VirusOnlyVersion.pl ./z3.update_alignPair/update.mergeBwa.interaction.sam gapped_reads_overlapped_with_exon_junction.list pairreads_distance_in_virus.list 600")
    os.system("perl ../../stat_of_class.pl class_of_AlignPair.log > num_of_class.log")

    #visual_by_JuiceBox
    os.system("perl ../../filter_too_mismatch.pl ./z3.update_alignPair/update.mergeBwa.interaction.Chimeric.sam 1 > %s_interaction.Chimeric.HQ.sam"%(title))
    os.system("perl ../../sam_to_loops_split.pl %s_interaction.Chimeric.HQ.sam NC_019843.3 0 40000 > %s_interaction.loops"%(title,title))#注意修改基因组名称
    command4='''awk -F '\t' 'BEGIN{OFS="\t"}{if(($3<=30119)&&($6<=30119)) print $0}' %s_interaction.loops > %s_interaction.filt.loops'''%(title,title)
    subprocess.call(command4,shell=True)
    os.system("perl ../../convert_for_circos_coor.pl %s_interaction.filt.loops > %s_virus.loops"%(title,title))
    os.system("perl ../../interaction_from_loop_to_JuiceBox.pl %s_virus.loops > %s_in_Virion.list"%(title,title))
    os.system("cat %s_in_Virion.list | sort -k3,3 -k7,7 > %s_in_Virion.sort.list"%(title,title))
    os.system("java -jar /storx/max/workspace/vRICseq/D.visual_by_JuiceBox/scripts/common/juicer_tools.jar pre -r 1,2,5,10,25,50 %s_in_Virion.sort.list %s_in_Virion.hic ../../MERS.chrom.sizes"%(title,title))



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="vRICseq data analysis pipeline-RNAseq. First-writtern by Max.")
    parser.add_argument("-r", "--refseq",type=str,required=True,help="Please type in the file of refseq.")
    parser.add_argument("-f","--fdata",type=str,required=True,help="Please type in the file of first sequencing data")
    parser.add_argument("-s", "--sdata",type=str,required=True,help="Please type in the file of second sequencing data")
    parser.add_argument("-t", "--title", type=str, required=True,help="Please type in the name of outputfile;for example:SARSCOV2-1")
    Args = parser.parse_args()
    vRICpipeline(os.path.abspath(Args.refseq),os.path.abspath(Args.fdata),os.path.abspath(Args.sdata),
                 os.path.abspath(Args.title))



