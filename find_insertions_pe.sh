#!/bin/bash
sample=$1       # sample name
te_fam=$2       # TE family name
fq1=`ls *_1.fastq |xargs| sed 's#\s\+#,#g'`
fq2=`ls *_2.fastq |xargs| sed 's#\s\+#,#g'`
genome=/public/home/wangling/tip_rice/ref/msu7      #reference genome
te_dir=/public/home/wangling/tip_rice/ref
maxdistance=75
echo "${sample} map2 ${te_fam}" >> pe_mapping_record
if [ ! -s ${sample}_map2_${te_fam}.sam ]
then
bowtie2 -p 12 -a --very-sensitive --local -x ${te_dir}/${te_fam} -1 $fq1 -2 $fq2 -S ${sample}_map2_${te_fam}.sam >> pe_mapping_record 2>&1
awk '$2==69 || $2==133 || $2==165 || $2==101 || $2==181 || $2==117'  ${sample}_map2_${te_fam}.sam > ${sample}_${te_fam}.used.sam
awk '$6 ~/S/' ${sample}_map2_${te_fam}.sam > ${sample}_${te_fam}.soft.sam
fi
#samtools faidx ${te_dir}/${te_fam}.fa
samtools view -t ${te_dir}/${te_fam}.fa.fai -bS ${sample}_${te_fam}.used.sam > ${sample}_${te_fam}.used.bam
bamToFastq -i ${sample}_${te_fam}.used.bam -fq ${sample}_${te_fam}.pair.unmap.fq
echo "${sample} ${te_fam} pair of mate map2 genome" >> pe_mapping_record
bowtie2 -p 12 --local -x $genome -U ${sample}_${te_fam}.pair.unmap.fq -S ${sample}_${te_fam}.insertion.sam >> pe_mapping_record 2>&1
#bowtie2 -p 12 -a --local -x $genome -U ${sample}_${te_fam}.pair.unmap.fq -S ${sample}_${te_fam}.insertion.all.sam
samtools view -bS ${sample}_${te_fam}.insertion.sam > ${sample}_${te_fam}.insertion.bam
samtools sort -o ${sample}_${te_fam}.insertion.sorted.bam ${sample}_${te_fam}.insertion.bam
mv ${sample}_${te_fam}.insertion.sorted.bam ${sample}_${te_fam}.insertion.bam
samtools index ${sample}_${te_fam}.insertion.bam
bamToBed -i ${sample}_${te_fam}.insertion.bam > ${sample}_${te_fam}.insertion.bed
mergeBed -s -i ${sample}_${te_fam}.insertion.bed -c 4,5,6 -o count,freqdesc,freqdesc  > ${sample}_${te_fam}.direction.bed
sort -k1,1 -k 2,2n ${sample}_${te_fam}.direction.bed -o ${sample}_${te_fam}.direction.bed
mergeBed -d $maxdistance -i ${sample}_${te_fam}.insertion.bed -c 4,5,6 -o count,freqdesc,freqdesc  > ${sample}_${te_fam}.pe.cloud.bed
awk '$4>2' ${sample}_${te_fam}.pe.cloud.bed > ${sample}_${te_fam}.pe.support3.cloud.bed
samtools view -t ${te_dir}/${te_fam}.fa.fai -bS ${sample}_${te_fam}.soft.sam > ${sample}_${te_fam}.soft.bam
bamToFastq -i ${sample}_${te_fam}.soft.bam -fq ${sample}_${te_fam}.soft.fq
#echo "${sample} ${te_fam} softcliping reads map2 genome" >> pe_mapping_record
#bowtie2 -p 12 -a --local -x $genome -U ${sample}_${te_fam}.soft.fq -S ${sample}_${te_fam}.soft.insertion.sam >> pe_mapping_record 2>&1
#bowtie2 -p 12 -a --local -x $genome -U ${sample}_${te_fam}.pair.unmap.fq -S ${sample}_${te_fam}.soft.insertion.all.sam
#samtools view -bS ${sample}_${te_fam}.soft.insertion.sam > ${sample}_${te_fam}.soft.insertion.bam
#samtools sort -o ${sample}_${te_fam}.soft.insertion.sorted.bam ${sample}_${te_fam}.soft.insertion.bam
#mv ${sample}_${te_fam}.soft.insertion.sorted.bam ${sample}_${te_fam}.soft.insertion.bam
#samtools index ${sample}_${te_fam}.soft.insertion.bam
#bamToBed -i ${sample}_${te_fam}.soft.insertion.bam > ${sample}_${te_fam}.soft.insertion.bed


awk -v OFS="\t" '{print $1,$2-500,$3+500,$4,$5,$6,$7,$8,$9}' new/${sample}_${te_fam}_insertion_summary.txt > ${sample}_${te_fam}_gin_around.bed
awk '{print $0"\t"$3-$2}' ${sample}_${te_fam}.pe.cloud.bed |awk '$7>700'|intersectBed -a - -b ${sample}_${te_fam}.direction.bed -wo |groupBy -i - -g 1-7 -c 13,14,14 -o collapse,collapse,sum |awk '{print $0"\t"$10-$7+1}' |awk '$11<20' > ${sample}_${te_fam}_pe_detected_new_insertion.txt
intersectBed -a ${sample}_${te_fam}.pe.cloud.bed -b ${sample}_${te_fam}_gin_around.bed -wo |awk '{print $0"\t"$3-$2}' |awk '$17 >700' > ${sample}_${te_fam}_gin_detected_new_insertion.txt
intersectBed -a ${sample}_${te_fam}.pe.cloud.bed -b ${te_dir}/${te_fam}_ltr+_1000.gtf -wo| groupBy -i - -g 7-9,12-14 -c 5,6 -o collapse,collapse -delim "_" > ${sample}_${te_fam}_gin_detected_original_insertion.txt
rm ${sample}_${te_fam}.soft.sam ${sample}_map2_${te_fam}.sam ${sample}_${te_fam}.soft.bam
#rm ${sample}_${te_fam}.soft.insertion.sam ${sample}_${te_fam}.soft.sam ${sample}_map2_${te_fam}.sam
