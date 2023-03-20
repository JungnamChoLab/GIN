#!/bin/bash

#########################################################################
#                       Global mapping of INsertions                    #
#                                  Ling Wang                            #
#                            lingwang@cemps.ac.cn                       #
#########################################################################

sample=$1       # sample name
te_fam=$2       # TE family name
ltr_end_len=15      # the LTR end length
min_read_len=20        # min read length
max_TSD_len=9           # TSD length -1
reference=/public/home/wangling/tip_rice/ref/msu7      #reference genome
te_dir=/public/home/wangling/tip_rice/ref #

#if [ ! -e ${sample}_*.fastq ]
#then
gunzip *.fastq.gz
#fi

if [ ! -s ${sample}_${te_fam}.fq ]
then
for i in `ls *_1.fastq`;do grep -A2 -B1 -f ${te_dir}/${te_fam}_Prime $i|sed "/^--$/d" >> ${sample}_${te_fam}_1.fq;done
for i in `ls *_2.fastq`;do grep -A2 -B1 -f ${te_dir}/${te_fam}_Prime $i|sed "/^--$/d" >> ${sample}_${te_fam}_2.fq;done
cat ${sample}_${te_fam}_1.fq  ${sample}_${te_fam}_2.fq > ${sample}_${te_fam}.fq
fi
###trim LTR sequences, keep insertion strand
if [ ! -d new ]
then
mkdir new;
fi
if [ ! -s ${sample}_${te_fam}_contain.fq ]
then
cutadapt -a file:${te_dir}/${te_fam}_Prime-start-up.fa --overlap $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_start_up.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cutadapt -a file:${te_dir}/${te_fam}_Prime-end-up.fa --overlap $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_end_up.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cutadapt -g file:${te_dir}/${te_fam}_Prime-start-down.fa --overlap $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_start_down.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cutadapt -g file:${te_dir}/${te_fam}_Prime-end-down.fa --overlap $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_end_down.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cutadapt -a file:${te_dir}/${te_fam}_Prime-up.fa --overlap  $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_up.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cutadapt -g file:${te_dir}/${te_fam}_Prime-down.fa --overlap  $ltr_end_len -m $min_read_len -e 0.01 --untrimmed-output tmp.fq -o ${sample}_${te_fam}_down.fq ${sample}_${te_fam}.fq 2>&1 >> cutadapt_record
cat ${sample}_${te_fam}_up.fq ${sample}_${te_fam}_down.fq > ${sample}_${te_fam}_contain.fq
fi
###
cd new;
bowtie2 -p 15 -a -x $reference -U ../${sample}_${te_fam}_contain.fq -S ${sample}_${te_fam}.sam 2>&1 >> mapping
_record
samtools view -bSF 4 ${sample}_${te_fam}.sam > ${sample}_${te_fam}.bam;samtools sort -o ${sample}_${te_fam}.sorted.bam ${sample}_${te_fam}.bam;samtools index ${sample}_${te_fam}.sorted.bam;samtools view ${sample}_${te_fam}.sorted.bam > ${sample}_${te_fam}.sorted.sam;bamToBed -i ${sample}_${te_fam}.sorted.bam > ${sample}_${te_fam}.bed;bamToBed -i ${sample}_${te_fam}.sorted.bam -ed > ${sample}_${te_fam}_editdistance.bed
###For insertion direction mapping

bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_start_up.fq -S ${sample}_${te_fam}_start_up.sam
samtools view -bSF 4 ${sample}_${te_fam}_start_up.sam > ${sample}_${te_fam}_start_up.bam
samtools sort -o ${sample}_${te_fam}_start_up.sorted.bam ${sample}_${te_fam}_start_up.bam
bamToBed -i ${sample}_${te_fam}_start_up.sorted.bam > ${sample}_${te_fam}_start_up.bed
samtools view ${sample}_${te_fam}_start_up.sorted.bam > ${sample}_${te_fam}_start_up.sam
mv ${sample}_${te_fam}_start_up.sorted.bam ${sample}_${te_fam}_start_up.bam

bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_start_down.fq -S ${sample}_${te_fam}_start_down.sam
samtools view -bSF 4 ${sample}_${te_fam}_start_down.sam > ${sample}_${te_fam}_start_down.bam
samtools sort -o ${sample}_${te_fam}_start_down.sorted.bam ${sample}_${te_fam}_start_down.bam
bamToBed -i ${sample}_${te_fam}_start_down.sorted.bam > ${sample}_${te_fam}_start_down.bed
samtools view ${sample}_${te_fam}_start_down.sorted.bam > ${sample}_${te_fam}_start_down.sam
mv ${sample}_${te_fam}_start_down.sorted.bam ${sample}_${te_fam}_start_down.bam
bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_end_up.fq -S ${sample}_${te_fam}_end_up.sam
samtools view -bSF 4 ${sample}_${te_fam}_end_up.sam > ${sample}_${te_fam}_end_up.bam
samtools sort -o ${sample}_${te_fam}_end_up.sorted.bam ${sample}_${te_fam}_end_up.bam
bamToBed -i ${sample}_${te_fam}_end_up.sorted.bam > ${sample}_${te_fam}_end_up.bed
samtools view ${sample}_${te_fam}_end_up.sorted.bam > ${sample}_${te_fam}_end_up.sam
mv ${sample}_${te_fam}_end_up.sorted.bam ${sample}_${te_fam}_end_up.bam

bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_end_down.fq -S ${sample}_${te_fam}_end_down.sam
samtools view -bSF 4 ${sample}_${te_fam}_end_down.sam > ${sample}_${te_fam}_end_down.bam
samtools sort -o ${sample}_${te_fam}_end_down.sorted.bam ${sample}_${te_fam}_end_down.bam
bamToBed -i ${sample}_${te_fam}_end_down.sorted.bam > ${sample}_${te_fam}_end_down.bed
samtools view ${sample}_${te_fam}_end_down.sorted.bam > ${sample}_${te_fam}_end_down.sam
mv ${sample}_${te_fam}_end_down.sorted.bam ${sample}_${te_fam}_end_down.bam


bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_up.fq -S ${sample}_${te_fam}_up.sam
samtools view -bSF 4 ${sample}_${te_fam}_up.sam > ${sample}_${te_fam}_up.bam
samtools sort -o ${sample}_${te_fam}_up.sorted.bam ${sample}_${te_fam}_up.bam
bamToBed -i ${sample}_${te_fam}_up.sorted.bam > ${sample}_${te_fam}_up_all.bed

awk '$6=="+"' ${sample}_${te_fam}_up_all.bed > ${sample}_${te_fam}_upstream.bed
awk '$6=="-"' ${sample}_${te_fam}_up_all.bed > ${sample}_${te_fam}_downstream.bed

bowtie2 -p 15 -a -x $reference  -U ../${sample}_${te_fam}_down.fq -S ${sample}_${te_fam}_down.sam
samtools view -bSF 4 ${sample}_${te_fam}_down.sam > ${sample}_${te_fam}_down.bam
samtools sort -o ${sample}_${te_fam}_down.sorted.bam ${sample}_${te_fam}_down.bam
bamToBed -i ${sample}_${te_fam}_down.sorted.bam > ${sample}_${te_fam}_down_all.bed

awk '$6=="-"' ${sample}_${te_fam}_down_all.bed >> ${sample}_${te_fam}_upstream.bed
sort -k1,1 -k2,2n ${sample}_${te_fam}_upstream.bed -o ${sample}_${te_fam}_upstream.bed
awk '$6=="+"' ${sample}_${te_fam}_down_all.bed >> ${sample}_${te_fam}_downstream.bed
sort -k1,1 -k2,2n ${sample}_${te_fam}_downstream.bed -o ${sample}_${te_fam}_downstream.bed

sort -k1,1 -k3,3n ${sample}_${te_fam}_upstream.bed |bedtools groupby -i - -g 1,3 -c 3,5 -o count,freqdesc |awk '{print $1"\t"$2-500"\t"$2"\t"$3"\t"$4}' > ${sample}_${te_fam}_upp.bed
sort -k1,1 -k3,3n ${sample}_${te_fam}_downstream.bed |bedtools groupby -i - -g 1,2 -c 2,5 -o count,freqdesc |awk '{print $1"\t"$2+1"\t"$2+501"\t"$3"\t"$4}' > ${sample}_${te_fam}_downn.bed
### filter reads mapping to original region
intersectBed -a ${sample}_${te_fam}_upp.bed -b ${te_dir}/${te_fam}_ltr.gtf -v > ${sample}_upp_${te_fam}_omited.bed
intersectBed -a ${sample}_${te_fam}_downn.bed -b ${te_dir}/${te_fam}_ltr.gtf -v > ${sample}_downn_${te_fam}_omited.bed
### filter reads stack
intersectBed -a ${sample}_upp_${te_fam}_omited.bed -b ${sample}_upp_${te_fam}_omited.bed -wo|cut -f 1,2,3,4,5|uniq -c |sort -k 1|awk '$1<4'|awk -v OFS="\t" '{print $2,$3,$4,$5,$6}'|sort -k1,1 -k2,2n > ${sample}_${te_fam}_up_cutting.bed
intersectBed -a ${sample}_downn_${te_fam}_omited.bed -b ${sample}_downn_${te_fam}_omited.bed -wo|cut -f 1,2,3,4,5 |uniq -c |sort -k 1|awk '$1<4'|awk -v OFS="\t" '{print $2,$3,$4,$5,$6}'|sort -k1,1 -k2,2n > ${sample}_${te_fam}_down_cutting.bed
###TSD get insertion position
intersectBed -a ${sample}_${te_fam}_up_cutting.bed -b ${sample}_${te_fam}_down_cutting.bed -wo |awk -v tsd=$max_TSD_len '$11<=tsd && $11>=2' |awk '{print $1"\t"$3"\t"$7"\t"$4"_"$9"\t"$3-$7+1}' > ${sample}_${te_fam}_insertion.txt
###Add mapping quality and insertion direction
intersectBed -a ${sample}_${te_fam}_up_cutting.bed -b ${sample}_${te_fam}_down_cutting.bed -wo |awk -v tsd=$max_TSD_len '$11<=tsd && $11>=2' | awk '{print $1"\t"$7-1"\t"$3"\t"$4"_"$9"\t"$3-$7+1"\t"$5"_"$10}' > ${sample}_${te_fam}_insertion_mappingquality
awk '{print $1"\t"$3-1"\t"$2"\t"$4"\t"$5}' ${sample}_${te_fam}_insertion.txt | intersectBed -a - -b ${sample}_${te_fam}_start_up.bed -wo |awk '$2==$7||$2==$8||$3==$7||$3==$8' | groupBy -i - -g 1,2,3 -c 11 -o mode|sed 's#\t#_#' |sed 's#\t#_#' > ${te_fam}_insertion_su
awk '{print $1"\t"$3-1"\t"$2"\t"$4"\t"$5}' ${sample}_${te_fam}_insertion.txt | intersectBed -a - -b ${sample}_${te_fam}_start_down.bed -wo |awk '$2==$7||$2==$8||$3==$7||$3==$8' | groupBy -i - -g 1,2,3 -c 11 -o mode|sed 's#\t#_#' |sed 's#\t#_#' > ${te_fam}_insertion_sd
awk '{print $1"\t"$3-1"\t"$2"\t"$4"\t"$5}' ${sample}_${te_fam}_insertion.txt | intersectBed -a - -b ${sample}_${te_fam}_end_down.bed -wo |awk '$2==$7||$2==$8||$3==$7||$3==$8' | groupBy -i - -g 1,2,3 -c 11 -o mode|sed 's#\t#_#' |sed 's#\t#_#' > ${te_fam}_insertion_ed
awk '{print $1"\t"$3-1"\t"$2"\t"$4"\t"$5}' ${sample}_${te_fam}_insertion.txt | intersectBed -a - -b ${sample}_${te_fam}_end_up.bed -wo |awk '$2==$7||$2==$8||$3==$7||$3==$8' | groupBy -i - -g 1,2,3 -c 11 -o mode|sed 's#\t#_#'|sed 's#\t#_#' > ${te_fam}_insertion_eu
sed 's#\t#_#' ${sample}_${te_fam}_insertion_mappingquality |sed 's#\t#_#' > ${te_fam}_insertion
perl ~/bin/merge_total.pl ${te_fam}_insertion_su 1 ${te_fam}_insertion 1 ${te_fam}_insertion_s
perl ~/bin/merge_total.pl ${te_fam}_insertion_sd 1 ${te_fam}_insertion_s 1 ${te_fam}_insertion_sa
perl ~/bin/merge_total.pl ${te_fam}_insertion_eu 1 ${te_fam}_insertion_sa 1 ${te_fam}_insertion_ea
perl ~/bin/merge_total.pl ${te_fam}_insertion_ed 1 ${te_fam}_insertion_ea 1 ${te_fam}_insertion_e
awk '{if ($6=="+"||$12=="+"||$8=="-"||$10=="-"){print$1"\t"$2"\t"$4"\t""+""\t"$3}else{print$1"\t"$2"\t"$4"\t""-""\t"$3}}' ${te_fam}_insertion_e |sed 's#_#\t#' |sed 's#_#\t#'> ${te_fam}_insertion_direction
rm ${te_fam}_insertion ${te_fam}_insertion_s* ${te_fam}_insertion_e*
#rm ${sample}_${te_fam}_end*.b* ${sample}_${te_fam}_end*.sam ${sample}_${te_fam}_start*.b* ${sample}_${te_fam}_start*.sam
###Add TSD sequence
intersectBed -a ${te_fam}_insertion_direction -b ${sample}_${te_fam}_upstream.bed -wo |awk -v OFS="\t" '{print$1,$2,$3,$4,$5,$6,$7,$8"_"$9+1"_"$11}' > ${sample}_${te_fam}_upstream_reads
awk '{print $3"_"$4"_"$1"\t"$10}' ${sample}_${te_fam}.sam > ${sample}_${te_fam}.seq
perl ~/bin/merge_total.pl ${sample}_${te_fam}.seq 1 ${sample}_${te_fam}_upstream_reads 8 tmp
cut -f 1-8,10 tmp |awk '{print $0"\t"substr($9,length($9)-$7+1,$7)}'|groupBy -i - -g 1-7 -c 10 -o mode > ${sample}_${te_fam}_upstream
rm tmp
intersectBed -a ${sample}_${te_fam}_upstream -b ${sample}_${te_fam}_downstream.bed -wo |awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9"_"$10+1"_"$12}'> ${sample}_${te_fam}_downstream_reads
perl ~/bin/merge_total.pl ${sample}_${te_fam}.seq 1 ${sample}_${te_fam}_downstream_reads 9 tmp
cut -f 1-8,11 tmp |awk '{print $0"\t"substr($9,0,$7)}'|groupBy -i - -g 1-8 -c 10 -o mode > ${sample}_${te_fam}_insertion_summary.txt
sort -k1,1 -k2,2n -o ${sample}_${te_fam}_insertion_summary.txt ${sample}_${te_fam}_insertion_summary.txt
rm tmp;rm ${sample}_${te_fam}_insertion_mappingquality;rm ${te_fam}_insertion_direction;
#rm ${sample}_${te_fam}_upstream*;rm ${sample}_${te_fam}_downstream*;rm ${sample}_${te_fam}_down.*am;rm ${sample}_${te_fam}_up.*am;
mv ${sample}_${te_fam}.sorted.bam ${sample}_${te_fam}.bam
mv ${sample}_${te_fam}.sorted.bam.bai ${sample}_${te_fam}.bam.bai
#rm ${sample}_${te_fam}.s*
