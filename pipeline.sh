#!/bin/bash
# This script takes fastq files as input and performs
# basic QC steps and then identifies reads characteristic
# of defective viral genomes using STAR. It then 
# tabulates the number of reads supporting individual DVGs
# as well as normalizing their support by the number of reads
# mapped to that gene.
# bash pipeline.sh -f REFERENCE_FASTA -o OUT_DIR INPUT_FILES  
# The pipeline includes the following steps:
#	1) Performs FASTQC on raw fastq files
#	2) Uses kraken with a viral database to identify reads
#	assinged to Influenza A
#	3) Trim	low quality bases with trimmomatic
# 	4) Aligns reads to reference using STAR
#	5) Calls consensus variants with samtools and creates
#	sample-wise reference
#	6) Realigns reads to sample-wise reference using 
#qqq	STAR two-pass mode (BASIC)
#	7) Removes duplicates using Picard
#	8) Identifies junction spanning reads using the filter_reads.py
#		script.
#	9) Tabulates DVG support and normalizes by number of mapped reads in that gene
#	10) Calculates overall as well as junction spanning read depth 
#	Author: Michael A. Martin (michael.martin2@emory.edu)

module load fastqc
module load staraligner

while getopts ":f:o:" opt; do
  case "$opt" in
    f ) ref=$OPTARG      ;;
    o ) out=$OPTARG	;;  
esac
done
shift $(( OPTIND - 1 ))

ref_dir=$(echo $ref | rev | cut -d'/' -f 2-|rev)

# Makes output directories
mkdir -p ${out}/kraken
mkdir -p ${out}/flu_fastq
mkdir -p ${out}/trimmed_fastq
mkdir -p ${out}/fastqc
mkdir -p ${out}/bam_1
mkdir -p ${out}/vcf
mkdir -p ${out}/sample_reference
mkdir -p ${out}/bam_2
mkdir -p ${out}/bam_ddup
mkdir -p ${out}/bam_ddup/metrics
mkdir -p ${out}/bam_junc
mkdir -p ${out}/depth
mkdir -p ${out}/read_length
mkdir -p ${out}/dvg_reads
mkdir -p ${out}/dvg_summary
mkdir -p ${out}/dvg_summary/raw

# Creates genome indexing files
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${ref_dir} \
--genomeFastaFiles ${ref} --genomeSAindexNbases 6
samtools faidx ${ref}

rm -f ${out}/depth/all_average_depth.tsv
rm -f ${out}/read_length/all_average_read_length.tsv

for i in "$@"
do
sample=$(echo $i | rev | cut -d '/' -f 1  | rev | cut -d'.' -f 1)
fq_dir=$(echo $i | rev | cut -d'/' -f 2-|rev)
echo $sample
echo $fq_dir

# --STEP 1: Fast QC--
fastqc -threads 4 -o ${out}/fastqc ${fq_dir}/${sample}.fastq

# --STEP 2: Kraken--
kraken2 --threads 4 --use-names --db \
/pylon5/ms5fpbp/mmartin4/utilities/kraken2_github/viral-db ${fq_dir}/${sample}.fastq \
> ${out}/kraken/${sample}.kraken

# Influenza A reads only
grep -i 'influenza A' ${out}/kraken/${sample}.kraken | \
awk '{print $2}' | \
grep -A 3 -Ff - ${fq_dir}/${sample}.fastq |\
grep -v -- "^--$" > ${out}/flu_fastq/${sample}.flu.fastq

# --STEP 4: Trimmomatic--
java -jar /home/mmartin4/.Trimmomatic-0.38/trimmomatic-0.38.jar \
SE -threads 4 -phred33 ${out}/flu_fastq/${sample}.flu.fastq \
${out}/trimmed_fastq/${sample}.trim.flu.fastq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50

# --STEP 5: STAR round 1--
mkdir -p ${out}/bam_1/${sample}
STAR --runThreadN 4 --genomeDir ${ref_dir} --readFilesIn ${out}/trimmed_fastq/${sample}.trim.flu.fastq \
--outFileNamePrefix ${out}/bam_1/${sample}/${sample}.trim.flu. --outSAMattributes NH HI AS nM jI

samtools view -bS -q 255 ${out}/bam_1/${sample}/${sample}.trim.flu.Aligned.out.sam | \
samtools sort > ${out}/bam_1/${sample}/${sample}.sort.trim.flu.bam 

rm ${out}/bam_1/${sample}/${sample}.trim.flu.Aligned.out.sam

# --STEP 5: Sample-wise Reference --
bcftools mpileup --threads 4 -d 1000 -q 255 -Q 20 -I \
-f ${ref} ${out}/bam_1/${sample}/${sample}.sort.trim.flu.bam |\

bcftools call -mv --threads 4 -Ov | bcftools norm -m- -Ov -o ${out}/vcf/${sample}.sort.trim.flu.vcf 

# Makes a copy of the vcf so we can have both compressed and uncompressed
scp ${out}/vcf/${sample}.sort.trim.flu.vcf ${out}/vcf/${sample}.sort.trim.flu.vcf2

# BGzip vcf 
bgzip ${out}/vcf/${sample}.sort.trim.flu.vcf2

mv ${out}/vcf/${sample}.sort.trim.flu.vcf2.gz ${out}/vcf/${sample}.sort.trim.flu.vcf.gz

tabix ${out}/vcf/${sample}.sort.trim.flu.vcf.gz

# Creates sample-wise consensus fasta
mkdir -p ${out}/sample_reference/${sample}

bcftools consensus -f ${ref} -i '(DP4[2]+DP4[3])/(DP4[0]+DP4[1]) > 1' \
${out}/vcf/${sample}.sort.trim.flu.vcf.gz > ${out}/sample_reference/${sample}/${sample}.fasta

# Creates sample-wise genome indexing files
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${out}/sample_reference/${sample} \
--genomeFastaFiles ${out}/sample_reference/${sample}/${sample}.fasta --genomeSAindexNbases 6

# --STEP 6: STAR round 2--
mkdir -p ${out}/bam_2/${sample}

STAR --runThreadN 4 --twopassMode Basic --genomeDir ${out}/sample_reference/${sample} \
--readFilesIn ${out}/trimmed_fastq/${sample}.trim.flu.fastq \
--outFileNamePrefix ${out}/bam_2/${sample}/${sample}.realign.trim.flu. --outSAMattributes NH HI AS nM jI

# Filters for primary alignments (MQ = 255) only
samtools view -bS -q 255 ${out}/bam_2/${sample}/${sample}.realign.trim.flu.Aligned.out.sam | \
samtools sort > ${out}/bam_2/${sample}/${sample}.sort.realign.trim.flu.bam 

rm ${out}/bam_2/${sample}/${sample}.realign.trim.flu.Aligned.out.sam


# --STEP 7: Picard--
java -jar $PICARD MarkDuplicates REMOVE_DUPLICATES=true \
I=${out}/bam_2/${sample}/${sample}.sort.realign.trim.flu.bam \
O=${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam \
M=${out}/bam_ddup/metrics/${sample}.txt


# Indexes BAM
samtools index ${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam 

# Calculates read depth across genome
samtools depth -aa ${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam \
> ${out}/bam_ddup/${sample}.depth.tsv 

# Tabulates read lengths of all reads aligned to each genes
samtools view ${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam | awk \
'{g[1]="PB2"; \
#g[2]="PB1"; \
#g[3]="PA"; \
#g[4]="HA"; \
#g[5]="NP"; \
#g[6]="NA"; \
#g[7]="M"; \
#g[8]="NS"}{split($3,a,"|"); b[a[1]]+=1; c[a[1]]+=length($10); d[a[1]]+=(length($10)*length($10))}END \
#{OFS=""; for (i in g){print g[i],"\t",b[g[i]],"\t",c[g[i]],"\t",d[g[i]],"\n"}}' | \
sed '/^[[:space:]]*$/d' - > ${out}/read_length/${sample}.read_length.tsv

# --STEP 8: Junction spanning reads--
# Identifies junction spanning reads
# We do this by reading into a python script
# and identifying any reads with:
#	-N in the cigar string
#	-at least 15 total alignment matches
#	-minimum # of consecutive alignment matches >=5
#	-no more than three indels
#	-minimum deletion size of 100bp
cat ${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam | \
python filter_reads.py > ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.sam

echo "junction spanning reads identified"

# Converts juction spanning read SAM to BAM
cat ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.sam | samtools view -bS > \
${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.bam

samtools view -bS ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.sam > \
${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.bam

echo "junction spanning sam converted to bam"

rm ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.sam

# --STEP 9: Tabulate DVG support--
samtools view ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.bam | \
awk '{OFS="\t"; start=split($15,a,","); gene=split($3,b,"|"); print b[1],$1,a[2],a[3]}' | \
awk '{ $5 = $4 - $3; if($5>100)print $1,$2,$3,$4}' | \
sort -k1,1 -k3n -k4n | sed 's/ /\t/g' > ${out}/dvg_reads/${sample}.dvg.tsv

echo "junction spanning reads tabulated"

# Tabulates individual DVG support
awk -F'\t' '{OFS="\t"; print $1,$3,$4}' ${out}/dvg_reads/${sample}.dvg.tsv | \
uniq -c | awk '{OFS="\t";print $2,$3,$4,$1}' > ${out}/dvg_summary/raw/${sample}.dvgCount.tsv

echo "individual dvg support tabulated"

# Changes filenames to subject_day format
sample_id=$(echo ${sample} | cut -f3 -d'_')
subject_day=$(grep ${sample_id} id_subjectday.tsv | \
awk '{print $2}')
subject_day=${subject_day::-1}

subject=$(echo ${subject_day} | cut -d'_' -f 1)

day=$(echo ${subject_day} | cut -d'_' -f 2)

mv ${out}/read_length/${sample}.read_length.tsv ${out}/read_length/${subject_day}.read_length.tsv
mv ${out}/dvg_summary/raw/${sample}.dvgCount.tsv ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv

sed -i $"s/^/${day}\t/g" ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv

sed -i $"s/^/${subject}\t/g" ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv

sed -i $"s/^/${subject_day}\t/g" ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv

echo "file names converted"

# Normalizes DVG support by number of reads in that gene
awk -F'\t' 'FNR==NR{split($1,g,"|"); a[g[1]]+=$3;next}{OFS="\t"; print $0,$7/a[$4]}' \
<(samtools idxstats ${out}/bam_ddup/${sample}.ddup.sort.realign.trim.flu.bam) \
<(cat ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv) \
> ${out}/dvg_summary/raw/${subject_day}.norm.dvgCount.tsv

echo "dvg support normalized"

rm ${out}/dvg_summary/raw/${subject_day}.dvgCount.tsv

# --STEP 10: Calculate read depth--
# Calculates junction spanning read depth
samtools depth -aa ${out}/bam_junc/${sample}.junc.ddup.sort.realign.trim.flu.bam > \
${out}/bam_junc/${sample}.junc.depth.tsv 

echo "junction spanning read depth calcualted"

# Combines read depths
join -j1 -o1.2,1.3,1.4,2.4 \
<(<${out}/bam_ddup/${sample}.depth.tsv awk '{print $1"-"$2" "$0}') \
<(<${out}/bam_junc/${sample}.junc.depth.tsv awk '{print $1"-"$2" "$0}') | \
sed 's/ /\t/g' > ${out}/depth/${sample}.all.depth.tsv

# Changes name to subject_day format
mv ${out}/depth/${sample}.all.depth.tsv ${out}/depth/${subject_day}.all.depth.tsv

# Averages read depths and adds to file of all read depths
awk -v subject_day="$subject_day" \
'{g[1]="WG";
g[2]="PB2"; 
g[3]="PB1";
g[4]="PA";
g[5]="HA";
g[6]="NP";
g[7]="NA";
g[8]="M";
g[9]="NS"}{split($1,a,"|"); \
b["WG"]+=1; c["WG"]+=$3; d["WG"]+=$3*3; e["WG"]+=$3; f["WG"]+=$4*$4; \
b[a[1]]+=1; c[a[1]]+=$3; d[a[1]]+=$3*$3; e[a[1]]+=$4; f[a[1]]+=$4*$4}END\
{OFS=""; for (i in g){print \
subject_day,"\t",g[i], "\t", b[g[i]], "\t", c[g[i]], "\t", d[g[i]], "\t", e[g[i]], "\t", f[g[i]], "\n"}}' \
${out}/depth/${subject_day}.all.depth.tsv | sed '/^[[:space:]]*$/d' - >> \
${out}/depth/all_average_depth.tsv

echo "read depths combined"

# Removes old depth files
rm ${out}/bam_ddup/${sample}.depth.tsv
rm ${out}/bam_junc/${sample}.junc.depth.tsv

# Combines individual read lengths into a single file
awk -v subject_day="$subject_day" \
'{OFS="\t"; print subject_day,$0}' ${out}/read_length/${subject_day}.read_length.tsv \
>> ${out}/read_length/all_average_read_length.tsv

echo "file name changed"
done







