#!/bin/bash

start=$SECONDS
# Requires paired end reads as command line input: forward read and reverse read in fastq format

#echo "started ...... "
bowtie2 -x /home/NGS/tools/hg38/hg38 -1 $1 -2 $2 -S output.sam
echo "bowtie completed"

samtools view output.sam -o output.bam
echo "bam file created"

echo "input bam recieved"
samtools sort output.bam -o output.sorted.bam
echo "sorted bam created"
samtools index output.sorted.bam -o output.sorted.bam.bai &
echo "samtools indexing done"

java -jar /home/NGS/ujjwal/picard/picard.jar MarkDuplicates I=output.sorted.bam O=output_marked_duplicates.bam M=marked_dup_metrics.txt
echo "picard tools completed"

samtools mpileup -B -q 1 -f /home/NGS/tools/hg38/hg38.fa output.sorted.bam > output.mpileup.bam 
echo "mpileup file created"

# Varscan
echo "Calling variants using Varscan .........."

java -jar /home/NGS/tools/varscan/VarScan.v2.4.5.jar mpileup2snp output.mpileup.bam --min-coverage 10 --min-reads2 4 --min-var-freq 0.20 --output-vcf 1 --strand-filter 0 --min-avg-qual 0 --p-value 0.99 > output.varscan.snp.vcf


java -jar /home/NGS/tools/varscan/VarScan.v2.4.5.jar mpileup2indel output.mpileup.bam --min-coverage 10 --min-reads2 4 --min-var-freq 0.20 --output-vcf 1 --strand-filter 0 --min-avg-qual 0 --p-value 0.99 > output.varscan.indel.vcf

echo "Varscan completed"

# Bcftools
echo "Calling variants using Bcftools .........."

bcftools mpileup -Ou -f /home/NGS/tools/hg38/hg38.fa output.sorted.bam | bcftools call -mv -Ob -o output.bcftools.bcf
bcftools view output.bcftools.bcf > output.bcftools.vcf

echo "Bcftools completed"


# VT
echo "Calling variants using VT .........."

/home/NGS/tools/vt/vt discover -b output.sorted.bam -s Systm_Genomics_lab -r /home/NGS/tools/hg38/hg38.fa -i 20 -o output.vt.vcf

echo "VT completed"


# Freebayes 
echo "Calling variants using Freebayes .........."

freebayes -f /home/NGS/tools/hg38/hg38.fa output.sorted.bam > output.freebayes.vcf

echo "Freebayes completed"


echo "Sorting varscan and bcftools chromosome wise"
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' output.varscan.snp.vcf > varscan.sorted.vcf
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' output.bcftools.vcf > bcftools.sorted.vcf
awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' output.vt.vcf > vt.sorted.vcf


# Consensus Calling

echo 'File input recived'
echo 'Count of varscan variants:'
wc -l output.varscan.snp.vcf
echo 'Count of bcftools variants:'
wc -l output.bcftools.vcf
echo 'Count of vt variants:'
wc -l output.vt.vcf
echo "Count of freebayes variants"
wc -l output.freebayes.vcf

mkdir output_files
cd output_files

echo 'Comparing Varscan and bcftools'

awk 'NR == FNR {a[$1 $2];next} $1 $2 in a' ../varscan.sorted.vcf ../bcftools.sorted.vcf > chr_varbcf.vcf
echo 'Number of chr loci same between varscan and bcftools is: '
wc -l chr_varbcf.vcf


awk 'NR == FNR {a[$5];next} $5 in a ' ../varscan.sorted.vcf chr_varbcf.vcf > consensus_varbcf.vcf
echo 'Number of common alt positions shared is: '
wc -l consensus_varbcf.vcf

echo 'Comparing consensus var/bcf output with vt'

awk 'NR == FNR {a[$1 $2];next} $1 $2 in a' ../vt.sorted.vcf consensus_varbcf.vcf > chr_varbcf_vt.vcf
echo 'Number of chr loci same between consensus var/bcf out  and vt is:'
wc -l chr_varbcf_vt.vcf

awk 'NR == FNR {a[$5];next} $5 in a' ../vt.sorted.vcf chr_varbcf_vt.vcf > consensus_varbcf_vt.vcf
echo 'Number of shared alt positions shared is: '
wc -l consensus_varbcf_vt.vcf

echo 'Comparing consensus var/bcf/vt with output of freebayes'

awk 'NR == FNR {a[$1 $2];next} $1 $2 in a' ../output.freebayes.vcf consensus_varbcf_vt.vcf > chr_varbcf_vt_freebayes.vcf
echo 'Number of chr loci same between consensus var/bcf, vt and freebayes is: '
wc -l chr_varbcf_vt_freebayes.vcf

awk 'NR == FNR {a[$5];next} $5 in a' ../output.freebayes.vcf chr_varbcf_vt_freebayes.vcf > final_consensus_variants.vcf
echo 'Final number of consensus variants are: '
wc -l final_consensus_variants.vcf


# Clinvar

echo "Finding valid Clinvar variant id's corresponding to the consensus"

mkdir clinvarResults
cd clinvarResults

awk 'NR == FNR {a[$1 $2 $5];next} (("chr")$1 $2 $5 in a)' ../final_consensus_variants.vcf /home/NGS/tools/clinvar_20230520.vcf > validated_clinvar_variants.vcf


echo "Total number of validated Clinvar entries corresponding to your list of variants is:"
wc -l validated_clinvar_variants.vcf

echo "Downstream Pathgenicity Analysis"

grep -w 'CLNSIG' validated_clinvar_variants.vcf > clnid_available.vcf

grep -w 'CLNSIG=Pathogenic' clnid_available.vcf > pathogenic_variants.vcf
echo "Number of pathogenic variants"
wc -l pathogenic_variants.vcf

grep -w 'CLNSIG=Likely_pathogenic' clnid_available.vcf > likely_pathogenic_variants.vcf
echo "Number of likely pathogenic variants"
wc -l likely_pathogenic_variants.vcf

grep -w 'CLNSIG=Pathogenic/Likely_pathogenic' clnid.vcf > pathogenic_and_likely_pathogenic.vcf
echo "Number of pathogenic_likely_pathogenic variants"
wc  -l pathogenic_and_likely_pathogenic.vcf

grep -w 'CLNSIG=Uncertain_significance' clnid.vcf > uncertain_significance.vcf
echo "Number of variants of uncertain significance"
wc  -l uncertain_significance_variants.vcf

grep -w 'CLNSIG=Conflicting_interpretations_of_pathogenicity' clnid.vcf > conflicting_variants.vcf
echo "Number of variants of uncertain significance"
wc  -l conflicting_variants.vcf


echo 'Congratulations! Pipeline ran successfully ....'

end=$SECONDS
timetaken=$(($end-$start))

echo "Elapsed time: $((($timetaken / 60) % 60))min $(($timetaken% 60))sec"


