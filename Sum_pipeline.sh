#!/bin/bash

set -ex # Enable error detection and command audit trail


CELL=$1

PROBAND_SAMPLE=$2
PROBAND_FASTQ1=$3
PROBAND_FASTQ2=$4

MOTHER_SAMPLE=$5
MOTHER_FASTQ1=$6
MOTHER_FASTQ2=$7

FATHER_SAMPLE=$8
FATHER_FASTQ1=$9
FATHER_FASTQ2=${10}

# Function to process a single sample
process_sample() {
    SAMPLE=$1
    CELL=$2
    FASTQ1=$3
    FASTQ2=$4

    # Check if the FASTQ file exists
    if ! test -e $FASTQ1 ; then
        echo "FASTQ file does not exist"
        exit
    fi

    # Alignment via BWA-MEM
    if test $FASTQ1 -nt temp_stage1_${SAMPLE}.sam && test $FASTQ1 -nt ${SAMPLE}.bam ; then
        /hpdm046/software/bwa/bwa mem -t 4 -M \
        -R "@RG\tID:${CELL}\tPL:ILLUMINA\tSM:${SAMPLE}\tLB:${SAMPLE}_${CELL}" \
        /hpdm046/resources/reference/human_g1k_v37.fasta $FASTQ1 $FASTQ2 \
        >temp_partial_stage1_${SAMPLE}.sam 2>log_${SAMPLE}_stage1_bwa.log
        mv temp_partial_stage1_${SAMPLE}.sam temp_stage1_${SAMPLE}.sam
    else
        echo "Skipping BWA MEM, as it has already been run"
    fi

    # Converting SAM file into BAM
    if test temp_stage1_${SAMPLE}.sam -nt temp_stage2_${SAMPLE}.bam && test temp_stage1_${SAMPLE}.sam -nt ${SAMPLE}.bam; then
        java -jar /hpdm046/software/picard/picard.jar FixMateInformation \
        I=temp_stage1_${SAMPLE}.sam \
        O=temp_partial_stage2_${SAMPLE}.bam VALIDATION_STRINGENCY=SILENT \
        COMPRESSION_LEVEL=0 CREATE_INDEX=true SORT_ORDER=coordinate TMP_DIR=. \
        2>log_${SAMPLE}_stage2_FixMateInformation.log
        mv temp_partial_stage2_${SAMPLE}.bam temp_stage2_${SAMPLE}.bam
        mv temp_partial_stage2_${SAMPLE}.bai temp_stage2_${SAMPLE}.bai
        rm temp_stage1_${SAMPLE}.sam
    fi

    # Removing PCR duplicates
    if test temp_stage2_${SAMPLE}.bam -nt temp_stage3_${SAMPLE}.bam && test temp_stage2_${SAMPLE}.bam -nt ${SAMPLE}.bam ; then
        java -jar /hpdm046/software/picard/picard.jar MarkDuplicates \
        I=temp_stage2_${SAMPLE}.bam \
        O=temp_partial_stage3_${SAMPLE}.bam \
        METRICS_FILE=temp_stage3_${SAMPLE}.duplicates \
        REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT \
        COMPRESSION_LEVEL=0 CREATE_INDEX=true TMP_DIR=. \
        2>log_${SAMPLE}_stage3_MarkDuplicates.log
        mv temp_partial_stage3_${SAMPLE}.bam temp_stage3_${SAMPLE}.bam
        mv temp_partial_stage3_${SAMPLE}.bai temp_stage3_${SAMPLE}.bai
        rm temp_stage2_${SAMPLE}.bam temp_stage2_${SAMPLE}.bai
    fi

    # Identifying Indels
    if test temp_stage3_${SAMPLE}.bam -nt temp_stage4_${SAMPLE}.intervals && test temp_stage3_${SAMPLE}.bam -nt ${SAMPLE}.bam; then
        java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T RealignerTargetCreator \
        -I temp_stage3_${SAMPLE}.bam -o temp_partial_stage4_${SAMPLE}.intervals -known \
        /hpdm046/resources/filtering_annotation/Mills_and_1000G_gold_standard.indels.b37.vcf \
        -R /hpdm046/resources/reference/human_g1k_v37.fasta \
        2>log_${SAMPLE}_stage4_RealignerTargetCreator.log
        mv temp_partial_stage4_${SAMPLE}.intervals temp_stage4_${SAMPLE}.intervals
    fi

    if test temp_stage4_${SAMPLE}.intervals -nt ${SAMPLE}.bam ; then
        java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T IndelRealigner \
        -I temp_stage3_${SAMPLE}.bam \
        -o ${SAMPLE}.bam -targetIntervals temp_stage4_${SAMPLE}.intervals -known \
        /hpdm046/resources/filtering_annotation/Mills_and_1000G_gold_standard.indels.b37.vcf \
        -R /hpdm046/resources/reference/human_g1k_v37.fasta \
        2>log_${SAMPLE}_stage5_IndelRealigner.log
        rm temp_stage3_${SAMPLE}.bam temp_stage3_${SAMPLE}.bai temp_stage4_${SAMPLE}.intervals
    fi

    # Base Quality Score Recalibration (BQSR)
    if test ${SAMPLE}.bam -nt ${SAMPLE}_bqsr.grp ; then
        java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T BaseRecalibrator \
        -I ${SAMPLE}.bam \
        -R /hpdm046/resources/reference/human_g1k_v37.fasta \
        -knownSites /hpdm046/resources/filtering_annotation/ExAC.r1.sites.vep.vcf.gz \
        -o ${SAMPLE}_bqsr.grp \
        2>log_${SAMPLE}_stage6_BaseRecalibrator
    fi

    # Variant Calling via HaplotypeCaller
    if test ${SAMPLE}_bqsr.grp -nt ${SAMPLE}.g.vcf ; then
        java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T HaplotypeCaller \
        -I ${SAMPLE}.bam --BQSR ${SAMPLE}_bqsr.grp \
        -o ${SAMPLE}.g.vcf --emitRefConfidence GVCF \
        -L v501_WGS.interval_list \
        -R /hpdm046/resources/reference/human_g1k_v37.fasta -stand_call_conf 10.0 \
        2>log_${SAMPLE}_stage7_HaplotypeCaller
    fi
}

# Process each sample in the trio by calling the process_sample function for each set of input parameters
process_sample $PROBAND_SAMPLE $CELL $PROBAND_FASTQ1 $PROBAND_FASTQ2
process_sample $MOTHER_SAMPLE $CELL $MOTHER_FASTQ1 $MOTHER_FASTQ2
process_sample $FATHER_SAMPLE $CELL $FATHER_FASTQ1 $FATHER_FASTQ2

# Joint genotyping
if [[ ${PROBAND_SAMPLE}.g.vcf -nt Family.g.vcf || ${MOTHER_SAMPLE}.g.vcf -nt Family.g.vcf || ${FATHER_SAMPLE}.g.vcf -nt Family.g.vcf ]]; then
    java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T CombineGVCFs \
        -R /hpdm046/resources/reference/human_g1k_v37.fasta \
        --variant ${PROBAND_SAMPLE}.g.vcf \
        --variant ${MOTHER_SAMPLE}.g.vcf \
        --variant ${FATHER_SAMPLE}.g.vcf \
        -o Family.g.vcf \
        2>log_Family_stage8_CombineGVCFs.log
fi

# Convert to VCF

if test Family.g.vcf -nt Family.vcf ; then

	java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T GenotypeGVCFs \
		-R /hpdm046/resources/reference/human_g1k_v37.fasta \
		-V Family.g.vcf \
		-o Family.vcf \
		2>log_Family_stage9_GenotypeGVCFs
fi

# Remove dbSNP Common and Artefacts

if test Family.vcf -nt Family_filtered.vcf ; then

	java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T SelectVariants \
		-R /hpdm046/resources/reference/human_g1k_v37.fasta \
		-V Family.vcf \
		-o temp_Family.vcf --discordance \
		/hpdm046/resources/filtering_annotation/common_no_known_medical_impact_20160302-edited.vcf\
		2>log_Family_stage10_SelectVariants_dbSNP.log

	java -jar /hpdm046/software/gatk/gatk-3.6.0/GenomeAnalysisTK.jar -T SelectVariants \
		-R /hpdm046/resources/reference/human_g1k_v37.fasta \
                -V temp_Family.vcf \
                -o Family_filtered.vcf --discordance \
                /hpdm046/resources/filtering_annotation/V501_19052014_common_artefacts.vcf \
		2>log_Family_stage10_SelectVariants_Artefacts.log
fi


if test Family_filtered.vcf -nt Alamut_Family.txt ; then

	/hpdm046/software/alamut/alamut-batch-standalone-1.11/alamut-batch --in Family_filtered.vcf \
	--ann Alamut_Family.txt --unann Alamut_Family_unnanotated.txt --ssIntronicRange 2 \
	--outputVCFInfo VariantType AC AF AN DP FS MQ QD --outputVCFGenotypeData GT AD DP GQ PL \
	--outputVCFQuality --outputVCFFIlter \
	2>log_FINAL_stage13_Alamut.log
fi



