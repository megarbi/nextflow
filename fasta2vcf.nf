#!/usr/bin/env nextflow

nextflow.enable.dsl=2

VERSION="0.0"

log.info "===================================================================="
log.info "GATK4 Best Practice Nextflow Pipeline (v${VERSION})                        "
log.info "===================================================================="

params.help = ""
if (params.help) {
  log.info " "
  log.info "USAGE: "
  log.info " "
  log.info "nextflow run my_main.nf --fastq1 read_R1.fastq.gz --fastq2 read_R2.fastq.gz"
  log.info " "
  log.info "Mandatory arguments:"
  log.info "    --fastq1        FILE               Fastq(.gz) file for read1"
  log.info "    --fastq2        FILE               Fastq(.gz) file for read2"
  log.info " "
  log.info "Optional arguments:"
  log.info "    --outdir        DIR                Output directory(default: ./Results)"
  log.info "    --samplename    STRING             Sample name(dafault: fastq1 basename)"
//  log.info "    --rg            STRING             Read group tag(dafault: fastq1 basename)"
  log.info " "
  log.info "===================================================================="
  exit 1
}


//params.out = "${params.outdir}/out"
//arams.tmpdir = "${params.outdir}/tmp"

// only works with paired ends

//  The default workflow
workflow {

    //  Input data is received through channels
    split_read(params.fastq1,params.fastq2)
    //split_read(params.fastq2)
    
    // quality control pre trimming
    fastqc_pre(split_read.out.split1,split_read.out.split2)
    multiqc_pre(fastqc_pre.out)

    // Trim reads
    trimmomatic(split_read.out.split1,split_read.out.split2)

    // quality control post trimming
    fastqc_post(trimmomatic.out.trim1,trimmomatic.out.trim2)
    multiqc_post(fastqc_post.out)
    
    // BWA
    bwa(trimmomatic.out.trim1,trimmomatic.out.trim2)
    
    // mark duplicates
    mark_duplicates(bwa.out)
    //NmMdAndUqTags

    sort_and_tag(mark_duplicates.out)
    // metrics
    alignment_metrics(sort_and_tag.out)
    
    //baserecalibrator
    base_recalibrator(sort_and_tag.out)
    
    //HaplotypeCaller
    haplotype_caller_vcf(base_recalibrator.out)
    
    
}


process multiqc_pre {


    //docker: gatk.4.0.3
    
    publishDir "${params.sample_ID}_PRE_logs"
    
    input:
    file("${params.sample_ID}_PRE_logs")

    //output:
    //file("${params.sample_ID}_PRE_multiqc")

    script:
    """
    multiqc -d ${params.sample_ID}_PRE_logs/ -o ${params.sample_ID}_PRE_logs/multiqc_R1 -x "*.R2.*" -f -q  &
    multiqc -d ${params.sample_ID}_PRE_logs/ -o ${params.sample_ID}_PRE_logs/multiqc_R2 -x "*.R1.*" -f 
    """
}

process multiqc_post {

    publishDir "${params.sample_ID}_POST_logs"
    
    input:
    file("${params.sample_ID}_POST_logs")

    //output:
    //file("${params.sample_ID}_POST_multiqc")

    script:
    """
    multiqc -d ${params.sample_ID}_POST_logs/ -o ${params.sample_ID}_POST_logs/multiqc_R1 -x "*._2P*" -f -q  &
    multiqc -d ${params.sample_ID}_POST_logs/ -o ${params.sample_ID}_POST_logs/multiqc_R2 -x "*._1P*" -f 
    """
}

process fastqc_pre{

    publishDir params.fastqcDir
    
    input:
    //tuple val(prefix), path(reads)
    path(split1)
    path(split2)
    
    output:
    file("${params.sample_ID}_PRE_logs")
    

    script:
    """
    mkdir ${params.sample_ID}_PRE_logs
    fastqc -o ${params.sample_ID}_PRE_logs -f fastq -q ${split1} -t 2 & 
	fastqc -o ${params.sample_ID}_PRE_logs -f fastq -q ${split2} -t 2 
    """  
}


process fastqc_post{

    publishDir params.fastqcDir
    
    input:
    //tuple val(prefix), path(reads)
    path(trim1)
    path(trim2)
    
    output:
    file("${params.sample_ID}_POST_logs")
    

    script:
    """
    mkdir ${params.sample_ID}_POST_logs
    fastqc -o ${params.sample_ID}_POST_logs -f fastq -q ${trim1} -t 2 & 
	fastqc -o ${params.sample_ID}_POST_logs -f fastq -q ${trim2} -t 2 
    """  
}


process split_read {

    publishDir params.splitDir

    input:
    //file(reads)
    path(fastq1)
    path(fastq2)

    output:
    path("${params.sample_ID}.*.L*.R1.fastq"), emit: split1
    path("${params.sample_ID}.*.L*.R2.fastq"), emit: split2
    //path("${params.sample_ID}.*.L*.R2.fastq"), emit: fnames2
    //tuple path("${params.sample_ID}.*.L*.R1.fastq"), path("${params.sample_ID}.*.L*.R2.fastq")
    
    script:
    """
    
    if [[ "${fastq1}" = *.gz ]]; then
    CAT=zcat 
    elif [[ "${fastq1}" = *.fastq ]]; then
    CAT=cat
    else
    	printf "Not a fastq or gzip file, exiting."
    false
    fi
    
    \${CAT}  < ${fastq1} | awk 'BEGIN {FS = ":| "} {
    fc=\$3; lane=\$4; 
    name="${params.sample_ID}."fc".L"lane".R1.fastq"; 
    print > name;  
    for (i = 1; i <= 3; i++) {getline; print > name }}'  & 


    if [[ "${fastq2}" = *.gz ]]; then
    CAT=zcat
    elif [[ "${fastq2}" = *.fastq ]]; then
    CAT=cat
    else
    	printf "Not a fastq or gzip file, exiting."
    false
    fi
    
    \${CAT} < ${fastq2} | awk 'BEGIN {FS = ":| "} {
    fc=\$3; lane=\$4; 
    name="${params.sample_ID}."fc".L"lane".R2.fastq"; 
    print > name;  
    for (i = 1; i <= 3; i++) {getline; print > name }}'    
  
    """
    
}


process trimmomatic {

    publishDir params.trimDir

    input:
    path(split_reads1)
    path(split_reads2)
    
    output:
    //path()
    //tuple path(fq_1_paired), path(fq_2_paired) into ch_out_trimmomatic
    path("${params.sample_ID}*._1P"), emit: trim1
    path("${params.sample_ID}*._2P"), emit: trim2
    //stdout
    
    script:

    // Remove Illumina adapters provided in the TruSeq3-PE.fa file (provided). Initially
    // Trimmomatic will look for seed matches (16 bases) allowing maximally 2
    // mismatches. These seeds will be extended and clipped if in the case of paired end
    // reads a score of 30 is reached (about 50 bases), or in the case of single ended reads a
    // score of 10
    // Headcrop: Removes the specified number of bases, regardless of quality, from the beginning of the read.
    // Minlen:   Removes reads that fall below the specified minimal length
    
    """
    #echo "AAA"
    #ls $split_reads1 $split_reads2 

    for r1 in ${split_reads1}
    do
    if [[ \$r1 =~ ".R1." ]]
    then
    r2=\${r1//".R1."/".R2."}
    basename=`echo \$r1 | awk -F"/" '{print \$NF}' | awk -F"R1.fastq" '{print \$1}'`
    trimmomatic PE -threads 2 \${r1} \${r2} \
	ILLUMINACLIP:$params.databaseDir/TruSeq3-PE.fa:2:30:10 HEADCROP:3 MINLEN:35 -baseout \${basename} &
	fi
    done
    wait
    """
    
}

//    echo trimmomatic PE -threads 2 \
//	${split_reads1} ${split_reads2} \
//	ILLUMINACLIP:$params.databaseDir/TruSeq3-PE.fa:2:30:10 \
//	HEADCROP:3 MINLEN:35


//
process bwa {

    input:
    path(trimmed_reads1)
    path(trimmed_reads2)
    //file reference
    //file bwa_index

    output:
    path("${params.sample_ID}/mapped_bam/")
    
    script:
    """
    mkdir -p ${params.sample_ID}/mapped_bam
    for r1 in ${trimmed_reads1}

    do

    if [[ \$r1 =~ "._1P" ]]
    then
    r2=\${r1//"_1P"/"_2P"}


    basename=`echo \$r1 | awk -F"/" '{print \$NF}' | awk -F"_1P" '{print \$1}'`
    fc_name=`echo \$basename | awk -F"." '{print \$(NF-2)}'`
    lane=`echo \$basename | awk -F"." '{print \$(NF-1)}'`
    id="\$fc_name.\${lane}"
    echo \${basename} \${fc_name} \$lane \${id}
    bwa mem -M -t 4 -R "@RG\\tID:\${id}\\tPL:Illumina\\tSM:${params.sample_ID}\\tLB:${params.sample_ID}" \
	${params.reference_genome} \${r1} \${r2} | \
	picard SortSam -I /dev/stdin -SORT_ORDER queryname \
	-O ${params.sample_ID}/mapped_bam/\${id}_map.bam 
    
    fi 
    done
    """

}


process mark_duplicates {

    input:
    path(mapped_bams)
    
    // make list with all mapped bams and feed the list to MarkDuplicates
    output:
    file("${params.sample_ID}/dedup_bam/${params.sample_ID}_merged_dedup.bam")
    

    script:
    """
    mkdir -p ${params.sample_ID}/dedup_bam

    ls -latr ${mapped_bams}/*.bam | awk '{print "-I "\$9}' > ${params.sample_ID}_mapped_bams.list.txt

    gatk MarkDuplicates \
        --arguments_file ${params.sample_ID}_mapped_bams.list.txt \
        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
        -MAX_FILE_HANDLES 1000 \
        --REMOVE_SEQUENCING_DUPLICATES false \
        --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES \
        --QUIET \
        --VERBOSITY INFO \
        -M ${params.sample_ID}/dedup_bam/${params.sample_ID}_duplicate_metrics \
        -O ${params.sample_ID}/dedup_bam/${params.sample_ID}_merged_dedup.bam
    """
}


process sort_and_tag{

    publishDir "."
    
    input:
    file(dedup_bam)

    output:
    path("${params.sample_ID}/bam/${params.sample_ID}_merged_dedup_sorted_tagged.bam")
    
    script:
    """
    mkdir -p ${params.sample_ID}/bam

    samtools sort ${dedup_bam} -o  ${params.sample_ID}/bam/${params.sample_ID}_dedup_sorted.bam 
    picard SetNmMdAndUqTags  R=${params.reference_genome} I=${params.sample_ID}/bam/${params.sample_ID}_dedup_sorted.bam \
        O=${params.sample_ID}/bam/${params.sample_ID}_merged_dedup_sorted_tagged.bam  
    
    rm ${params.sample_ID}/bam/${params.sample_ID}_dedup_sorted.bam 
    """

}

process alignment_metrics{

    publishDir "."
    
    input:
    file(tagged_bam)

    output:
    path("${params.sample_ID}/metrics/${params.sample_ID}_merged_AlignmentMetrics.txt")

    script:
    """
    mkdir -p ${params.sample_ID}/metrics/
    picard CollectAlignmentSummaryMetrics R=${params.reference_genome} \
        I=${tagged_bam} \
        O=${params.sample_ID}/metrics/${params.sample_ID}_merged_AlignmentMetrics.txt
    """
    

}


process base_recalibrator{

    publishDir "."

    input:
    file(tagged_bam)

    output:
    //path("${params.sample_ID}/bam/${params.sample_ID}_recal_data.table")
    path("${params.sample_ID}/bam/${params.sample_ID}_bqsr.bam")
    
    script:
    """
    # create BAI file
    samtools index ${tagged_bam} 

    # ip = interval padding (0 should be default anyway)
    
    mkdir -p ${params.sample_ID}/bam

    gatk --java-options "-Xmx16g -Xms16g" \
         BaseRecalibrator \
         -I ${tagged_bam} \
         -R ${params.reference_genome} \
         --known-sites ${params.recalibrator_known_sites} \
         -L ${params.recalibrator_interval_list} \
         --ip 0 \
         -O ${params.sample_ID}/bam/${params.sample_ID}_recal_data.table

    gatk --java-options "-Xmx16g -Xms16g" \
        ApplyBQSR \
        -R ${params.reference_genome} \
	-I ${tagged_bam} \
        --bqsr-recal-file ${params.sample_ID}/bam/${params.sample_ID}_recal_data.table \
	-L ${params.recalibrator_interval_list} \
	--ip 0 -O ${params.sample_ID}/bam/${params.sample_ID}_bqsr.bam
 
    rm ${tagged_bam}
    
    """
}

process haplotype_caller_vcf {
    
    publishDir "."

    input:
    file(bqsr_bam)

    output:
    //path("${params.sample_ID}/bam/${params.sample_ID}_recal_data.table")
    path("${params.sample_ID}/vcf/${params.sample_ID}.vcf.gz")
    
    script:
    """
    mkdir -p ${params.sample_ID}/vcf
    gatk --java-options "-Xmx16g -Xms16g" \
	HaplotypeCaller \
        -R ${params.reference_genome} \
        -I ${bqsr_bam} --native-pair-hmm-threads 4 \
        -L ${params.recalibrator_interval_list} --ip 0 \
        -O ${params.sample_ID}/vcf/${params.sample_ID}.vcf.gz
    """

}
process haplotype_caller_gvcf {
    
    publishDir "."

    input:
    file(bqsr_bam)

    output:
    path("${params.sample_ID}/vcf/${params.sample_ID}.g.vcf.gz")
    
    script:
    """
    mkdir -p ${params.sample_ID}/vcf
    gatk --java-options "-Xmx16g -Xms16g" \
	HaplotypeCaller \
        -R ${params.reference_genome} \
        -I ${bqsr_bam} --native-pair-hmm-threads 4 \
        -L ${params.recalibrator_interval_list} --ip 0 \
        -O ${params.sample_ID}/gvcf/${params.sample_ID}.g.vcf.gz \
        -ERC GVCF
    """

}


//   gatk --java-options "-Xmx16g -Xms16g" \
 //       ApplyBQSR \
//        -R ${params.reference_genome} \
//	-I ${tagged_bam} \
//        --bqsr-recal-file ${params.sample_ID}_recal_data.table \
//	-L ${params.recalibrator_interval_list} \
//	--ip 0 -O ${params.sample_ID}/bam/${params.sample_ID}_bqsr.bam


