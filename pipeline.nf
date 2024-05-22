
params.sra_file = "SRR_Acc_List.txt"
params.outdir = "/mnt/c/Users/javi_/Escritorio/entrenamiento/resultados"
params.reference = "/mnt/c/Users/javi_/Escritorio/entrenamiento/GCF_000005845.2_ASM584v2_genomic.fna"
params.annotation = "/mnt/c/Users/javi_/Escritorio/entrenamiento/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.57.gff3"
params.script_R = ""

// SRA IDs
process get_ids{

    input:
        path sra_file
    output:
        stdout
    script:
        """
        ids=\$(cat ${sra_file})
        echo -n \$ids
        """
}


// Download SRA files
process descargaSRA {

    publishDir("sra_files", mode: 'copy')

    input:
        each id_sample

    output:
        path "${id_sample}.sra"

    script:
        """
        prefetch  ${id_sample} -O "."
        mv ${id_sample}/*.sra . ; rm -r ${id_sample}
        """
}


// Convert SRA to FASTQ
process convert_to_fastq {

    publishDir("fastq_files",mode:'copy')

    input:
        each sra_file
    output:
        path "SRR*_{1,2}.fastq"
    script:
        """
        fastq-dump --split-files ${sra_file}
        """
}


// Filters out low quality reads from FASTQ
process filter_reads {

    publishDir("fastp_files", mode: 'copy')

    input:
        each fastq_files
    output:
        path "clean_SRR*_{1,2}.fastq"

    script:
        """
        id=\$(basename ${fastq_files[0]} | awk -F "/" '{print \$NF}' | awk -F "_" '{print \$1}')
        fastp -i ${fastq_files[0]} -I ${fastq_files[1]} -o clean_\${id}_1.fastq -O clean_\${id}_2.fastq --cut_tail --cut_window_size 5 --cut_mean_quality 30
        """
}


// Mapping reads to reference to generate BAM files

process mapping_to_reference{
    publishDir("bam_files",mode:"copy")

    input:
        path reference
        each fastp_files

    output:
        path "*sorted.bam"

    script:
        """
        id=\$(basename ${fastp_files[0]} | awk -F "_" '{print \$2}')
        bwa index ${reference} -p genome
        bwa mem genome ${fastp_files[0]} ${fastp_files[1]} > \${id}.sam
        samtools view \${id}.sam -q 50 -o \${id}.bam
        samtools sort \${id}.bam -o \${id}_sorted.bam
        samtools index \${id}_sorted.bam
        """

}

// Get the Count Matrix and DEG analysis

process perform_DEG_analysis{
    input:
        path script_R
        path annotation
        path fastp_files
    
    output:
        stdout
    
    script:
        """
        Rscript ${script_R} ${fastp_files} ${annotation}
        """



}
// Main workflow
workflow {
    // Input Channels
    ref_ch = Channel.fromPath(params.reference)
    ids_ch = Channel.fromPath(params.sra_file).splitText().map{it -> it.trim()}.collect()
    script_R_ch = Channel.fromPath(params.script_R)
    annotation_ch = Channel.fromPath(params.annotation)
    // Start processing
    sra_ch = descargaSRA(ids_ch)
    fastq_ch = convert_to_fastq(sra_ch)
    fastp_ch = filter_reads(fastq_ch)
    map_ch = mapping_to_reference(ref_ch,fastp_ch)
    perform_DEG_analysis(script_R_ch,fastp_ch,annotation_ch)

}

