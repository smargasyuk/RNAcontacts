localrules: filterPass1Junctions

rule star_index:
    input:
        fasta="resources/star_genome_input/{genome}/genome.fasta",
        annotation="resources/star_genome_input/{genome}/genome.gtf",
    output:
        directory("resources/star_genome/{genome}/")
    params:
        extra = lambda wildcards, input: f"--sjdbGTFfile {input.annotation}"
    cache: True
    threads: 4
    log:
        "logs/{genome}/star_index_genome.log",
    wrapper:
        "0.84.0/bio/star/index"

   
rule fastp_pe:
    input:
        sample=get_pass1_fq_for_adapterremoval,
    output:
        trimmed=["results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}_0.fastq.gz",
         "results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}_1.fastq.gz"],
        unpaired1="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}.u1.fastq",
        unpaired2="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}.u2.fastq",
        failed="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}.failed.fastq",
        html="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}.html",
        json="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}.json"
    log:
        "results/{genome}/{project}/fastq_trimmed/logs/fastp/pe/{sample}.log"
    params:
        extra='--detect_adapter_for_pe'
    threads: 16
    resources:
        mem_mb=5000
    wrapper:
        "0.84.0/bio/fastp"


rule align_pass1:
    input:
        fq1="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}_0.fastq.gz",
        fq2="results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}_1.fastq.gz",
        index="resources/star_genome/{genome}",
    output:
        bam="results/{genome}/{project}/bam/pass1/{sample}/Aligned.sortedByCoord.out.bam",
        sj="results/{genome}/{project}/bam/pass1/{sample}/SJ.out.tab",
    log:
        "logs/{genome}/{project}/bam/pass1/{sample}/Log.txt",
    params:
        index=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate"
    threads: 16
    resources:
        mem_mb=45000
    wrapper:
        "0.84.0/bio/star/align"  

rule filterPass1Junctions:
    input:
        "results/{genome}/{project}/bam/pass1/{sample}/SJ.out.tab"
    output:
        "results/{genome}/{project}/bam/pass1/{sample}/SJ.out.filtered.tab"
    conda: "../envs/postprocess.yaml"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule align_pass2:
    input:
        fq1 = "results/{genome}/{project}/fastq_trimmed/trimmed/pe/{sample}_{mate}.fastq.gz",
        sjdb = get_pass2_sj,
        index="resources/star_genome/{genome}"
    output:
        bam = "results/{genome}/{project}/bam/pass2/{sample}_{mate}/Aligned.sortedByCoord.out.bam",
        chim_junc = "results/{genome}/{project}/bam/pass2/{sample}_{mate}/Chimeric.out.junction",
    log:
        "logs/{genome}/{project}/bam/pass2/{sample}_{mate}/Log.txt"
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input:"--outSAMtype BAM SortedByCoordinate --chimOutType Junctions " + config['params']['star2pass'] 
        + f" --sjdbFileChrStartEnd {input['sjdb']}"
    threads: 16
    resources:
        mem_mb=45000
    wrapper:
        "0.84.0/bio/star/align"   
