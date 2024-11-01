localrules: filterPass1Junctions

# rule star_index:
#     input:
#         fasta="resources/star_genome_input/{genome}/genome.fasta",
#         annotation="resources/star_genome_input/{genome}/genome.gtf",
#     output:
#         directory("resources/star_genome/{genome}/")
#     params:
#         extra = lambda wildcards, input: f"--sjdbGTFfile {input.annotation}"
#     cache: True
#     threads: 4
#     log:
#         "logs/{genome}/star_index_genome.log",
#     wrapper:
#         "0.84.0/bio/star/index"

   
rule fastp_pe:
    input:
        sample=[PREFIX + "/fastq/{sample}_0.fastq.gz",
            PREFIX + "/fastq/{sample}_1.fastq.gz"],
    output:
        trimmed=[
            PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_0.fastq.gz",
            PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_1.fastq.gz"],
        unpaired1=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.u1.fastq",
        unpaired2=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.u2.fastq",
        failed=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.failed.fastq",
        html=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.html",
        json=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.json"
    log:
        PREFIX + "/fastq_trimmed/logs/fastp/pe/{sample}.log"
    params:
        extra='--detect_adapter_for_pe'
    threads: 16
    resources:
        mem_mb=5000
    wrapper:
        "0.84.0/bio/fastp"


rule align_pass1_pe:
    input:
        fq1=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_0.fastq.gz",
        fq2=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_1.fastq.gz",
        index=config["genomes_dir"] + "/{assembly}/star_genome/"
    output:
        bam=PREFIX + "/{assembly}/bam/pe/pass1/{sample}/Aligned.sortedByCoord.out.bam",
        sj=PREFIX + "/{assembly}/bam/pe/pass1/{sample}/SJ.out.tab",
    log:
        PREFIX + "/{assembly}/bam/pe/pass1/{sample}/Log.txt",
    params:
        index=lambda wc, input: input.index,
        extra="--outSAMtype BAM SortedByCoordinate"
    threads: 16
    resources:
        mem_mb=45000
    wrapper:
        "0.84.0/bio/star/align"  

rule fastp_se:
    input:
        sample=[PREFIX + "/fastq/{sample}.fastq.gz"]
    output:
        trimmed=PREFIX + "/fastq_trimmed/trimmed/se/{sample}_0.fastq.gz",
        failed=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.failed.fastq",
        html=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.html",
        json=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.json"
    log:
        PREFIX + "/fastq_trimmed/logs/fastp/se/{sample}.log"
    threads: 16
    resources:
        mem_mb=5000
    wrapper:
        "0.84.0/bio/fastp"

rule align_pass1_se:
    input:
        fq1=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.fastq.gz",
        index=config["genomes_dir"] + "/{assembly}/star_genome/"
    output:
        bam=PREFIX + "/{assembly}/bam/se/pass1/{sample}/Aligned.sortedByCoord.out.bam",
        sj=PREFIX + "/{assembly}/bam/se/pass1/{sample}/SJ.out.tab",
    log:
        PREFIX + "/{assembly}/bam/se/pass1/{sample}/Log.txt",
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
        sj=PREFIX + "/{assembly}/bam/{layout}/pass1/{sample}/SJ.out.tab"
    output:
        PREFIX + "/{assembly}/bam/{layout}/pass1/{sample}/SJ.out.filtered.tab"
    conda: "../envs/postprocess.yaml"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule align_pass2:
    input:
        fq1 = PREFIX + "/fastq_trimmed/trimmed/{layout}/{file_id}.fastq.gz",
        sjdb = get_pass2_sj,
        index=config["genomes_dir"] + "/{assembly}/star_genome/"
    output:
        bam = PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}/Aligned.sortedByCoord.out.bam",
        chim_junc = PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}/Chimeric.out.junction",
    log:
        PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}//Log.txt"
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input:"--outSAMtype BAM SortedByCoordinate --chimOutType Junctions " + config['params']['star2pass'] 
        + f" --sjdbFileChrStartEnd {input['sjdb']}"
    threads: 16
    resources:
        mem_mb=45000
    wrapper:
        "0.84.0/bio/star/align"   
