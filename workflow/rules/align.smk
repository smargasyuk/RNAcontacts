localrules: filterPass1Junctions

rule star_index:
    input:
        fasta=config["genomes_dir"] + "/{assembly}/{assembly}.fa",
        annotation=config["genomes_dir"] + "/{assembly}/{assembly}.annotation.gtf",
    output:
        directory(config["genomes_dir"] + "/{assembly}/star_genome/")
    params:
        extra = lambda wildcards, input: f"--sjdbGTFfile {input.annotation}"
    cache: True
    threads: 4
    resources:
        mem_mb=10000,
        runtime=600
    log:
        "logs/{assembly}/star_index_genome.log",
    wrapper:
        "0.84.0/bio/star/index"

   
rule fastp_pe:
    input:
        sample=[PREFIX + "/fastq/pe/{sample}/{sample}_1.fastq.gz",
            PREFIX + "/fastq/pe/{sample}/{sample}_2.fastq.gz"],
    output:
        trimmed=[
            PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_1.fastq.gz",
            PREFIX + "/fastq_trimmed/trimmed/pe/{sample}_2.fastq.gz"],
        unpaired1=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.u1.fastq",
        unpaired2=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.u2.fastq",
        failed=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.failed.fastq",
        html=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.html",
        json=PREFIX + "/fastq_trimmed/trimmed/pe/{sample}.json"
    log:
        PREFIX + "/fastq_trimmed/logs/fastp/pe/{sample}.log"
    params:
        extra='--detect_adapter_for_pe'
    threads: 4
    resources:
        mem_mb=5000,
        runtime=600
    wrapper:
        "0.84.0/bio/fastp"

rule fastp_se:
    input:
        sample=[PREFIX + "/fastq/se/{sample}/{sample}.fastq.gz"]
    output:
        trimmed=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.fastq.gz",
        failed=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.failed.fastq",
        html=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.html",
        json=PREFIX + "/fastq_trimmed/trimmed/se/{sample}.json"
    log:
        PREFIX + "/fastq_trimmed/logs/fastp/se/{sample}.log"
    threads: 4
    resources:
        mem_mb=5000,
        runtime=600
    wrapper:
        "0.84.0/bio/fastp"

rule align_pass1:
    input:
        fq1=PREFIX + "/fastq_trimmed/trimmed/{layout}/{file_id}.fastq.gz",
        index=config["genomes_dir"] + "/{assembly}/star_genome/"
    output:
        bam=PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/Aligned.sortedByCoord.out.bam",
        sj=PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/SJ.out.tab",
    log:
        PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/Log.txt",
    params:
        index=lambda wc, input: input.index,
        extra=" --outSAMtype BAM SortedByCoordinate --limitOutSJcollapsed 100000000 --limitIObufferSize=3000000000 "
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 45000 + 45000 * (attempt - 1),
        runtime=600
    wrapper:
        "0.84.0/bio/star/align"  

rule filterPass1Junctions:
    input:
        sj=PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/SJ.out.tab"
    output:
        PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/SJ.out.filtered.tab"
    conda: "../envs/postprocess.yaml"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule align_pass2:
    input:
        fq1 = PREFIX + "/fastq_trimmed/trimmed/{layout}/{file_id}.fastq.gz",
        sjdb = PREFIX + "/{assembly}/bam/{layout}/pass1/{file_id}/SJ.out.filtered.tab",
        index=config["genomes_dir"] + "/{assembly}/star_genome/"
    output:
        bam = PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}/Aligned.sortedByCoord.out.bam",
        chim_junc = PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}/Chimeric.out.junction",
    log:
        PREFIX + "/{assembly}/bam/{layout}/pass2/{file_id}/Log.txt"
    params:
        index=lambda wc, input: input.index,
        extra=lambda wc, input: " --limitOutSJcollapsed 100000000 --limitIObufferSize=3000000000 --limitSjdbInsertNsj 10000000 " + " --outSAMtype BAM SortedByCoordinate --chimOutType Junctions " + config['params']['star2pass'] 
        + f" --sjdbFileChrStartEnd {input['sjdb']}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: 45000 + 45000 * (attempt - 1),
        runtime=600
    wrapper:
        "0.84.0/bio/star/align"   


rule bam_p2_all:
    input: get_pass2_bam
