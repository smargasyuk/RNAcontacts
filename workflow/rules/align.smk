rule align_pass1:
    input:
        unpack(get_pass1_fq),
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
        fq1 = get_pass2_fq,
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

