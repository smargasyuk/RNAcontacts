rule mergeKnownJunctions:
    input: get_known_junctions
    output:
        "results/{genome}/{project}/all_sj.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input} | sort -u > {output}
"""

rule extractChimericJunctions:
    input:
        mate0 = "results/{genome}/{project}/bam/pass2/{id}_0/Chimeric.out.junction",
        mate1 = "results/{genome}/{project}/bam/pass2/{id}_1/Chimeric.out.junction"
    output:
        "results/{genome}/{project}/junctions/{id}/Chimeric.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort --parallel=8 -S4G -k9,9 <(cut -f1-6,10 {input.mate0} | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') \
           <(cut -f1-6,10 {input.mate1} | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | \
           grep -v GL > {output}
"""

rule extractNeoJunctions:
    input:
        mate0 = "results/{genome}/{project}/bam/pass2/{id}_0/Aligned.sortedByCoord.out.bam",
	mate1 = "results/{genome}/{project}/bam/pass2/{id}_1/Aligned.sortedByCoord.out.bam",
        junctions = "results/{genome}/{project}/all_sj.tsv"
    output:
        "results/{genome}/{project}/junctions/{id}/Neo.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9  -S4G --parallel=8 <(samtools view {input.mate0} | perl workflow/scripts/neo.pl {input.junctions} 0) \
           <(samtools view {input.mate1} | perl workflow/scripts/neo.pl {input.junctions} 1) | grep -v GL |\
            awk -v 'OFS=\\t' '{{print $1,$2,$2+1,$3,$4,$5,$5+1,$6,$7,$8}}' > {output}       
"""

rule clusterDonorsAcceptors:
    input: get_all_junction_files
    output:
        donors = "results/{genome}/{project}/junctions/donors.bed",
        acceptors = "results/{genome}/{project}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    conda: "../envs/postprocess.yaml"
    shell:
        """
cut -f1-4 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.donors}
cut -f5-8 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.acceptors}
"""

rule clusterJunctions:
    input:
        junctions = "results/{genome}/{project}/junctions/{id}/{jtype}.tsv",
	donors = "results/{genome}/{project}/junctions/donors.bed",
	acceptors = "results/{genome}/{project}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "results/{genome}/{project}/clusters/{id}/{jtype}.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort-bed {input.junctions} | \
intersectBed -a stdin -b {input.donors} -wa -wb | cut -f5- | sort-bed - | \
intersectBed -a stdin -b {input.acceptors} -wa -wb | cut -f5- | sort -k1,1 > {output}
"""

rule extractContacts:
    input:
        "results/{genome}/{project}/clusters/{id}/{jtype}.tsv",
    output:
        "results/{genome}/{project}/contacts/{id}/{jtype}.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]++}}END{{for(j in n){{print j"\t"n[j]}}}}' {input} | sort > {output}
"""