# configfile: "config.yaml"

localrules: filterPass1Junctions, mergeKnownJunctions, extractChimericJunctions, extractNeoJunctions, clusterDonorsAcceptors, clusterJunctions, extractContacts

rule alignControlPE:
    input:
        fastqs = lambda wildcards: config['samples'][wildcards["id"]],
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"data/{wildcards['genome']}/bam/pass1/{wildcards['id']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate',
	runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        sj = "data/{genome}/bam/pass1/{id}/SJ.out.tab"
    conda: "envs/env.yaml"
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --readFilesIn {input.fastqs} --outFileNamePrefix {params.out_prefix}
"""

rule filterPass1Junctions:
    input:
        "data/{genome}/bam/pass1/{id}/SJ.out.tab"
    output:
        "data/{genome}/bam/pass1/{id}/SJ.out.filtered.tab"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule alignRICSE:
    input:
        fastq = lambda wildcards: f'{config["samples"][wildcards["id"]][int(wildcards["mate"])]}',
        junctions = lambda wildcards: expand("data/{genome}/bam/pass1/{id}/SJ.out.filtered.tab", id=config['samples_control'], genome=wildcards["genome"]),
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"data/{wildcards['genome']}/bam/pass2/{wildcards['id']}_{wildcards['mate']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate --chimOutType Junctions --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreJunctionNonGTAG  -1 --scoreGapNoncan -1 --scoreGapATAC -1 --scoreGapGCAG -1 --chimSegmentReadGapMax 3 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5',
        runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        bam = "data/{genome}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        chm = "data/{genome}/bam/pass2/{id}_{mate}/Chimeric.out.junction"
    conda: "envs/env.yaml"
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --sjdbFileChrStartEnd {input.junctions} --readFilesIn {input.fastq} --outFileNamePrefix {params.out_prefix}
"""

rule mergeKnownJunctions:
    input:
        lambda wildcards: expand("data/{genome}/bam/pass1/{id}/SJ.out.tab", genome = config['genome'], id=config['samples_control']),
	lambda wildcards: f"{config['star']['genome_path']}sjdbList.out.tab"
    output:
        "data/{genome}/all_sj.tsv"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input} | sort -u > {output}
"""

rule extractChimericJunctions:
    input:
        mate0 = "data/{genome}/bam/pass2/{id}_0/Chimeric.out.junction",
        mate1 = "data/{genome}/bam/pass2/{id}_1/Chimeric.out.junction"
    output:
        "data/{genome}/junctions/{id}/Chimeric.tsv"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9 <(cut -f1-6,10 {input.mate0} | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') \
           <(cut -f1-6,10 {input.mate1} | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | \
           grep -v GL > {output}
"""

rule extractNeoJunctions:
    input:
        mate0 = "data/{genome}/bam/pass2/{id}_0/Aligned.sortedByCoord.out.bam",
	mate1 = "data/{genome}/bam/pass2/{id}_1/Aligned.sortedByCoord.out.bam",
        junctions = "data/{genome}/all_sj.tsv"
    output:
        "data/{genome}/junctions/{id}/Neo.tsv"
    conda: "envs/env.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9 <(samtools view {input.mate0} | perl scripts/neo.pl {input.junctions} 0) \
           <(samtools view {input.mate1} | perl scripts/neo.pl {input.junctions} 1) | grep -v GL > {output}        
"""

rule clusterDonorsAcceptors:
    input:
        lambda wildcards: expand("data/{genome}/junctions/{id}/{jtype}.tsv",
         id=config['samples'].keys(), jtype=["Neo", "Chimeric"], genome = config['genome'])
    output:
        donors = "data/{genome}/junctions/donors.bed",
        acceptors = "data/{genome}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    conda: "envs/env.yaml"
    shell:
        """
cut -f1-4 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.donors}
cut -f5-8 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.acceptors}
"""

rule clusterJunctions:
    input:
        junctions = "data/{genome}/junctions/{id}/{jtype}.tsv",
	donors = "data/{genome}/junctions/donors.bed",
	acceptors = "data/{genome}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "data/{genome}/clusters/{id}/{jtype}.tsv"
    conda: "envs/env.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort-bed {input.junctions} | \
intersectBed -a stdin -b {input.donors} -wa -wb | cut -f5- | sort-bed - | \
intersectBed -a stdin -b {input.acceptors} -wa -wb | cut -f5- | sort -k1,1 > {output}
"""

rule extractContacts:
    input:
        "data/{genome}/clusters/{id}/{jtype}.tsv",
    output:
        "data/{genome}/contacts/{id}/{jtype}.tsv",
    shell:
        """
mkdir -p $(dirname {output})
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]++}}END{{for(j in n){{print j"\t"n[j]}}}}' {input} | sort > {output}
"""

rule expandStructures:
    input: "data/input/known_structures.bed"
    output: "data/{genome}/tmp/known_structures_{range}.bed"
    shell: """
mkdir -p $(dirname {output})
bedops -u --range {wildcards.range} {input} > {output}
"""

rule readsNearStructures:
    input: 
        bam="data/{genome}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        structures="data/{genome}/tmp/known_structures_{range}.bed"
    output:
        bam = "data/{genome}/views/near_known_structure/{id}_{mate}_{range}.bam",
        bed = "data/{genome}/views/near_known_structure/{id}_{mate}_{range}_split.bed"
    shell: """
mkdir -p $(dirname {output.bam})
samtools view -@4 -h -b -L {input.structures} {input.bam} > {output.bam}
paste <(bedtools bamtobed -bed12 -i {output.bam}) <(bedtools bamtobed -cigar -i {output.bam})\
 | awk -v 'OFS=\\t' '($3 - $2)<{wildcards.range}' | awk -v 'OFS=\\t' '$19~/N/' > {output.bed}
"""

rule readsWNeoJunctions:
    input: 
        bam="data/{genome}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        junctions = "data/{genome}/all_sj.tsv"
    output:
        "data/{genome}/views/reads_w_neo/{id}_{mate}_{range}.bam"
    shell:"""
mkdir -p $(dirname {output})
python scripts/reads_w_new_junctions.py --max-length {wildcards.range} {input.bam} {input.junctions} | samtools view -h -b - > {output}
"""

rule mergePairedNeoReads:
    input:
        bam0 = "data/{genome}/views/reads_w_neo/{id}_0_{range}.bam",
        bam1 = "data/{genome}/views/reads_w_neo/{id}_1_{range}.bam"
    output:
        "data/{genome}/views/reads_w_neo_bed/{id}_{range}.bed"
    shell: """
mkdir -p $(dirname {output})
sort-bed \
<(bedtools bamtobed -bed12 -i {input.bam0}) \
<(bedtools bamtobed -bed12 -i {input.bam1} | awk -v 'OFS=\\t' '$9="0,0,255"') > {output}
"""

rule junctionsToIntronicBed:
    input: config['star']['genome_path'] + "/sjdbList.out.tab"
    output: "data/{genome}/introns.bed"
    shell:"""
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$3=$3+1' {input} | sort-bed - > {output}    
"""

rule readsBed12ExtractJunctions:
    input:
        bed="data/{genome}/views/reads_w_neo_bed/{id}_-1.bed",
        junctions='data/{genome}/introns.bed'
    output:
        all = "data/{genome}/views/neo_junctions_bed/{id}.bed",
        intronic = "data/{genome}/views/neo_junctions_bed/{id}_intronic.bed"
    shell: """
mkdir -p $(dirname {output})
cat {input.bed} | python scripts/bed12_split_junctions.py | python scripts/bed12_filter_junctions.py --sjdb {input.junctions} |\
sort-bed - > {output.all}
bedops -e {output.all} {input.junctions} > {output.intronic}
"""


rule allContacts:
    input:
        expand("data/{genome}/contacts/{id}/{jtype}.tsv", genome = config['genome'], id=config['samples'].keys(), jtype=["Neo", "Chimeric"]) 


rule getRICAligned:
    input:
        expand("data/{genome}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam", genome = config['genome'], id=config['samples'].keys(), mate=[0,1])


rule allReadsNearStructures:
    input: 
        expand("data/{genome}/views/near_known_structure/{id}_{mate}_{range}.bam", genome = config['genome'], id=config['samples'].keys(), mate=[0,1], range=[1000, 10000, 30000])


rule allReadsWNeo:
    input: 
        expand("data/{genome}/views/reads_w_neo_bed/{id}_{range}.bed", genome = config['genome'], id=config['samples'].keys(), range=[1000, 10000, 30000])


rule AllNeoJunctions:
    input: expand("data/{genome}/views/neo_junctions_bed/{id}_intronic.bed", genome = config['genome'], id=config['samples'].keys())

