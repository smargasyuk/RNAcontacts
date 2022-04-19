# configfile: "config.yaml"

localrules: filterPass1Junctions, mergeKnownJunctions, extractChimericJunctions, extractNeoJunctions, clusterDonorsAcceptors, clusterJunctions, extractContacts
wildcard_constraints:
    jtype="Neo|Chimeric"

import subprocess
from pathlib import Path
import distinctipy
PROJECT_PREFIX = f"{config['genome']}/{config['project']}"

rule alignControlPE:
    input:
        fastqs = lambda wildcards: config['samples'][wildcards["id"]],
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"{wildcards['project_prefix']}/bam/pass1/{wildcards['id']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat',
	runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        sj = "{project_prefix}/bam/pass1/{id}/SJ.out.tab"
    conda: "envs/env.yaml"
    threads: 16
    resources:
        mem_mb=45000
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --readFilesIn {input.fastqs} --outFileNamePrefix {params.out_prefix}
"""

rule filterPass1Junctions:
    input:
        "{project_prefix}/bam/pass1/{id}/SJ.out.tab"
    output:
        "{project_prefix}/bam/pass1/{id}/SJ.out.filtered.tab"
    shell:
        """
awk -v 'OFS="\t"' '$5 == 1' {input} > {output}
"""

rule alignRICSE:
    input:
        fastq = lambda wildcards: f'{config["samples"][wildcards["id"]][int(wildcards["mate"])]}',
        junctions = lambda wildcards: expand("{project_prefix}/bam/pass1/{id}/SJ.out.filtered.tab", project_prefix = wildcards['project_prefix'], id=config['samples_control']),
        genome = lambda wildcards: config['star']['genome_path']
    params:
        STAR = lambda wildcards: config['star']['program'],
        out_prefix = lambda wildcards: f"{wildcards['project_prefix']}/bam/pass2/{wildcards['id']}_{wildcards['mate']}/",
        STAR_params = '--runMode alignReads --outSAMtype BAM SortedByCoordinate --chimOutType Junctions --chimSegmentMin 15 --chimJunctionOverhangMin 15 --chimScoreJunctionNonGTAG  -1 --scoreGapNoncan -1 --scoreGapATAC -1 --scoreGapGCAG -1 --chimSegmentReadGapMax 3 --outFilterMatchNminOverLread 0.5 --outFilterScoreMinOverLread 0.5 --readFilesCommand zcat',
        runThreadN = lambda wildcards: config['star']['runThreadN']
    output:
        bam = "{project_prefix}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        chm = "{project_prefix}/bam/pass2/{id}_{mate}/Chimeric.out.junction"
    conda: "envs/env.yaml"
    threads: 16
    resources:
        mem_mb=45000
    shell: 
        """
mkdir -p $(dirname {output})
{params.STAR} {params.STAR_params} --runThreadN {params.runThreadN} --genomeDir {input.genome} --sjdbFileChrStartEnd {input.junctions} --readFilesIn {input.fastq} --outFileNamePrefix {params.out_prefix}
"""

rule mergeKnownJunctions:
    input:
        lambda wildcards: expand("{project_prefix}/bam/pass1/{id}/SJ.out.tab", project_prefix = wildcards['project_prefix'], id=config['samples_control']),
	lambda wildcards: f"{config['star']['genome_path']}sjdbList.out.tab"
    output:
        "{project_prefix}/all_sj.tsv"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input} | sort -u > {output}
"""

rule extractChimericJunctions:
    input:
        mate0 = "{project_prefix}/bam/pass2/{id}_0/Chimeric.out.junction",
        mate1 = "{project_prefix}/bam/pass2/{id}_1/Chimeric.out.junction"
    output:
        "{project_prefix}/junctions/{id}/Chimeric.tsv"
    shell:
        """
mkdir -p $(dirname {output})
sort --parallel=8 -S4G -k9,9 <(cut -f1-6,10 {input.mate0} | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') \
           <(cut -f1-6,10 {input.mate1} | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | \
           grep -v GL > {output}
"""

rule extractNeoJunctions:
    input:
        mate0 = "{project_prefix}/bam/pass2/{id}_0/Aligned.sortedByCoord.out.bam",
	mate1 = "{project_prefix}/bam/pass2/{id}_1/Aligned.sortedByCoord.out.bam",
        junctions = "{project_prefix}/all_sj.tsv"
    output:
        "{project_prefix}/junctions/{id}/Neo.tsv"
    conda: "envs/env.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9  -S4G --parallel=8 <(samtools view {input.mate0} | perl scripts/neo.pl {input.junctions} 0) \
           <(samtools view {input.mate1} | perl scripts/neo.pl {input.junctions} 1) | grep -v GL |\
            awk -v 'OFS=\\t' '{{print $1,$2,$2+1,$3,$4,$5,$5+1,$6,$7,$8}}' > {output}       
"""

rule clusterDonorsAcceptors:
    input:
        lambda wildcards: expand("{project_prefix}/junctions/{id}/{jtype}.tsv",
         id=config['samples'].keys(), project_prefix=wildcards['project_prefix'], jtype=["Neo", "Chimeric"], genome = config['genome'])
    output:
        donors = "{project_prefix}/junctions/donors.bed",
        acceptors = "{project_prefix}/junctions/acceptors.bed"
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
        junctions = "{project_prefix}/junctions/{id}/{jtype}.tsv",
	donors = "{project_prefix}/junctions/donors.bed",
	acceptors = "{project_prefix}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "{project_prefix}/clusters/{id}/{jtype}.tsv"
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
        "{project_prefix}/clusters/{id}/{jtype}.tsv",
    output:
        "{project_prefix}/contacts/{id}/{jtype}.tsv",
    shell:
        """
mkdir -p $(dirname {output})
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10]++}}END{{for(j in n){{print j"\t"n[j]}}}}' {input} | sort > {output}
"""

rule expandStructures:
    input: "data/input/known_structures.bed"
    output: "{project_prefix}/tmp/known_structures_{range}.bed"
    shell: """
mkdir -p $(dirname {output})
bedops -u --range {wildcards.range} {input} > {output}
"""

rule readsNearStructures:
    input: 
        bam="{project_prefix}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        structures="{project_prefix}/tmp/known_structures_{range}.bed"
    output:
        bam = "{project_prefix}/views/near_known_structure/{id}_{mate}_{range}.bam",
        bed = "{project_prefix}/views/near_known_structure/{id}_{mate}_{range}_split.bed"
    shell: """
mkdir -p $(dirname {output.bam})
samtools view -@4 -h -b -L {input.structures} {input.bam} > {output.bam}
paste <(bedtools bamtobed -bed12 -i {output.bam}) <(bedtools bamtobed -cigar -i {output.bam})\
 | awk -v 'OFS=\\t' '($3 - $2)<{wildcards.range}' | awk -v 'OFS=\\t' '$19~/N/' > {output.bed}
"""

rule readsWNeoJunctions:
    input: 
        bam="{project_prefix}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam",
        junctions = "{project_prefix}/all_sj.tsv"
    output:
        "{project_prefix}/views/reads_w_neo/{id}_{mate}_{range}.bam"
    shell:"""
mkdir -p $(dirname {output})
python scripts/reads_w_new_junctions.py --max-length {wildcards.range} {input.bam} {input.junctions} | samtools view -h -b - > {output}
"""

rule mergePairedNeoReads:
    input:
        bam0 = "{project_prefix}/views/reads_w_neo/{id}_0_{range}.bam",
        bam1 = "{project_prefix}/views/reads_w_neo/{id}_1_{range}.bam"
    output:
        "{project_prefix}/views/reads_w_neo_bed/{id}_{range}.bed"
    shell: """
mkdir -p $(dirname {output})
sort-bed \
<(bedtools bamtobed -bed12 -i {input.bam0}) \
<(bedtools bamtobed -bed12 -i {input.bam1} | awk -v 'OFS=\\t' '$9="0,0,255"') > {output}
"""

rule junctionsToIntronicBed:
    input: config['star']['genome_path'] + "/sjdbList.out.tab"
    output: "{project_prefix}/introns.bed"
    shell:"""
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$3=$3+1' {input} | sort-bed - > {output}    
"""

rule readsBed12ExtractJunctions:
    input:
        bed="{project_prefix}/views/reads_w_neo_bed/{id}_-1.bed",
        junctions='{project_prefix}/introns.bed'
    output:
        all = "{project_prefix}/views/neo_junctions_bed/{id}.bed",
        intronic = "{project_prefix}/views/neo_junctions_bed/{id}_intronic.bed"
    shell: """
mkdir -p $(dirname {output.all})
cat {input.bed} | python scripts/bed12_split_junctions.py | python scripts/bed12_filter_junctions.py --sjdb {input.junctions} |\
sort-bed - > {output.all}
bedops -e {output.all} {input.junctions} > {output.intronic}
"""

rule junctionsToBedJ:
    input:
        tsv = "{project_prefix}/junctions/{id}/{jtype}.tsv"
    output:
        bed = "{project_prefix}/contacts_v2/{id}/{jtype}.bed"
    shell: """
mkdir -p $(dirname {output.bed})
cat {input.tsv} | awk -v 'OFS=\\t' '$1==$5' | \
cut -f1,2,4,7,9,10 |  \
awk -v 'OFS=\\t' -v 'name_prefix={wildcards.id}_{wildcards.jtype}' '{{print $1,$2,$4,name_prefix"_"NR,1,$3,$5,$6}}' | \
awk -v 'OFS=\\t' '{{if ($2 > $3){{t=$3;$3=$2;$2=t}}; print}}' |\
awk -v 'OFS=\\t' '$2<$3' |\
sort-bed - > {output.bed}
"""

rule filterContactsBothNeo:
    input: 
        contacts="{project_prefix}/contacts_v2/{id}/{jtype_long}.bed",
        junctions = "{project_prefix}/introns.bed"
    output: "{project_prefix}/contacts_v2/{id}/{jtype_long}_neoR.bed"
    shell: """
mkdir -p $(dirname {output})
cat {input.contacts} | python scripts/bedj_filter_junctions.py --sjdb {input.junctions} > {output}
"""



rule filterContactsByLength:
    input: "{project_prefix}/contacts_v2/{id}/{jtype_long}.bed",
    output: "{project_prefix}/contacts_v2/{id}/{jtype_long}_len{range,\d+}.bed"
    shell: """
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$3-$2<{wildcards.range}'  {input} > {output}
"""

rule filterContactsWIntrons:
    input: 
        contacts ="{project_prefix}/contacts_v2/{id}/{jtype_long}.bed",
        junctions = "{project_prefix}/introns.bed"
    output: "{project_prefix}/contacts_v2/{id}/{jtype_long}_inIntrons.bed"
    shell: """
mkdir -p $(dirname {output})    
bedops -e {input.contacts} {input.junctions} > {output}
"""

rule mergeChimericNeo:
    input:
        neo = "{project_prefix}/contacts_v2/{id}/Neo_{other}.bed",
        chim = "{project_prefix}/contacts_v2/{id}/Chimeric_{other}.bed"
    output: "{project_prefix}/contacts_v2/{id}/All_{other}_merged.bed"
    shell: """
bedops -u {input.neo} {input.chim} > {output}
"""

rule mergeContacts:
    input: "{project_prefix}/contacts_v2/{id}/All_{other}.bed"
    output: "{project_prefix}/contacts_v2/{id}/All_{other}_aggregated.bed"
    shell: """
cut -f 1,2,3 {input} | sort | uniq -c | sed -r 's/([0-9]) /\\1\\t/' |\
awk -v 'OFS=\\t' -v 'name_prefix=id' '{{print $2,$3,$4,name_prefix"_"NR,$1,"+"}}' |\
sort-bed - > {output}
"""

rule mergeProjectContacts:
    input: lambda wildcards: expand("{project_prefix}/contacts_v2/{id}/All_len40000_inIntrons_neoR_merged_aggregated.bed", project_prefix=wildcards['project_prefix'], id=config["samples_RIC"])
    output: "{project_prefix}/contacts_v2/AllContacts.bed"
    shell: """
sort-bed {input} > {output}
"""

rule prettyShowContacts:
    input: "{project_prefix}/contacts_v2/{id}/{jtype_longer}.bed"
    output: "{project_prefix}/contacts_v2/{id}/{jtype_longer}_view.bed"
    shell: """
cat {input} | python scripts/bedj_generate_view.py > {output}
"""


rule mergeReplicatesForView:
    input: 
        bed_files = lambda wildcards: expand("data/{project_prefix}/contacts_v2/{id}/All_len40000_inIntrons_neoR_merged_aggregated_view.bed", project_prefix=wildcards['project_prefix'], id=config['samples_RIC'])
    output: "data/RIC-hub/input/{project_prefix}.bed"
    run:
        Path(str(output)).parent.mkdir(parents=True, exist_ok=True)
        colors = distinctipy.get_colors(len(input.bed_files), pastel_factor=0.7)
        color_strings = [','.join([str(int(i * 255)) for i in rgb_frac_tuple]) for rgb_frac_tuple in colors]
        for idx, file in enumerate(input.bed_files):
            subprocess.run(f"awk -v 'OFS=\\t' '$9=\"{color_strings[idx]}\"' {file} >> {output}", shell=True)


rule allContacts:
    input:
        expand("data/{project_prefix}/contacts/{id}/{jtype}.tsv", project_prefix=PROJECT_PREFIX, id=config['samples'].keys(), jtype=["Neo", "Chimeric"]) 


rule getRICAligned:
    input:
        expand("data/{project_prefix}/bam/pass2/{id}_{mate}/Aligned.sortedByCoord.out.bam", project_prefix=PROJECT_PREFIX, id=config['samples'].keys(), mate=[0,1])


rule allReadsNearStructures:
    input: 
        expand("data/{project_prefix}/views/near_known_structure/{id}_{mate}_{range}.bam", project_prefix=PROJECT_PREFIX, id=config['samples'].keys(), mate=[0,1], range=[1000, 10000, 30000])


rule allReadsWNeo:
    input: 
        expand("data/{project_prefix}/views/reads_w_neo_bed/{id}_{range}.bed", project_prefix=PROJECT_PREFIX, id=config['samples'].keys(), range=[1000, 10000, 30000])


rule AllNeoJunctions:
    input: expand("data/{project_prefix}/views/neo_junctions_bed/{id}_intronic.bed", project_prefix=PROJECT_PREFIX, id=config['samples'].keys())

rule contactsV2:
    input: expand("data/{project_prefix}/contacts_v2/{id}/All_len40000_inIntrons_neoR_merged_aggregated_view.bed", project_prefix=PROJECT_PREFIX, id=config['samples'].keys()) + [f"data/{PROJECT_PREFIX}/contacts_v2/AllContacts.bed"]

rule viewForHub:
    input: f"data/RIC-hub/input/{PROJECT_PREFIX}.bed"

