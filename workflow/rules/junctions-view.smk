conda: "workflow/envs/postprocess.yaml"


rule junctionsToIntronicBed:
    input: "results/{genome}/{project}/all_sj.tsv"
    output: "results/{genome}/{project}/introns.bed"
    conda: "../envs/postprocess.yaml"
    shell:"""
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$3=$3+1' {input} | sort-bed - > {output}    
"""


rule junctionsToBedJ:
    input:
        tsv = "results/{genome}/{project}/junctions/{id}/{jtype}.tsv"
    output:
        bed = "results/{genome}/{project}/junctions-view/{id}/{jtype}.bed"
    conda: "../envs/postprocess.yaml"
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
        contacts="results/{genome}/{project}/junctions-view/{id}/{jtype_long}.bed",
        junctions = "results/{genome}/{project}/introns.bed"
    output: "results/{genome}/{project}/junctions-view/{id}/{jtype_long}_neoR.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
mkdir -p $(dirname {output})
cat {input.contacts} | python workflow/scripts/bedj_filter_junctions.py --sjdb {input.junctions} > {output}
"""


rule filterContactsByLength:
    input: "results/{genome}/{project}/junctions-view/{id}/{jtype_long}.bed",
    output: "results/{genome}/{project}/junctions-view/{id}/{jtype_long}_len{range,\d+}.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$3-$2<{wildcards.range}'  {input} > {output}
"""


rule filterContactsWIntrons:
    input: 
        contacts ="results/{genome}/{project}/junctions-view/{id}/{jtype_long}.bed",
        junctions = "results/{genome}/{project}/introns.bed"
    output: "results/{genome}/{project}/junctions-view/{id}/{jtype_long}_inIntrons.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
mkdir -p $(dirname {output})    
bedops -e {input.contacts} {input.junctions} > {output}
"""


rule mergeChimericNeo:
    input:
        neo = "results/{genome}/{project}/junctions-view/{id}/Neo_{other}.bed",
        chim = "results/{genome}/{project}/junctions-view/{id}/Chimeric_{other}.bed"
    output: "results/{genome}/{project}/junctions-view/{id}/All_{other}_merged.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
bedops -u {input.neo} {input.chim} > {output}
"""


rule mergeContacts:
    input: "results/{genome}/{project}/junctions-view/{id}/All_{other}.bed"
    output: "results/{genome}/{project}/junctions-view/{id}/All_{other}_aggregated.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
cut -f 1,2,3 {input} | sort | uniq -c | sed -r 's/([0-9]) /\\1\\t/' |\
awk -v 'OFS=\\t' -v 'name_prefix=id' '{{print $2,$3,$4,name_prefix"_"NR,$1,"+"}}' |\
sort-bed - > {output}
"""


rule bedjToBed12:
    input: "results/{genome}/{project}/junctions-view/{id}/All_len40000_inIntrons_neoR_merged_aggregated.bed"
    output:  "results/{genome}/{project}/junctions-view/{id}/All_final.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
mkdir -p $(dirname {output})  
cat {input} | python workflow/scripts/bedj_to_bed12.py > {output}
"""


rule mergeProjectContacts:
    input: get_all_junctions
    output: "results/{genome}/{project}/views/global/junctions.bed"
    conda: "../envs/postprocess.yaml"
    shell: """
python workflow/scripts/bed12_merge_color.py {input} | sort-bed - | python workflow/scripts/bed_enumerate.py > {output}
"""
