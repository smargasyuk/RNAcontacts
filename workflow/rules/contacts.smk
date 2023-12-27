rule mergeKnownJunctions:
    input: unpack(get_known_junctions)
    output:
        "results/{genome}/{project}/all_sj.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input.control_jxn} {input.star_ref_dir}/sjdbList.out.tab | sort -u > {output}
"""

rule extractChimericJunctions:
    input:
        mate0 = "results/{genome}/{project}/bam/pass2/{id}_0/Chimeric.out.junction",
        mate1 = "results/{genome}/{project}/bam/pass2/{id}_1/Chimeric.out.junction"
    output:
        "results/{genome}/{project}/junctions/{id}/Chimeric.tsv.gz"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort --parallel=8 -S4G -k9,9 \
           <(cut -f1-6,10 {input.mate0} | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') \
           <(cut -f1-6,10 {input.mate1} | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | \
           grep -v GL | pigz > {output}
"""

rule extractNeoJunctions:
    input:
        mate0 = "results/{genome}/{project}/bam/pass2/{id}_0/Aligned.sortedByCoord.out.bam",
	mate1 = "results/{genome}/{project}/bam/pass2/{id}_1/Aligned.sortedByCoord.out.bam",
        junctions = "results/{genome}/{project}/all_sj.tsv"
    output:
        "results/{genome}/{project}/junctions/{id}/Neo.tsv.gz"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9  -S4G --parallel=8 <(samtools view {input.mate0} | perl workflow/scripts/neo.pl {input.junctions} 0) \
           <(samtools view {input.mate1} | perl workflow/scripts/neo.pl {input.junctions} 1) | grep -v GL |\
            awk -v 'OFS=\\t' '{{print $1,$2,$2+1,$3,$4,$5,$5+1,$6,$7,$8}}' | pigz > {output}       
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
unpigz -c {input} | cut -f1-4 | sort-bed - | bedops -m --range {params.radius} - > {output.donors}
unpigz -c {input} | cut -f5-8 | sort-bed - | bedops -m --range {params.radius} - > {output.acceptors}
"""

rule clusterJunctions:
    input:
        junctions = "results/{genome}/{project}/junctions/{id}/{jtype}.tsv.gz",
        donors = "results/{genome}/{project}/junctions/donors.bed",
        acceptors = "results/{genome}/{project}/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "results/{genome}/{project}/clusters/{id}/{jtype}.tsv.gz"
    conda: "../envs/postprocess.yaml"
    resources:
        mem_mb=20000
    shell:
        """
mkdir -p $(dirname {output})
unpigz -c {input.junctions} |\
sort-bed - |\
intersectBed -a stdin -b {input.donors} -wa -wb | cut -f5- | sort-bed - |\
intersectBed -a stdin -b {input.acceptors} -wa -wb | cut -f5- | sort -k1,1 | pigz > {output}
"""

rule extractContactsPerExperiment:
    input:
        "results/{genome}/{project}/clusters/{id}/{jtype}.tsv.gz",
    output:
        "results/{genome}/{project}/contacts/{id}/{jtype}.tsv.gz"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
unpigz -c {input} |\
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8]++}}END{{for(j in n){{print j"\t"n[j]}}}}' |\
sort | pigz > {output}
"""

rule viewContactsPerSample:
    input: 
        "results/{genome}/{project}/clusters/{id}/Neo.tsv.gz",
        "results/{genome}/{project}/clusters/{id}/Chimeric.tsv.gz"
    output:
       "results/{genome}/{project}/views/per_sample/{id}/contacts.bed"  
    conda: "../envs/postprocess.yaml"
    params:
        max_size = config['view_contacts_max_range']
    shell:
        """
mkdir -p $(dirname {output[0]})
unpigz -c {input} |\
cut -f 3- |\
awk -v 'OFS=\t' '$1==$4' |\
awk -v 'OFS=\t' '{{if ($2>$6){{s1=$2;s2=$3;$2=$5;$3=$6;$5=s1;$6=s2}};print}}' |\
awk -v 'OFS=\t' '$2<$6' |\
awk -v 'OFS=\t' '($3<$5)' |\
awk -v 'OFS=\t' '$6-$2<{params.max_size}' |\
awk -v 'OFS=\t' '{{n[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]++}}END{{for(j in n){{print j"\t"n[j]}}}}' |\
awk -v 'OFS=\t' '{{print $1,$2,$6,"id_"NR",reads="$7,$7,"+",$2,$6,"0,0,0",2,$3-$2","$6-$5,0","$5-$2}}' | sort-bed - > {output}
"""

rule viewContactsGlobal:
    input: get_project_contact_bed_files
    output:
       "results/{genome}/{project}/views/global/contacts.bed"    
    conda: "../envs/postprocess.yaml"
    params:
        max_size = config['view_contacts_max_range']
    shell:
        """
mkdir -p $(dirname {output})
echo 'track name="{wildcards.project} contacts" visibility=2 itemRgb="On"' > {output}
python workflow/scripts/bed12_merge_color.py {input} | \
sort-bed - |\
awk -v 'OFS=\t' '$4="id_"NR",count="$5' >> {output}
"""