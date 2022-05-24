rule mergeKnownJunctionsPooled:
    input: unpack(get_pooled_known_junctions)
    output:
        "results/{genome}/global/pooled/all_sj.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
cut -f1,2,3 {input.control_jxn} {input.star_ref_dir}/sjdbList.out.tab | sort -u > {output}
"""

rule extractChimericJunctionsPooled:
    input: unpack(get_pooled_chim_files)
    output:
        "results/{genome}/global/pooled/junctions/Chimeric.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort --parallel=8 -S4G -k9,9 <(cat {input.mate0} | cut -f1-6,10 | awk -v OFS="\t" 'BEGIN{{s["+"]="-";s["-"]="+";}}{{print $4, $5, $5+1, s[$6], $1, $2, $2+1, s[$3], $7, 0}}') <(cat {input.mate1} | cut -f1-6,10 | awk -v OFS="\t" '{{print $1, $2, $2+1, $3, $4, $5, $5+1, $6, $7, 1}}') | grep -v GL > {output}
"""

rule extractNeoJunctionsPooled:
    input: unpack(get_pooled_bam_files)
    output:
        "results/{genome}/global/pooled/junctions/Neo.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort -k9,9  -S4G --parallel=8 <(echo -n {input.mate0} | parallel --verbose -d ' ' -r --keep-order "samtools view {{}}" | perl workflow/scripts/neo.pl {input.junctions} 0) \
           <(echo -n {input.mate1} | parallel --verbose -d ' ' -r --keep-order "samtools view {{}}" | perl workflow/scripts/neo.pl {input.junctions} 1) | grep -v GL |\
            awk -v 'OFS=\\t' '{{print $1,$2,$2+1,$3,$4,$5,$5+1,$6,$7,$8}}' > {output}       
"""

rule clusterDonorsAcceptorsPooled:
    input: "results/{genome}/global/pooled/junctions/Neo.tsv", "results/{genome}/global/pooled/junctions/Chimeric.tsv"
    output:
        donors = "results/{genome}/global/pooled/junctions/donors.bed",
        acceptors = "results/{genome}/global/pooled/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    conda: "../envs/postprocess.yaml"
    shell:
        """
cut -f1-4 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.donors}
cut -f5-8 {input} | sort-bed - | bedops -m --range {params.radius} - > {output.acceptors}
"""

rule clusterJunctionsPooled:
    input:
        junctions = "results/{genome}/global/pooled/junctions/{jtype}.tsv",
	    donors = "results/{genome}/global/pooled/junctions/donors.bed",
	    acceptors = "results/{genome}/global/pooled/junctions/acceptors.bed"
    params:
        radius = lambda wildcards: config['junctions_merge_radius']
    output:
        "results/{genome}/global/pooled/clusters/{jtype}.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
sort-bed {input.junctions} | \
intersectBed -a stdin -b {input.donors} -wa -wb | cut -f5- | sort-bed - | \
intersectBed -a stdin -b {input.acceptors} -wa -wb | cut -f5- | sort -k1,1 > {output}
"""

rule extractContactsPerExperimentPooled:
    input:
        "results/{genome}/global/pooled/clusters/{jtype}.tsv",
    output:
        "results/{genome}/global/pooled/contacts/{jtype}.tsv"
    conda: "../envs/postprocess.yaml"
    shell:
        """
mkdir -p $(dirname {output})
awk '{{n[$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8]++}}END{{for(j in n){{print j"\t"n[j]}}}}' {input} | sort > {output}
"""

rule viewContactsPerSamplePooled:
    input: 
        "results/{genome}/global/pooled/clusters/Neo.tsv",
        "results/{genome}/global/pooled/clusters/Chimeric.tsv"
    output:
       "results/{genome}/global/pooled/views/per_sample/contacts.bed"  
    conda: "../envs/postprocess.yaml"
    params:
        max_size = config['view_contacts_max_range']
    shell:
        """
mkdir -p $(dirname {output[0]})
cat {input} |\
cut -f 3- |\
awk -v 'OFS=\t' '$1==$4' |\
awk -v 'OFS=\t' '{{if ($2>$6){{s1=$2;s2=$3;$2=$5;$3=$6;$5=s1;$6=s2}};print}}' |\
awk -v 'OFS=\t' '$2<$6' |\
awk -v 'OFS=\t' '($3<$5)' |\
awk -v 'OFS=\t' '$6-$2<{params.max_size}' |\
awk -v 'OFS=\t' '{{n[$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6]++}}END{{for(j in n){{print j"\t"n[j]}}}}' |\
awk -v 'OFS=\t' '{{print $1,$2,$6,"id_"NR",reads="$7,$7,"+",$2,$6,"0,0,0",2,$3-$2","$6-$5,0","$5-$2}}' | sort-bed - > {output}
"""

rule viewContactsGlobalPooled:
    input: "results/{genome}/global/pooled/views/per_sample/contacts.bed" 
    output:
       "results/{genome}/global/pooled/views/global/contacts.bed"    
    conda: "../envs/postprocess.yaml"
    params:
        max_size = config['view_contacts_max_range']
    shell:
        """
mkdir -p $(dirname {output})
echo 'track name="pooled contacts" visibility=2 itemRgb="On"' > {output}
python workflow/scripts/bed12_merge_color.py {input} | \
sort-bed - |\
awk -v 'OFS=\t' '$4="id_"NR",count="$5' >> {output}
"""