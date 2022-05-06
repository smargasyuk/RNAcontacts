wildcard_constraints:
    genome="[^\/]+",
    fname="[^\/]+"

HUB_PATH = "results/trackhub"

hub_template = """
track {0}
type bigBed 12 +
bigDataUrl {0}.bb
itemRgb	on
visibility 3
shortLabel {0}
longLabel {0}
"""

genome_template = """
genome {0}
trackDb {0}/tracks.txt
"""


rule copyJunctions:
    input: "results/{genome}/{project}/views/global/junctions.bed"
    output: '{hub_prefix}/input/raw/{genome}/{project}-junctions.bed'
    shell: """
mkdir -p $(dirname {output})
cp {input} {output}
"""


rule copyContacts:
    input: "results/{genome}/{project}/views/global/contacts.bed"
    output: '{hub_prefix}/input/raw/{genome}/{project}-contacts.bed'
    shell: """
mkdir -p $(dirname {output})
tail -n +2 {input} > {output}
"""

rule bedClip:
    input: 
        bed='{hub_prefix}/input/raw/{genome}/{fname}.bed',
        chromsizes='{hub_prefix}/input/raw/{genome}/chromSizes',
    output:
        bed='{hub_prefix}/input/clipped/{genome}/{fname}.bed'
    conda: "../envs/hub.yaml"
    shell: """
mkdir -p $(dirname {output})
bedClip {input.bed} {input.chromsizes} {output}
"""


rule fetchChromSizes:
    output: '{hub_prefix}/input/raw/{genome}/chromSizes'
    conda: "../envs/hub.yaml"
    shell: """
mkdir -p $(dirname {output})
fetchChromSizes {wildcards.genome} > {output}
"""


rule sortBed:
    input: '{hub_prefix}/input/clipped/{genome}/{fname}.bed',
    output: '{hub_prefix}/input/sorted/{genome}/{fname}.bed'
    conda: "../envs/hub.yaml"
    shell: """
mkdir -p $(dirname {output})
awk -v 'OFS=\\t' '$5=$5<999?$5:999' {input} | sort-bed - > {output}    
"""


rule bed2bigbed:
    input: 
        bed='{hub_prefix}/input/sorted/{genome}/{fname}.bed',
        chromsizes='{hub_prefix}/input/raw/{genome}/chromSizes',
    output: '{hub_prefix}/{genome}/{fname}.bb'
    conda: "../envs/hub.yaml"
    shell: """
mkdir -p $(dirname {output})
bedToBigBed {input.bed} {input.chromsizes} {output}
"""		
        
                    
rule hub_files:
    input: get_genome_hub_files
    output: '{hub_prefix}/{genome}/tracks.txt'
    run:
        with open(output[0], 'w') as handle:
            for i in input:
                handle.write(hub_template.format(Path(i).stem))


rule genomes_description:
    input: config['samples']
    output: f'{HUB_PATH}/genomes.txt'
    run:
        with open(output[0], 'w') as handle:
            for g in get_all_genomes():
                handle.write(genome_template.format(g))


rule hub_description:
    input: 'resources/trackhub/hub.txt'
    output: f'{HUB_PATH}/hub.txt'
    shell: """
mkdir -p $(dirname {output})
cp {input} {output}
"""
