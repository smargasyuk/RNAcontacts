rule get_genome:
    output:
        "resources/star_genome_input/{genome}/genome.fasta",
    log:
        "logs/{genome}/get-genome.log",
    params:
        species=lambda wildcards:config["genomes"][wildcards["genome"]]["species"],
        datatype="dna",
        build=lambda wildcards:config["genomes"][wildcards["genome"]]["build"],
        release=lambda wildcards:config["genomes"][wildcards["genome"]]["release"],
    cache: True
    wrapper:
        "0.77.0/bio/reference/ensembl-sequence"


rule get_annotation:
    output:
        "resources/star_genome_input/{genome}/genome.gtf",
    params: 
        species=lambda wildcards: config["genomes"][wildcards["genome"]]["species"],
        fmt="gtf",
        build=lambda wildcards:config["genomes"][wildcards["genome"]]["build"],
        release=lambda wildcards:config["genomes"][wildcards["genome"]]["release"],
        flavor="",    
    cache: True
    log:
        "logs/{genome}/get_annotation.log",
    wrapper:
        "0.77.0/bio/reference/ensembl-annotation"


rule star_index:
    input:
        fasta="resources/star_genome_input/{genome}/genome.fasta",
        annotation="resources/star_genome_input/{genome}/genome.gtf",
    output:
        directory("resources/star_genome/{genome}/")
    params:
        extra = lambda wildcards, input: f"--sjdbGTFfile {input.annotation}"
    cache: True
    threads: 4
    log:
        "logs/{genome}/star_index_genome.log",
    wrapper:
        "0.84.0/bio/star/index"