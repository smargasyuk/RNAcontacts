rule star_index:
    input:
        fasta="resources/star_genome_input/{genome}/genome.fasta",
        annotation="resources/star_genome_input/{genome}/genome.gtf",
    output:
        directory("resources/star_genome/{genome}/")
    params:
        extra = lambda wildcards, input: f"--sjdbGTFfile {input.annotation}"
    threads: 4
    log:
        "logs/{genome}/star_index_genome.log",
    wrapper:
        "0.84.0/bio/star/index"