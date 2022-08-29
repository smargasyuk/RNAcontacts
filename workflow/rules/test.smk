import filecmp

rule download_test_files:
    output: "RICseq_toy_data.tgz"
    shell: "wget https://zenodo.org/record/6475703/files/RICseq_toy_data.tgz"


rule extract_test_data:
    input: "RICseq_toy_data.tgz"
    output:
        fastq = expand("resources/fastq/test/{ftype}_toy_rep{i}r{j}.fastq.gz",
            ftype=["RNAseq", "RICseq"], i=[1,2], j=[0,1]),
        genome = expand("resources/star_genome_input/test_hg19/genome.{ext}", ext=["gtf", "fasta"]),
        ref_out = expand("resources/test_results/test_hg19/test/contacts/{sample}/{ftype}.tsv",
         sample=["RIC-seq_rep1", "RIC-seq_rep2"], ftype=["Chimeric", "Neo"])
    shell: "tar -xf RICseq_toy_data.tgz"


rule test:
    input:
        pipeline_out = expand("results/test_hg19/test/contacts/{sample}/{ftype}.tsv",
         sample=["RIC-seq_rep1", "RIC-seq_rep2"], ftype=["Chimeric", "Neo"]),
        ref_out = expand("resources/test_results/test_hg19/test/contacts/{sample}/{ftype}.tsv",
         sample=["RIC-seq_rep1", "RIC-seq_rep2"], ftype=["Chimeric", "Neo"])
    run:
        for f1, f2 in zip(input["pipeline_out"], input["ref_out"]):
            if not filecmp.cmp(f1, f2, shallow = False):
                print(f"Files {f1} and {f2} are different")
                raise ValueError
        print("Test passed")