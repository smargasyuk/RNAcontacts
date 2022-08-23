rule download_test_files:
    output: "RICseq_toy_data.tgz"
    shell: "wget https://zenodo.org/record/6475703/files/RICseq_toy_data.tgz"


rule extract_test_data:
    input: "RICseq_toy_data.tgz"
    output:
        fastq = expand("resources/fastq/test/{ftype}_toy_rep{i}r{j}.fastq.gz",
            ftype=["RNAseq", "RICseq"], i=[1,2], j=[0,1]),
        genome = expand("resources/star_genome_input/test_hg19/genome.{ext}", ext=["gtf", "fasta"]),
        ref_output = directory("resources/test_results")
    shell: "tar -xf RICseq_toy_data.tgz"


rule test:
    input:
        expand("results/test_hg19/test/contacts/{sample}/{ftype}.tsv",
         sample=["RIC-seq_rep1", "RIC-seq_rep2"], ftype=["Chimeric", "Neo"])
    shell: """
for i in "${{{input}[@]}}"; do 
  echo "$i"
done
# for f in `ls results/test_hg19/test/contacts/`; do cmp results/test_hg19/test/contacts/$$f/Neo.tsv resources/test_results/test_hg19/test/contacts/$$f/Neo.tsv; done
# for f in `ls results/test_hg19/test/contacts/`; do cmp results/test_hg19/test/contacts/$$f/Chimeric.tsv resources/test_results/test_hg19/test/contacts/$$f/Chimeric.tsv; done
"""