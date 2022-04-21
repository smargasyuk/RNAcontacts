.PHONY : test download

download: config/samples.tsv resources/star_genome_input/test_hg19/genome.gtf resources/star_genome_input/test_hg19/genome.fasta
	### all test data are on place ###

test:
	snakemake --cores 8 --use-conda
	for f in `ls results/test_hg19/test/contacts/`; do cmp results/test_hg19/test/contacts/$$f/Neo.tsv resources/test_results/test_hg19/test/contacts/$$f/Neo.tsv; done
	for f in `ls results/test_hg19/test/contacts/`; do cmp results/test_hg19/test/contacts/$$f/Chimeric.tsv resources/test_results/test_hg19/test/contacts/$$f/Chimeric.tsv; done
	### all tests completed successfully ###


config/samples.tsv resources/star_genome_input/test_hg19/genome.gtf resources/star_genome_input/test_hg19/genome.fasta:
	wget http://arkuda.skoltech.ru/~smargasyuk/RIC-contacts/RICseq_toy_data_22.04.2022.tgz
	tar -xf RICseq_toy_data_22.04.2022.tgz
