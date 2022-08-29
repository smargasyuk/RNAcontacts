# General settings

To configure this workflow, modify ``config/config.yaml`` according to your needs, following the explanations provided in the file.

# Sample sheet

Add samples to `config/samples.tsv`. For each sample, the following columns have to be defined:

+ `sample_name`: unique sample name;
+ `project`: unique project name. Control junctions are merged across all project samples, experiment contact clustering and final view generation are also performed for the whole project; 
+ `genome`: genome name. Put your genome fasta sequence to `resources/star_genome_input/{genome}/genome.fasta` and GTF transcript annotation to `resources/star_genome_input/{genome}/annotation.gtf`, or define the assembly in `config.yaml` block `genomes.{genome}` for automated download with [genomepy](https://github.com/vanheeringen-lab/genomepy).
+ `treatment`: 'experiment' or 'control'. Currently we use RNA-Seq data from the same cell line as a control;
+ `fq1`, `fq2`: paired-end fastq read files.
