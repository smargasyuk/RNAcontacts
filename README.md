# RIC-contacts

Prediction of RNA-RNA contacts from RIC-seq data. Developed by Sergei Margasyuk (smargasyuk@gmail.com) and Dmitri Pervouchine (pervouchine@gmail.com).

## Description

This package contains a pipeline for prediction of RNA-RNA contacts from RIC-seq data (Cai et al., [2020](https://doi.org/10.1038/s41586-020-2249-1)). 

## Usage

### Step 1: Obtain a copy of this workflow

[Clone](https://help.github.com/en/articles/cloning-a-repository) this repository to your local system, into the place where you want to perform the data analysis.
    
    git clone https://github.com/pervouchine/RIC-contacts.git
    cd RIC-contacts

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the files in the `config/` folder. Adjust `config.yaml` to configure the workflow execution, and `samples.tsv` to specify your sample setup.

### Step 3: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).

### Step 4: Execute workflow

Activate the conda environment:

    conda activate snakemake

Test your configuration by performing a dry-run via

    snakemake --use-conda -n

Execute the workflow locally via

    snakemake --use-conda --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --use-conda --cluster qsub --jobs 100

or

    snakemake --use-conda --drmaa --jobs 100

See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

#### Test run

To make a test run, type

```
make download
make test
```

The script will download a toy dataset (sample sheet, truncated fastq files, genome, and genome annotation confined to the first 100MB of chr1), unpack, 
update the config file, and execute the pipeline. The output files in `results/test_hg19/test/contacts` will be compared to those provided in the archive.

#### Run on full RIC-seq data

Download the RIC-seq files for HeLa cell line from GEO repository [GSE127188](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127188). Download the control RNA-seq files from ENCODE consortium [webpage](https://www.encodeproject.org/).

The files are as follows:
```
  RNASeq_HeLa_total_rep1:
    - fastq/ENCFF000FOM.fastq
    - fastq/ENCFF000FOV.fastq
  RNASeq_HeLa_total_rep2:
    - fastq/ENCFF000FOK.fastq
    - fastq/ENCFF000FOY.fastq
  RIC-seq_HeLa_rRNA_depleted_rep1:
    - fastq/SRR8632820_1.fastq
    - fastq/SRR8632820_2.fastq
  RIC-seq_HeLa_rRNA_depleted_rep2:
    - fastq/SRR8632821_1.fastq
    - fastq/SRR8632821_2.fastq
```


### Step 5: Investigate results

The output of the pipeline consists of the following files:
 + `results/{genome}/{project}/{sample}/contacts` is the list of contacts and their respective read counts in tsv format (columns 1-3 and 4-6 are the contacting coordinates, column 7 is read count). 
 + `results/{genome}/{project}/views/global/contacts.bed` is the BED12 file with all contacts with length less than 






