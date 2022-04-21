import glob

import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t")
    .applymap(lambda x: x.strip() if isinstance(x, str) else x)
    .set_index("sample", drop=False)
    .sort_index()
)


def get_pass1_fq(wildcards):
    fq = samples.loc[wildcards["sample"], ["fq1", "fq2"]]
    return {
        "fq1": f"{fq.fq1}",
        "fq2": f"{fq.fq2}" 
    }


def get_pass2_fq(wildcards):
    fq = samples.loc[wildcards["sample"], ["fq1", "fq2"]]
    if wildcards["mate"] == "0":
        return f"{fq.fq1}"
    elif wildcards["mate"] == "1": 
        return f"{fq.fq2}" 
    else:
        raise Exception(f"Mate {wildcards['mate']} for sample {wildcards['sample']} not found")


def get_pass2_sj(wildcards):
    rows = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) 
        & (samples.genome == wildcards["genome"])]
    return [f"results/{row.genome}/{row.project}/bam/pass1/{row.sample}/SJ.out.tab" for row in rows.itertuples()]


def get_all_outputs(wildcards):
    bam = [f"results/{row.genome}/{row.project}/bam/pass2/{row.sample}_{mate}/Aligned.sortedByCoord.out.bam" for row in samples.itertuples() for mate in [0,1]] 
    contacts = [f"results/{row.genome}/{row.project}/contacts/{row.sample}/{jtype}.tsv" for row in samples.itertuples() for jtype in ["Neo", "Chimeric"]]
    return bam + contacts


def get_known_junctions(wildcards):
    control_samples = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample"].to_list()
    return [f"results/{wildcards.genome}/{wildcards.project}/bam/pass1/{sample}/SJ.out.tab" for sample in control_samples] + [f"resources/star_genome/{wildcards.genome}/sjdbList.out.tab"]

def get_all_junction_files(wildcards):
    relevant_samples = samples.loc[(samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample"].to_list()
    return [f"results/{wildcards['genome']}/{wildcards['project']}/junctions/{id}/{jtype}.tsv" for id in relevant_samples for jtype in ["Neo", "Chimeric"]]