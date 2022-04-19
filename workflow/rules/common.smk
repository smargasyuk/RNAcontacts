import glob

import pandas as pd

samples = (
    pd.read_csv(config["samples"], sep="\t")
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
    rows = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]
    return [f"results/{row.genome}/{row.project}/bam/pass1/{row.sample}/SJ.out.tab" for row in rows.itertuples()]