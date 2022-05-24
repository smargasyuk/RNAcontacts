import glob

import pandas as pd
HUB_PATH = "results/trackhub"

samples = (
    pd.read_csv(config["samples"], sep="\t")
    .applymap(lambda x: x.strip() if isinstance(x, str) else x)
    .set_index("sample_name", drop=False)
)


def get_pass1_fq(wildcards):
    fq = samples.loc[(samples["sample_name"]  == wildcards["sample"]) & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"]), ["fq1", "fq2"]].iloc[0]
    return {
        "fq1": f"{fq.fq1}",
        "fq2": f"{fq.fq2}" 
    }


def get_pass2_fq(wildcards):
    fq = samples.loc[(samples["sample_name"] == wildcards["sample"]) & (samples.project == wildcards["project"]) &
     (samples.genome == wildcards["genome"]), ["fq1", "fq2"]].iloc[0]
    if wildcards["mate"] == "0":
        return f"{fq.fq1}"
    elif wildcards["mate"] == "1": 
        return f"{fq.fq2}" 
    else:
        raise Exception(f"Mate {wildcards['mate']} for sample {wildcards['sample']} not found")


def get_pass2_sj(wildcards):
    rows = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) 
        & (samples.genome == wildcards["genome"])]
    return [f"results/{row.genome}/{row.project}/bam/pass1/{row.sample_name}/SJ.out.tab" for row in rows.itertuples()]


def get_all_outputs(wildcards):
    bam = [f"results/{row.genome}/{row.project}/bam/pass2/{row.sample_name}_{mate}/Aligned.sortedByCoord.out.bam" for row in samples.itertuples() for mate in [0,1]] 
    contacts = [f"results/{row.genome}/{row.project}/contacts/{row.sample_name}/{jtype}.tsv" for row in samples.itertuples() for jtype in ["Neo", "Chimeric"]]
    global_contacts_view = [f"results/{row.genome}/{row.project}/views/global/contacts.bed" for row in samples.itertuples()]
    return bam + contacts + global_contacts_view + all_hub_files() 


def get_known_junctions(wildcards):
    control_samples = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()
    control_junctions = [f"results/{wildcards.genome}/{wildcards.project}/bam/pass1/{sample}/SJ.out.tab" for sample in control_samples]
    return {
        "control_jxn": control_junctions,
        "star_ref_dir": f"resources/star_genome/{wildcards.genome}"
    }


def get_all_junction_files(wildcards):
    relevant_samples = samples.loc[(samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()
    return [f"results/{wildcards['genome']}/{wildcards['project']}/junctions/{id}/{jtype}.tsv" for id in relevant_samples for jtype in ["Neo", "Chimeric"]]


def get_project_samples(wildcards):
    return samples.loc[(samples.treatment == "experiment") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()


def get_all_clusters(wildcards):
    relevant_samples = get_project_samples(wildcards)
    return [f"results/{wildcards['genome']}/{wildcards['project']}/clusters/{id}/{jtype}.tsv" for id in relevant_samples for jtype in ["Neo", "Chimeric"]]


def get_project_contact_bed_files(wildcards):
    relevant_samples = samples.loc[(samples.treatment == "experiment") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()
    return [f"results/{wildcards['genome']}/{wildcards['project']}/views/per_sample/{id}/contacts.bed" for id in relevant_samples]


def get_all_junctions(wildcards):
    relevant_samples = get_project_samples(wildcards)
    return [f"results/{wildcards['genome']}/{wildcards['project']}/junctions-view/{sample}/All_final.bed" for sample in relevant_samples]


def get_all_genomes():
    return samples.genome.unique().tolist()


def get_genome_hub_files(wildcards):
    relevant_projects = samples.loc[(samples.genome == wildcards["genome"])]["project"].unique().tolist()
    junctions = [f"{wildcards['hub_prefix']}/{wildcards['genome']}/{project}-junctions.bb" for project in relevant_projects]
    contacts = [f"{wildcards['hub_prefix']}/{wildcards['genome']}/{project}-contacts.bb" for project in relevant_projects]
    return contacts + junctions


def all_hub_files():
    tracks = [f'{HUB_PATH}/{g}/tracks.txt' for g in get_all_genomes()]
    static = [f'{HUB_PATH}/genomes.txt', f'{HUB_PATH}/hub.txt']
    return tracks + static

def get_genome_by_assembly(assembly):
    for g_name, g_dict in config['genomes'].items():
        if g_dict['assembly'] == assembly:
            return g_name

def get_assembly_by_genome(genome):
    return config['genomes'][genome]['assembly']
