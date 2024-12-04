import glob

import pandas as pd
HUB_PATH = "results/trackhub"
PREFIX = config["root_dir"]

# samples = (
#     pd.read_csv(config["samples"], sep="\t")
#     .applymap(lambda x: x.strip() if isinstance(x, str) else x)
#     .set_index("sample_name", drop=False)
# )

samples_table = pd.read_csv(config["samples"])

def get_pass1_fq(wildcards):
    fq = samples.loc[(samples["sample_name"]  == wildcards["sample"]) & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"]), ["fq1", "fq2"]].iloc[0]
    return {
        "fq1": f"{fq.fq1}",
        "fq2": f"{fq.fq2}" 
    }

def get_pass1_fq_for_adapterremoval(wildcards):
    fq = samples.loc[(samples["sample_name"]  == wildcards["sample"]) & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"]), ["fq1", "fq2"]].iloc[0]
    return [str(fq.fq1), str(fq.fq2)]

def get_pass2_fq(wildcards):
    fq = samples.loc[(samples["sample_name"] == wildcards["sample"]) & (samples.project == wildcards["project"]) &
     (samples.genome == wildcards["genome"]), ["fq1", "fq2"]].iloc[0]
    if wildcards["mate"] == "0":
        return f"{fq.fq1}"
    elif wildcards["mate"] == "1": 
        return f"{fq.fq2}" 
    else:
        raise Exception(f"Mate {wildcards['mate']} for sample {wildcards['sample']} not found")

def get_all_outputs(wildcards):
    bam = get_all_alignments(wildcards)
    contacts = [f"results/{row.genome}/{row.project}/contacts/{row.sample_name}/{jtype}.tsv.gz" for row in samples.itertuples() for jtype in ["Neo", "Chimeric"]]
    global_contacts_view = [f"results/{row.genome}/{row.project}/views/global/contacts.bed" for row in samples.itertuples()]
    return bam + contacts + global_contacts_view + all_hub_files() 

def get_all_alignments(wildcards):
    return [f"results/{row.genome}/{row.project}/bam/pass2/{row.sample_name}_{mate}/Aligned.sortedByCoord.out.bam" for row in samples.itertuples() for mate in [0,1]] 


def get_known_junctions(wildcards):
    control_samples = samples.loc[(samples.treatment == "control") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()
    control_junctions = [f"results/{wildcards.genome}/{wildcards.project}/bam/pass1/{sample}/SJ.out.tab" for sample in control_samples]
    return {
        "control_jxn": control_junctions,
        "star_ref_dir": f"resources/star_genome/{wildcards.genome}"
    }


def get_all_junctions_target(wildcards):
    return [f"results/{row.genome}/{row.project}/junctions/{row.sample_name}/{jtype}.tsv.gz" for row in samples.loc[(samples.treatment == "experiment")].itertuples() for jtype in ["Neo", "Chimeric"]]


def get_all_contacts_target(wildcards):
    return [f"results/{row.genome}/{row.project}/contacts/{row.sample_name}/{jtype}.tsv.gz" for row in samples.loc[(samples.treatment == "experiment")].itertuples() for jtype in ["Neo", "Chimeric"]]


def get_all_junction_files(wildcards):
    relevant_samples = samples.loc[(samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()
    return [f"results/{wildcards['genome']}/{wildcards['project']}/junctions/{id}/{jtype}.tsv.gz" for id in relevant_samples for jtype in ["Neo", "Chimeric"]]


def get_project_samples(wildcards):
    return samples.loc[(samples.treatment == "experiment") & (samples.project == wildcards["project"]) & (samples.genome == wildcards["genome"])]["sample_name"].to_list()


def get_all_clusters(wildcards):
    relevant_samples = get_project_samples(wildcards)
    return [f"results/{wildcards['genome']}/{wildcards['project']}/clusters/{id}/{jtype}.tsv.gz" for id in relevant_samples for jtype in ["Neo", "Chimeric"]]


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


def get_pass2_bam(wildcards):
    st2 = samples_table.loc[samples_table["Organism"] == "Homo sapiens"].reset_index()
    samples_pe = list(st2.loc[st2["LibraryLayout"] == "PAIRED"]['Run'].values)
    files_pe = [PREFIX + f"/{assembly}/bam/pe/pass2/{sra_id}_{mate}/Aligned.sortedByCoord.out.bam" for mate in [1,2] for sra_id in samples_pe for assembly in [config["assembly"]]]
    samples_se = list(st2.loc[st2["LibraryLayout"] == "SINGLE"]['Run'].values)
    files_se = [PREFIX + f"/{assembly}/bam/se/pass2/{sra_id}/Aligned.sortedByCoord.out.bam" for sra_id in samples_se for assembly in [config["assembly"]]]
    return files_pe + files_se

def get_pass2_sj(wildcards):
    rows = samples_table.copy()
    rows["LibraryLayout"] = rows["LibraryLayout"].replace({"PAIRED" : "pe", "SINGLE": 'se'})
    return [PREFIX + f"/{assembly}/bam/{row.LibraryLayout}/pass1/{row.Run}/SJ.out.tab"  for row in rows.itertuples() for assembly in [config["assembly"]]]

def split_samples_by_layout():
    st2 = samples_table.loc[samples_table["Organism"] == "Homo sapiens"].reset_index()
    samples_pe = list(st2.loc[st2["LibraryLayout"] == "PAIRED"]['Run'].values)
    samples_se = list(st2.loc[st2["LibraryLayout"] == "SINGLE"]['Run'].values)
    return samples_pe, samples_se

def get_all_sample_id_w_layout_idx():
    samples_pe, samples_se = split_samples_by_layout()
    return samples_se + [f"{s}_{l}" for s in samples_pe for l in [1,2]]

def get_junction_files(wildcards):
    sample_ids = get_all_sample_id_w_layout_idx()
    return [PREFIX + f"/{assembly}/RNAcontacts/junctions/{file_id}/{jtype}.tsv.gz"  for file_id in sample_ids for assembly in [config["assembly"]] for jtype in ["Neo", "Chimeric"]]

def get_sample_layout(file_id):

    samples_pe, samples_se = split_samples_by_layout()
    samples_pe = [f"{s}_{l}" for s in samples_pe for l in [1,2]]
    
    if file_id in samples_pe:
        return "PE"
    if file_id in samples_se:
        return "SE"
    raise KeyError("sample not found")
        
    
def get_bam_file_by_junction_file(wildcards):

    sample_layout = get_sample_layout(wildcards.file_id)
    
    if sample_layout == "PE":
        return PREFIX + f"/{wildcards['assembly']}/bam/pe/pass2/{wildcards['file_id']}/Aligned.sortedByCoord.out.bam"
    if sample_layout == "SE":
        return PREFIX + f"/{wildcards['assembly']}/bam/se/pass2/{wildcards['file_id']}/Aligned.sortedByCoord.out.bam"


def get_chim_file_by_junction_file(wildcards):

    sample_layout = get_sample_layout(wildcards["file_id"])
    
    if sample_layout == "PE":
        return PREFIX + f"/{wildcards['assembly']}/bam/pe/pass2/{wildcards['file_id']}/Chimeric.out.junction"
    if sample_layout == "SE":
        return PREFIX + f"/{wildcards['assembly']}/bam/se/pass2/{wildcards['file_id']}/Chimeric.out.junction"

    
