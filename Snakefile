import pandas as pd
configfile: "./config.yaml"
from pathlib import Path

#snakemake -c 15 --use-conda
SAMPLES = pd.read_csv(config["samples"])

def get_genomes(wildcards):
    pid = wildcards.project_id
    primary = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["primary_ref"].item()
    alt = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["alternative_ref"].item()
    
    primary_path = expand("projects/{pid}/genome/hap1/{primary}.fa", pid=pid, primary=primary)
    alt_path = expand("projects/{pid}/genome/hap2/{alt}.fa", pid=pid, alt=alt)
     

    return {"hap1": primary_path, "hap2": alt_path}

def get_mutation_gen(wildcards):
    pid = wildcards.project_id
    mut = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["mutation_rate"].item()
    gen = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["generations"].item()

    return {"mutation_rate": mut, "generation": gen}

rule all:
    input:
        expand("projects/{project_id}/results/msmc2/{project_id}_msmc.png", project_id=SAMPLES['ccgp_project_id']),
        expand("projects/{project_id}/results/msmc2/{project_id}_msmc_split.png", project_id=SAMPLES['ccgp_project_id']),
        expand("projects/{project_id}/results/summary_statistics/summary.txt", project_id=SAMPLES['ccgp_project_id']),
        expand("projects/{project_id}/results/summary_statistics/{project_id}_mut_histogram.png", project_id=SAMPLES['ccgp_project_id']),
        expand("projects/{project_id}/results/msmc2/msmc_split.final.txt", project_id=SAMPLES['ccgp_project_id'])
        ###expand("projects/{project_id}/done_flags/hist_done.txt", project_id=SAMPLES['ccgp_project_id']),
        ###expand("projects/{project_id}/done_flags/rule_msmc2_done.txt",  project_id=SAMPLES['ccgp_project_id'])

rule download_refgenome:
    output:
        ref = "projects/{project_id}/genome/{hap}/{refGenome}.fa"
    params:
        dataset = "projects/{project_id}/genome/{hap}/{refGenome}_dataset.zip",
        outdir = "projects/{project_id}/genome/{hap}/{refGenome}"
    conda:
        "envs/refGenome.yml"
    log:
        "projects/{project_id}/logs/{refGenome}/{hap}/download_ref/log.txt"
    benchmark:
        "benchmarks/{project_id}/{refGenome}/{hap}/download_ref/benchmark.txt"
    shell:
        """
        mkdir -p {params.outdir}
        datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} \
        && 7z x {params.dataset} -aoa -o{params.outdir} || unzip -o {params.dataset} -d {params.outdir} \
        && cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref} && samtools faidx {output.ref}
        """


# rule download_refgenome:
#     output:
#         dataset = "projects/{project_id}/genome/{hap}/{refGenome}_dataset.zip"
#     params:
#         outdir = "projects/{project_id}/genome/{hap}/{refGenome}"
#     conda:
#         "envs/refGenome.yml"
#     log:
#         "projects/{project_id}/logs/{refGenome}/{hap}/download_ref/log.txt"
#     benchmark:
#         "benchmarks/{project_id}/{refGenome}/{hap}/download_ref/benchmark.txt"
#     shell:
#         """
#         mkdir -p {params.outdir}
#         datasets download genome accession --filename {output.dataset} {wildcards.refGenome}
#         """

# rule unzip_refgenome:
#     input:
#         dataset = "projects/{project_id}/genome/{hap}/{refGenome}_dataset.zip"
#     output:
#         ref = "projects/{project_id}/genome/{hap}/{refGenome}.fa",
#     params:
#         #dataset = "projects/{project_id}/genome/{hap}/{refGenome}_dataset.zip",
#         outdir = "projects/{project_id}/genome/{hap}/{refGenome}"
#     conda:
#         "envs/refGenome.yml"
#     log:
#         "projects/{project_id}/logs/{refGenome}/{hap}/unzip_ref/log.txt"
#     benchmark:
#         "benchmarks/{project_id}/{refGenome}/{hap}/unzip_ref/benchmark.txt"
#     shell:
#         """
#         7z x {input.dataset} -aoa -o{params.outdir} || unzip -o {input.dataset} -d {params.outdir} \
#         && cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref} && samtools faidx {output.ref}
#         """

rule prepare_pav_input:
    input:
        unpack(get_genomes)
    output:
        pav_config = "projects/{project_id}/config.json",
        assemblies = "projects/{project_id}/assemblies.tsv"
    params:
        project_id = "{project_id}",
        
    script:
        "scripts/prepare_pav.py"


rule run_pav:
    input:
        pav_config = "projects/{project_id}/config.json",
        assemblies = "projects/{project_id}/assemblies.tsv"
    output:
        pav_vcf = "projects/{project_id}/pav_{project_id}.vcf.gz",
        bed_zip = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed.gz"
    params:
        absdir = "$(realpath projects/{project_id})",
        homedir = "projects/{project_id}",
    log:
        "projects/{project_id}/logs/rule_run_pav.log"
    shell:
        '''
        docker run --rm -v {params.absdir}:{params.absdir} --user "$(id -u):$(id -g)" --workdir {params.absdir} -e XDG_CACHE_HOME={params.absdir} becklab/pav:latest -c 16 &> {log}
        '''

rule gunzip_bed:
    input:
        pav_vcf = "projects/{project_id}/pav_{project_id}.vcf.gz",
        bed_zip = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed.gz"
    output:
        bed = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed"
    # params:
    #     bed_zip = "projects/{project_id}/results/callable/callable_regions_h1_500.bed.gz"
    shell:
        '''
        gunzip -c {input.bed_zip} | tail -n +2 > {output.bed}
        '''

rule clean_vcf:
    input:
        pav_vcf = "projects/{project_id}/pav_{project_id}.vcf.gz"
    output:
        cleaned_pav_txt = temp("projects/{project_id}/pav_{project_id}_cleaned_data.txt"),
        pav_header = temp("projects/{project_id}/pav_{project_id}_cleaned_header.txt"),
        cleaned_pav_vcf_temp = temp("projects/{project_id}/pav_{project_id}_cleaned_temp.vcf"),

    shell:
    ## sets custom genome to 0/1 | updates 0/1 to 0|1 | queries all columns in vcf and makes all REF & ALT uppercase letters | removes lines that start with '#' | output cleaned_pav_vfc.
        '''
        bcftools +setGT {input.pav_vcf} -- -t a -n c:0/1 | \
        bcftools +setGT -- -t a -n p | \
        bcftools view -v snps | \
        bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT\n' | \
        awk -F'\t' '{{OFS="\t"; $4=toupper($4); $5=toupper($5); print }}' > {output.cleaned_pav_txt}

        zgrep -e "#" {input.pav_vcf} > {output.pav_header}

        cat {output.pav_header} {output.cleaned_pav_txt} > {output.cleaned_pav_vcf_temp}

        '''

    
rule index_pav:
    input:
        cleaned_pav_vcf_temp = "projects/{project_id}/pav_{project_id}_cleaned_temp.vcf"
    output:
        cleaned_pav_vcf = "projects/{project_id}/pav_{project_id}_SNVs_cleaned.vcf.gz"
    conda:
        "envs/bcftools.yml"
    shell:
        '''
        bcftools view -Oz -o {output.cleaned_pav_vcf} {input.cleaned_pav_vcf_temp}
        bcftools index {output.cleaned_pav_vcf}
        '''

def get_primary_fai(wildcards):
    pid = wildcards.project_id
    primary = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["primary_ref"].item()
    primary_path = Path("projects", pid, "genome/hap1", primary + ".fa.fai")
    # primary_path = expand("projects/{pid}/genome/hap1/{primary}.fa", pid=pid, primary=primary)

    return primary_path

checkpoint generate_chromosome_list:
    input:
        cleaned_pav_vcf = "projects/{project_id}/pav_{project_id}_SNVs_cleaned.vcf.gz"
    output:
        chromosome_list = "projects/{project_id}/genome/chromosome_list.txt"
    params:
        primary = get_primary_fai
    script:
        "scripts/get_chromosomes.py"


def get_chromosomes(wc):
    pid = wc.project_id
    chromosomes = checkpoints.generate_chromosome_list.get(**wc).output[0]

    chromosome_list = Path(chromosomes)


    all_chromosomes = []
    with open(chromosome_list, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                all_chromosomes.append(line)
    
    print(len(all_chromosomes))
    out = [*expand("projects/{project_id}/chromosome_vcfs/{chromosome}.vcf.gz", **wc, chromosome=all_chromosomes),
    *expand("projects/{project_id}/results/msmc_input/{chromosome}_hetsep.txt", **wc, chromosome=all_chromosomes)]

    return out

rule split_vcf:
    input:
        cleaned_pav_vcf = "projects/{project_id}/pav_{project_id}_SNVs_cleaned.vcf.gz",
    output:
        vcf_by_chromosome = temp("projects/{project_id}/chromosome_vcfs/{chromosome}.vcf.gz")
    conda:
        "envs/bcftools.yml"
    shell:
        '''
        bcftools view -Oz -o {output.vcf_by_chromosome} {input.cleaned_pav_vcf} {wildcards.chromosome}
        '''

rule msmc_input:
    input:
        vcf_by_chromosome = "projects/{project_id}/chromosome_vcfs/{chromosome}.vcf.gz",
        bed = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed"
    output:
        msmc_output = "projects/{project_id}/results/msmc_input/{chromosome}_hetsep.txt"
    shell:
        """
        msmc-tools/generate_multihetsep.py --chr {wildcards.chromosome} --mask {input.bed} {input.vcf_by_chromosome} > {output.msmc_output}
        """


rule get_chromosome_files:
    input: get_chromosomes
    output: touch("projects/{project_id}/genome/chrom.done")


    
checkpoint get_files:
    input:
        chromosome_list = "projects/{project_id}/genome/chromosome_list.txt",
        done = "projects/{project_id}/genome/chrom.done"
    output:
        done_flag = "projects/{project_id}/results/hetseps_used.txt"
    params:
        msmc_path = "projects/{project_id}/results/msmc_input"
    script:
        "scripts/check_file_size.py"

def get_files_with_data(wc):
    pid = wc.project_id
    msmc_files = checkpoints.get_files.get(**wc).output[0]

    msmc_file_list = Path(msmc_files)


    all_files = []
    with open(msmc_file_list, "r") as file:
        for line in file:
            line = line.strip()
            if line:
                all_files.append(line)
    
    print(len(all_files))

    return all_files

rule msmc2:
    input: get_files_with_data
    output:
        msmc2 = "projects/{project_id}/results/msmc2/msmc.final.txt",
        msmc2_split = "projects/{project_id}/results/msmc2/msmc_split.final.txt",
        # rule_msmc2_done = "projects/{project_id}/done_flags/rule_msmc2_done.txt"
    params:
        prefix = "projects/{project_id}/results/msmc2/msmc"
    threads: 18
    log:
        msmc2_log = "projects/{project_id}/logs/rule_msmc2.log",
        msmc2_log_split = "projects/{project_id}/logs/rule_msmc2_split.log"
    shell:
        '''
        ./msmc2_Linux -t {threads} -p 1*4+25*2+1*4+1*6 --outFilePrefix {params.prefix} {input} &> {log.msmc2_log}
        ./msmc2_Linux -t {threads} -p 1*2+1*2+25*2+1*4+1*6 --outFilePrefix {params.prefix}_split {input} &> {log.msmc2_log_split}
        '''
        #touch {output.rule_msmc2_done}

rule plot_msmc2:
    input:
        msmc2 = "projects/{project_id}/results/msmc2/msmc.final.txt",
        msmc2_split = "projects/{project_id}/results/msmc2/msmc_split.final.txt",
        # rule_msmc2_done = "projects/{project_id}/done_flags/rule_msmc2_done.txt"
    output:
        plot = "projects/{project_id}/results/msmc2/{project_id}_msmc.png",
        plot_split = "projects/{project_id}/results/msmc2/{project_id}_msmc_split.png",
    params:
        pid = '{project_id}',
        mut = lambda wc: SAMPLES[SAMPLES["ccgp_project_id"] == wc.project_id]["mutation_rate"].item(),
        gen = lambda wc: SAMPLES[SAMPLES["ccgp_project_id"] == wc.project_id]["generations"].item()
    conda:
        "envs/plotting.yml"
    script:
        "scripts/plot_msmc.R"

rule generate_summary:
    input:
        plot = "projects/{project_id}/results/msmc2/{project_id}_msmc.png"
    output:
        summary_file = "projects/{project_id}/results/summary_statistics/summary.txt"
    params:
        indel_del = "projects/{project_id}/results/{project_id}/bed/indel_del.bed.gz",
        indel_ins = "projects/{project_id}/results/{project_id}/bed/indel_ins.bed.gz",
        snv_snv = "projects/{project_id}/results/{project_id}/bed/snv_snv.bed.gz",
        sv_del = "projects/{project_id}/results/{project_id}/bed/sv_del.bed.gz",
        sv_ins = "projects/{project_id}/results/{project_id}/bed/sv_ins.bed.gz",
        sv_inv = "projects/{project_id}/results/{project_id}/bed/sv_inv.bed.gz",
        hetsep = "projects/{project_id}/results/msmc_input"
    script:
        "scripts/generate_summary.py"

rule generate_histogram:
    input:
        summary_file = "projects/{project_id}/results/summary_statistics/summary.txt",
        chromosome_list = "projects/{project_id}/genome/chromosome_list.txt"
    output:
        #done_file = "projects/{project_id}/done_flags/hist_done.txt",
        histogram = "projects/{project_id}/results/summary_statistics/{project_id}_mut_histogram.png"
    params:
        indel_del = "projects/{project_id}/results/{project_id}/bed/indel_del.bed.gz",
        indel_ins = "projects/{project_id}/results/{project_id}/bed/indel_ins.bed.gz",
        snv_snv = "projects/{project_id}/results/{project_id}/bed/snv_snv.bed.gz",
        sv_del = "projects/{project_id}/results/{project_id}/bed/sv_del.bed.gz",
        sv_ins = "projects/{project_id}/results/{project_id}/bed/sv_ins.bed.gz",
        sv_inv = "projects/{project_id}/results/{project_id}/bed/sv_inv.bed.gz",
        out_dir = "projects/{project_id}/results/summary_statistics/",
        pid = '{project_id}'
    script:
        "scripts/generate_histograms.py"
    