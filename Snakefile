import pandas as pd
configfile: "./config.yaml"
from pathlib import Path
#snakemake --workflow-profile profiles/pheonix_slurm/ --use-conda --local-cores 10 --rerun-incomplete
#snakemake --use-conda --workflow-profile profiles/pheonix_slurm/
SAMPLES = pd.read_csv(config["samples"])
localrules: download_refgenome, prepare_pav_input

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
        #expand("projects/{project_id}/results/msmc2/{project_id}_msmc.png", project_id=SAMPLES['ccgp_project_id']),
        #expand("projects/{project_id}/results/msmc2/{project_id}_msmc_split.png", project_id=SAMPLES['ccgp_project_id']),
        #expand("projects/{project_id}/results/summary_statistics/mutation_summary.tsv", project_id=SAMPLES['ccgp_project_id']),
        expand("projects/{project_id}/results/summary_statistics/{project_id}_mut_histogram.png", project_id=SAMPLES['ccgp_project_id']),
        #expand("projects/{project_id}/results/msmc2/msmc_split.final.txt", project_id=SAMPLES['ccgp_project_id']),
        #expand("projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed.gz", project_id=SAMPLES['ccgp_project_id'])

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
        datasets download genome accession {wildcards.refGenome} --filename {params.dataset}
        7z x {params.dataset} -aoa -o{params.outdir} || unzip -o {params.dataset} -d {params.outdir}
        cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fna > {output.ref} && samtools faidx {output.ref} 
        """


rule prepare_pav_input:
    input:
        unpack(get_genomes)
    output:
        pav_config = "projects/{project_id}/config.json",
        assemblies = "projects/{project_id}/assemblies.tsv"
    params:
        project_id = "{project_id}",
        assembly_path = "projects/{project_id}/assemblies.tsv"
    script:
        "scripts/prepare_pav.py"

rule run_pav:
    input:
        pav_config = "projects/{project_id}/config.json",
        assemblies = "projects/{project_id}/assemblies.tsv",
    output:
        pav_vcf = "projects/{project_id}/{project_id}.vcf.gz",
        bed_zip = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed.gz",

    params:
        absdir = "$(realpath projects/{project_id})",
        homedir = "projects/{project_id}",
        native_pav_dir = "/private/groups/corbettlab/ryan/stairway/pav/Snakefile",
        pav_absdir = "projects/{project_id}/Snakefile",
        project_id = "{project_id}",

    log:
        "projects/{project_id}/logs/rule_run_pav.log"
    conda:
        "envs/pav.yml"
    benchmark:
        "benchmarks/{project_id}/run_pav.txt"
    shell:
        '''
        docker run --rm -v {params.absdir}:{params.absdir} --user "$(id -u):$(id -g)" --workdir {params.absdir} -e XDG_CACHE_HOME={params.absdir} becklab/pav:latest -c {threads} --notemp &> {log}
        '''


rule gunzip_bed:
    input:
        pav_vcf = "projects/{project_id}/{project_id}.vcf.gz",
        bed_zip = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed.gz",
    output:
        bed = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed"
    shell:
        '''
        gunzip -c {input.bed_zip} | tail -n +2 > {output.bed}
        '''

rule clean_vcf:
    input:
        pav_vcf = "projects/{project_id}/{project_id}.vcf.gz",
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
        bcftools view -Oz -o {output.cleaned_pav_vcf} {input.cleaned_pav_vcf_temp} --threads {threads}
        bcftools index {output.cleaned_pav_vcf} --threads {threads}
        '''

def get_primary_fai(wildcards):
    pid = wildcards.project_id
    primary = SAMPLES[SAMPLES["ccgp_project_id"] == pid]["primary_ref"].item()
    primary_path = Path("projects", pid, "genome/hap1", primary + ".fa.fai")

    return primary_path

checkpoint generate_chromosome_list:
    # Makes list of chromosomes WITH callable variants.
    input:
        cleaned_pav_vcf = "projects/{project_id}/pav_{project_id}_SNVs_cleaned.vcf.gz"
    output:
        chromosome_list = "projects/{project_id}/genome/chromosome_list.txt"
    shell:
        """
        bcftools query -f '%CHROM\n' {input.cleaned_pav_vcf} | sort | uniq -c | awk '$1 > 0 {{print $2}}' > {output.chromosome_list}
        """


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
        bcftools view -Oz -o {output.vcf_by_chromosome} {input.cleaned_pav_vcf} {wildcards.chromosome} --threads {threads}
        '''

rule msmc_input:
    input:
        vcf_by_chromosome = "projects/{project_id}/chromosome_vcfs/{chromosome}.vcf.gz",
        bed = "projects/{project_id}/results/{project_id}/callable/callable_regions_h1_500.bed"
    output:
        chromosome_hetsep_unfilt = "projects/{project_id}/results/unfiltered_msmc_input/{chromosome}_hetsep_unfilt.txt"
    shell:
        """
        msmc-tools/generate_multihetsep.py --chr {wildcards.chromosome} --mask {input.bed} {input.vcf_by_chromosome} > {output.chromosome_hetsep_unfilt}
        """

rule validate_hetsep:
    # Some msmc-tools/generate_multihetsep.py outputs create duplicate callable sites, with 0 callable regions.
    # This calls msmc2 to fail, so by adding this rule any regions with 0 callable sites are excluded.
    input:
        chromosome_hetsep_unfilt = "projects/{project_id}/results/unfiltered_msmc_input/{chromosome}_hetsep_unfilt.txt"
    output:
        hetsep_filtered = "projects/{project_id}/results/msmc_input/{chromosome}_hetsep.txt"
    log:
        filtered_out = "projects/{project_id}/logs/excluded_hetsep_sites/{chromosome}_excluded_sites.log"
    shell:
        """
        touch {output.hetsep_filtered}

        if [ -s {input.chromosome_hetsep_unfilt} ]; then
            awk '$3 > 0 {{ print > "{output.hetsep_filtered}" }} $3 <= 0 {{ print > "{log.filtered_out}" }}' {input.chromosome_hetsep_unfilt}
        fi
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
    
    print(f'# of chromosome vcf files with data: {len(all_files)}')

    return all_files

rule msmc2:
    input: get_files_with_data
    output:
        msmc2 = "projects/{project_id}/results/msmc2/msmc.final.txt",
        msmc2_split = "projects/{project_id}/results/msmc2/msmc_split.final.txt",
    params:
        prefix = "projects/{project_id}/results/msmc2/msmc"
    log:
        msmc2_log = "projects/{project_id}/logs/rule_msmc2.log",
        msmc2_log_split = "projects/{project_id}/logs/rule_msmc2_split.log"
    benchmark:
        "benchmarks/{project_id}/msmc2.txt"
    shell:
        '''
        ./msmc2_Linux -t {threads} -p 1*2+1*2+25*2+1*4+1*6 --outFilePrefix {params.prefix}_split {input} &> {log.msmc2_log_split}
        ./msmc2_Linux -t {threads} -p 1*4+25*2+1*4+1*6 --outFilePrefix {params.prefix} {input} &> {log.msmc2_log}
        '''

rule plot_msmc2:
    input:
        msmc2_split = "projects/{project_id}/results/msmc2/msmc_split.final.txt",
        msmc2 = "projects/{project_id}/results/msmc2/msmc.final.txt",
    output:
        plot_split = "projects/{project_id}/results/msmc2/{project_id}_msmc_split.png",
        plot = "projects/{project_id}/results/msmc2/{project_id}_msmc.png",
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
        plot = "projects/{project_id}/results/msmc2/{project_id}_msmc.png",
        plot_split = "projects/{project_id}/results/msmc2/{project_id}_msmc_split.png"
    output:
        summary_file = "projects/{project_id}/results/summary_statistics/mutation_summary.tsv"
    params:
        snv_snv = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/snv_snv.bed.gz",
        sv_inv = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/sv_inv.bed.gz",
        indel_del = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/svindel_del.bed.gz",
        indel_ins = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/svindel_ins.bed.gz",
        hetsep = "projects/{project_id}/results/msmc_input"
    script:
        "scripts/generate_summary.py"

rule generate_histogram:
    input:
        summary_file = "projects/{project_id}/results/summary_statistics/mutation_summary.tsv",
        chromosome_list = "projects/{project_id}/genome/chromosome_list.txt"
    output:
        histogram = "projects/{project_id}/results/summary_statistics/{project_id}_mut_histogram.png"
    params:
        snv_snv = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/snv_snv.bed.gz",
        sv_inv = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/sv_inv.bed.gz",
        indel_del = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/svindel_del.bed.gz",
        indel_ins = "projects/{project_id}/results/{project_id}/bed_hap/pass/h1/svindel_ins.bed.gz",
        out_dir = "projects/{project_id}/results/summary_statistics/",
        pid = '{project_id}'
    script:
        "scripts/generate_histograms.py"
    