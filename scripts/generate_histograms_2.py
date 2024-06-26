import pandas as pd
import matplotlib.pyplot as plt
import gzip
import os



def get_chromosomes(chromosome_file):
    chromosome_list = []
    with open(chromosome_file, 'r') as file:
        for chromosome in file:
            chromosome_list.append(chromosome.strip())

    return chromosome_list

def generate_histograms(files, done_file, file_names, chromosomes, out_dir):
    with open(done_file, 'w') as out:
        for chromosome in chromosomes:
            combined_data = pd.DataFrame(columns=['#CHROM', 'SVLEN'])
            for file, file_name in zip(files, file_names):
                with gzip.open(file, 'rt') as f:
                    df = pd.read_csv(f, delimiter='\t')
                    chrom_data = df[df['#CHROM'] == chromosome] 
                    # chrom_data['File'] = file_name
                    combined_data_frame = combined_data._append(chrom_data, ignore_index=True)
            
            print(combined_data)
    
def main():
    indel_del_f = snakemake.params["indel_del"]
    indel_ins_f = snakemake.params["indel_ins"]
    snv_snv_f = snakemake.params["snv_snv"]
    sv_del_f = snakemake.params["sv_del"]
    sv_ins_f = snakemake.params["sv_ins"]
    sv_inv_f = snakemake.params["sv_inv"]
    done_file = snakemake.output["done_file"]
    chromosome_file = snakemake.input["chromosome_list"]
    out_dir = snakemake.params["out_dir"]
    file_list = [indel_del_f, indel_ins_f, snv_snv_f, sv_del_f, sv_ins_f, sv_inv_f]
    file_name = ["indel_del", "indel_ins", "snv_snv", "sv_del", "sv_ins", "sv_inv"]

    chromosomes = get_chromosomes(chromosome_file)
    generate_histograms(file_list, done_file, file_name, chromosomes, out_dir)

if __name__ == "__main__":
    main()