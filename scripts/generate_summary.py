import pandas as pd
import gzip
import os



def tabulate_mutations(files, output, file_name, mutation_types):
    #with open(output, 'w') as out:
    data_output = []
    total_structural_variants = 0
    for file, name in zip(files, file_name):
        print(file)
        with gzip.open(file, 'rt') as f:
            df = pd.read_csv(f, delimiter='\t')
            total_mutations = len(df) # total mutations
            mut_bp_total = df['SVLEN'].sum()
            mut_bp_mean = df['SVLEN'].mean()
            mut_bp_std = df['SVLEN'].std()
            data_output.append([mutation_types.get(name, name), total_mutations, mut_bp_total, mut_bp_mean, mut_bp_std])
            if name in ["sv_inv", "indel_del", "indel_ins"]:
                total_structural_variants += total_mutations
    df_output = pd.DataFrame(data_output, columns=['Mutation_Type', 'Total_Mutations', 'Mutation_BP_Total', 'Mean', 'Std_Dev'])
    df_output.to_csv(output, sep='\t', index=False)

    return output, total_structural_variants

def hetsep_count(hetsep, output_file, structural_variants):
    
    total_lines = 0
    for filename in os.listdir(hetsep):
        filepath = os.path.join(hetsep, filename)
        if os.path.isfile(filepath):
            with open(filepath, 'r') as file:
                num_lines = sum(1 for line in file)
                print(f'{filename}: {num_lines}')
                total_lines += num_lines
        
    with open(output_file, 'a') as out:
        out.write(f"\nTotal lines in hetsep directory: {total_lines}\n")
        out.write(f"#Structural Variants: {structural_variants}\n")
                


def main():
    indel_del_f = snakemake.params["indel_del"]
    indel_ins_f = snakemake.params["indel_ins"]
    snv_snv_f = snakemake.params["snv_snv"]
    sv_inv_f = snakemake.params["sv_inv"]

    output = snakemake.output["summary_file"]

    hetsep = snakemake.params["hetsep"]

    file_list = [indel_del_f, indel_ins_f, snv_snv_f, sv_inv_f]
    file_name = ["indel_del", "indel_ins", "snv_snv", "sv_inv"]
    mutation_types = {"indel_del": "Deletion", "indel_ins": "Insertion", "sv_inv": "Inversion", "snv_snv": "SNV"}

    output_file, structural_variants = tabulate_mutations(file_list, output, file_name, mutation_types)
    hetsep_count(hetsep, output_file, structural_variants)

if __name__ == "__main__":
    main()

