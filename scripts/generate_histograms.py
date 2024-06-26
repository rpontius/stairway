import pandas as pd
import matplotlib.pyplot as plt
import gzip
import os


def generate_histograms(files, done_file, file_names, out_dir, project_id, histogram_out):
    combined_data = pd.DataFrame(columns=['#CHROM', 'SVLEN', 'File']) 
    for file, file_name in zip(files, file_names):
        with gzip.open(file, 'rt') as f:
            df = pd.read_csv(f, delimiter='\t', usecols=['#CHROM', 'SVLEN'])
            df['File'] = file_name  # Add a 'File' column to identify the source file
            combined_data = combined_data._append(df, ignore_index=True)  # Append filtered data to combined_data
    # Plotting the histograms
    if not combined_data.empty:
        plt.figure(figsize=(10, 6))
        bins = range(0, 201, 2)
        plt.hist([combined_data[combined_data['File'] == file_name]['SVLEN'] for file_name in file_names],
                bins=bins, stacked=True, alpha=0.7, label=file_names)
        plt.xlabel('SVLEN')
        plt.ylabel('# Mutations')
        plt.legend()
        plt.title(f'{project_id}')
        
        # fig_name = f'{project_id}_histogram.png'
        # output_hist = os.path.join(out_dir, fig_name)
        plt.savefig(histogram_out)   # Save the plot to the specified file
        print(f'Finished {project_id}')
    #### with open(done_file, 'w') as out:
    ####     out.write('')
    
def main():
    indel_del_f = snakemake.params["indel_del"]
    indel_ins_f = snakemake.params["indel_ins"]
    #snv_snv_f = snakemake.params["snv_snv"]
    sv_del_f = snakemake.params["sv_del"]
    sv_ins_f = snakemake.params["sv_ins"]
    sv_inv_f = snakemake.params["sv_inv"]
   ##### done_file = snakemake.output["done_file"]
    done_file = "Nothing"
    histogram_out = snakemake.output["histogram"]
    chromosome_file = snakemake.input["chromosome_list"]
    out_dir = snakemake.params["out_dir"]
    project_id = snakemake.params["pid"]
    file_list = [indel_del_f, indel_ins_f, sv_del_f, sv_ins_f, sv_inv_f]
    file_name = ["indel_del", "indel_ins", "sv_del", "sv_ins", "sv_inv"]

    # chromosomes = get_chromosomes(chromosome_file)

    generate_histograms(file_list, done_file, file_name, out_dir, project_id, histogram_out)

if __name__ == "__main__":
    main()