import pandas as pd

def generate_list(fai, output_list):
    fai_file_path = str(fai)
    output_path = str(output_list)
    print(f'FAI: {fai}')
    fai_file = pd.read_csv(fai_file_path, sep="\t", header=None, usecols=[0])

    chromosome_unique = set()
    for i, row in fai_file.iterrows():
        chrom = row[0]
        chromosome_unique.add(chrom)

    chromosome_final = list(chromosome_unique)


    with open(output_path, "w") as file:
        for chrom in chromosome_final:
            file.write(chrom + "\n")


def main():
    fai_file = snakemake.params["primary"]
    output_list = snakemake.output["chromosome_list"]

    generate_list(fai_file, output_list)

if __name__ == "__main__":
    main()