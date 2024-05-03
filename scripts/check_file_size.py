import os
from pathlib import Path


def check_msmcs(chromosome_list, output, msmc_path):
    chroms = []
    with open(chromosome_list, 'r') as file:
        for chrom in file:
            print(chrom)
            chroms.append(chrom.strip())

    
    print(msmc_path)
    with open(output, "w") as file:
        for chromosome in chroms:
            het_path = Path(msmc_path, chromosome + "_hetsep.txt")
            het_name = str(chromosome) + "_hetsep.txt"
            #print(het_path)
            if het_path.is_file():
                if os.path.getsize(het_path) != 0:
                    file.write(str(het_path) + "\n")
                else:
                    print(f'{het_name} is EMPTY !!')
            else:
                print(f'{het_path} is NOT a directory !!')



def main():
    chromosome_list = snakemake.input["chromosome_list"]
    output = snakemake.output["done_flag"]
    msmc_path = snakemake.params["msmc_path"]
    check_msmcs(chromosome_list, output, msmc_path)


if __name__ == "__main__":
    main()