import os



def check_msmcs(msmc_output, msmc_checked):
    if os.path.getsize(msmc_output) != 0:
        with open(msmc_checked, 'w'):
            pass
        print(f'{msmc_output} passed')
    else:
        print(f'WARNING: {msmc_output} is empty!! Not making {msmc_checked}')


def main():
    vcf_by_chrom = snakemake.input["vcf_by_chrom"]
    msmc_checked = snakemake.output["vcf_checked"]
    msmc_data = snakemake.params["msmc_data"]
    check_msmcs(msmc_output, msmc_checked, msmc_data)


if __name__ == "__main__":
    main()