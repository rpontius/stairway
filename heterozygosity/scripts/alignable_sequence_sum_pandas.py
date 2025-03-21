import pandas as pd

def sum_alignable_sites(bed, output_file):


    df = pd.read_csv(bed, sep="\t", header=None, usecols=[1, 2], names=["start", "stop"])
    length = (df["stop"] - df["start"]).sum()

    with open(output_file, "w") as output_file:
        output_file.write(str(length))

def main():
    bed = snakemake.input["bed_file"]
    output_file = snakemake.output["callable_sites"]

    sum_alignable_sites(bed, output_file)

if __name__ == "__main__":
    main()