import pandas as pd
import os
import json

def generate_config(pav_config, project_id, ref_path):
    # config_output = f"{pav_config}"
    # ref_path = f"/scratch2/erik/stairway/{ref_path}"
    ref_path_string = str(ref_path)
    genome_path_split = ref_path_string.split("/genome/", 1)
    genome_path = "genome/" + genome_path_split[1]
    print(genome_path)
    json_dict = {"reference": genome_path}
    json_string = json.dumps(json_dict)


    with open(pav_config, "w") as file:
        file.write(json_string)
        

def generate_assemblies(assemblies, project_id, ref_path):
    # assemblies_output = f"/scratch2/erik/stairway/{assemblies}"
    # ref_path = f"/scratch2/erik/stairway/{ref_path}"
    ref_path_string = str(ref_path)
    genome_path_split = ref_path_string.split("/genome/", 1)
    genome_path = "genome/" + genome_path_split[1]

    header = "NAME\tHAP1\tHAP2"
    data = f"{project_id}\t{genome_path}\t"

    with open(assemblies, "w") as file:
        file.write(header + "\n")
        file.write(data)


def main():
    pav_config = snakemake.output["pav_config"]
    assemblies = snakemake.output["assemblies"]

    ref_path_primary = snakemake.input["hap1"]
    ref_path_alt = snakemake.input["hap2"]

    project_id = snakemake.params["project_id"]

    generate_config(pav_config, project_id, ref_path_primary)
    generate_assemblies(assemblies, project_id, ref_path_alt)


    

if __name__ == "__main__":
    main()