## Note: This README is still a work in progress.

## Overview

This workflow is designed to infer an organisms effective population size through time using the Multiple Sequentially Markovian Coalescent 2 (MSMC2) method. This software requires access to an NCBI reference genome (both primary and alternative accessions are needed).

## Dependencies

- Python 3.10
- Snakemake (7.32.3 reccomended)
- Mamba (for environment management)
- Access to an NCBI reference genome. Both haplotypes (primary & alternative) accession ID's are needed.

## Inputs and Outputs
- Inputs:
	i) 'samples.csv'
- Outputs:
	i) Text file of the MSMC2 analysis results.
	ii) Plot of the MSMC2 analysis results displaying effective population size estimates over time.
	iii) List of chromosomes used in the analysis.
	

## Installation

1. **Clone this repository:**

   ```bash
   git clone https://github.com/rpontius/stairway.git

## How to Run the Workflow

1. **Update input file ('samples.csv') with required information:**

  - 'organism_name' = Identifier for the organism of interest.
  - 'primary_ref' = NCBI Primary reference genome ID.
  - 'alternative_ref' = NCBI Alternative reference genome ID.
  - 'mutation_rate' = The average rate at which new mutations occur in the organism's genome per generation.
  - 'generations' = The average time (in years) it takes for one generation of the organism to occur.
