#!/usr/bin/env python3

# Need to generate file where each line is a command for running input generator
import pandas as pd
import numpy as np

# Reading in genes of interest (those w/ better/worse performance w/ centromere haplotypes included
better_w_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_w_cent.tsv", sep="\t")
better_wo_df = pd.read_csv("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/weights/tables/better_wo_cent.tsv", sep="\t")

# Generating list of commands
base_command = "module load miniconda; conda activate /home/sl2749/.conda/envs/tensorqtl; python /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/tensorqtl.v2.GATK_1KGSV_Pangenie2final_EDGE_TR.bp.rescue_mas5.nonredundant.py "

# Better w/ centromeres
better_w_commands = list(better_w_df.apply(lambda row: base_command + row["chr"] + " " + row["gene"], axis=1))
# Better w/out centromeres
better_wo_commands = list(better_wo_df.apply(lambda row: base_command + row["chr"] + " " + row["gene"], axis=1))
# All commands
commands = better_w_commands + better_wo_commands

# Writing to file line by line
with open("/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/full_variant_tests/input_generator_commands.sh", "w") as file:
    for line in commands:
        file.write(f"{line}\n")
