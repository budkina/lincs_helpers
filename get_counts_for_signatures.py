import argparse
import os.path
import pkg_resources
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse

parser = argparse.ArgumentParser()
parser.add_argument('--lincs_data', help='LINCS data directory')
parser.add_argument('--compounds', nargs='+', help='Сompounds to get counts, default to test: sorafenib bortezomib nelfinavir saracatinib alisertib',
        default=['sorafenib', 'bortezomib', 'nelfinavir', 'saracatinib', 'alisertib'])
parser.add_argument('--output', help='Output directory')
args = parser.parse_args()
pkg_resources.get_distribution("cmapPy").version

# read metadata for each experiment 
print("Reading GSE92742_Broad_LINCS_inst_info.txt")
inst_info = pd.read_csv(os.path.join(args.lincs_data, "GSE92742_Broad_LINCS_inst_info.txt"), sep="\t", dtype={'inst_id': str,
    'rna_plate': str,
    'rna_well': str,
    'pert_id': str,
    'pert_iname': str,
    'pert_type': str,
    'pert_dose': str,
    'pert_dose_unit': str,
    'pert_time': int,
    'pert_time_unit': str,
    'cell_id': str})
inst_info.set_index("inst_id", inplace=True)

# read metadata for each signature
print("Reading GSE92742_Broad_LINCS_sig_info.txt")
sig_info = pd.read_csv(os.path.join(args.lincs_data, "GSE92742_Broad_LINCS_sig_info.txt"), sep="\t", dtype={'sig_id': str,
    'pert_id': str,
    'pert_iname': str,
    'cell_id': str,
    'pert_dose': str,
    'pert_dose_unit': str,
    'pert_idose': str,
    'pert_time': int,
    'pert_time_unit': str,
    'pert_itime': str,
    'distil_id': str})

print("Reading GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx")
level2_epsilon = parse(os.path.join(args.lincs_data, "GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx"))
counts = pd.DataFrame(data=level2_epsilon.data_df,
    index=level2_epsilon.row_metadata_df.index,
    columns=level2_epsilon.col_metadata_df.index)

# read counts table
print("Saving counts tables")
df_info = pd.DataFrame(columns=['drug', 'replicates'])
for comp in args.compounds:
    # select compound signatures
    print(f"Reading {comp} counts")
    comp_sig_info = sig_info[sig_info["pert_iname"] == comp]
    comp_inst_ids = inst_info[inst_info["pert_iname"] == comp].index
    for index, sig in comp_sig_info.iterrows():

        # get all replicates for the signature
        drug_replicates = sig['distil_id']
        df_info = df_info.append({'drug':comp, 'replicates':drug_replicates}, ignore_index=True)
        replicates = drug_replicates.split('|')

        # get corresponding counts for each replicate
        for replicate in replicates:
            rna_plate = inst_info.loc[replicate]['rna_plate']
            plate_population = inst_info[(inst_info['rna_plate'] == rna_plate)]

            # write counts table
            counts[list(plate_population.index)].to_csv(os.path.join(args.output, f'{comp}_{replicate}'), sep ='\t')

# write replicates table
df_info.to_csv(os.path.join(args.output, "replicates_info.csv"), sep ='\t')
