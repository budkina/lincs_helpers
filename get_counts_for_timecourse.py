import argparse
import os.path
import pkg_resources
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import multiprocessing as mp
import random
from random import sample

def save_counts(counts, comp_insts_filtered, control_insts_filtered, sample_num):
    for comp_id in comp_insts_filtered:
        counts_df = pd.DataFrame()
        with open(os.path.join(args.output,comp_id),'w') as f:
            samples = comp_insts_filtered[comp_id]
            comp_sample_num = len(samples)
            for inst_time in samples:
                inst = inst_time[0]
                time = inst_time[1]
                f.write(f"{inst}\tcp\t{time}\n")
                counts_df[inst] = counts[inst]
            samples = control_insts_filtered[comp_id]
            for i in range(sample_num):
                if comp_sample_num < len(samples):
                    samples = random.sample(samples, comp_sample_num)
            for inst_time in samples:
                inst = inst_time[0]
                time = inst_time[1]
                f.write(f"{inst}\tcontrol\t{time}\n")           
                counts_df[inst] = counts[inst]
        counts_df['index'] = counts.index
        counts_df = counts_df.set_index('index')
        counts_df.to_csv(os.path.join(args.output, f"{comp_id}_counts.csv"), sep = '\t', header = True)

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

parser = argparse.ArgumentParser()
parser.add_argument('--data', help='LINCS data directory')
parser.add_argument('--output', help='Output directory')
parser.add_argument('--cores_num', type = int, help='Number of processes')
parser.add_argument('--all_dmso', help='Use all DMSO from RNA plate, default disabled', action='store_true')

args = parser.parse_args()

# read metadata for each experiment 
inst_info = pd.read_csv(os.path.join(args.data, "GSE92742_Broad_LINCS_inst_info.txt"), sep="\t", dtype={'inst_id': str,
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


print("Reading GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx")
level2_epsilon = parse(os.path.join(args.data, "GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx"))
counts = pd.DataFrame(data=level2_epsilon.data_df,
    index=level2_epsilon.row_metadata_df.index,
    columns=level2_epsilon.col_metadata_df.index)

controls = inst_info[inst_info['pert_id'] == 'DMSO']
compounds = inst_info[inst_info['pert_type'] == 'trt_cp']
compounds_with_dmso = pd.DataFrame(columns = compounds.columns)
used_rna_plates =set()
for inst_id, inst in controls.iterrows():
    rna_plate  = inst['rna_plate']
    if rna_plate in used_rna_plates:
        continue
    compounds_with_dmso = compounds_with_dmso.append(compounds[compounds['rna_plate'] == rna_plate])
    used_rna_plates.add(rna_plate)

comp_insts = {}
contol_insts = {}
# get perturbation times and controls
for inst_id, inst in compounds_with_dmso.iterrows():
    if inst_id not in counts.columns:
        continue
    pert_id = inst['pert_id']
    cell_id = inst['cell_id']
    pert_time = inst['pert_time']
    pert_dose = inst['pert_dose']
    rna_plate = inst['rna_plate']
    comp_id = f'{pert_id}_{cell_id}_{pert_dose}'
    if comp_id not in comp_insts:
        comp_insts[comp_id] = set()
    if comp_id not in contol_insts:
        contol_insts[comp_id] = set()
    comp_insts[comp_id].add((inst_id, pert_time))
    dmso_same_rna_plate = controls[controls['rna_plate'] == rna_plate]
    contol_inst_id = list(dmso_same_rna_plate.index)
    control_pert_time = dmso_same_rna_plate['pert_time']
    contol_insts[comp_id].update(list(zip(contol_inst_id, control_pert_time)))

# filter compounds with more than 2 time conditions
comp_insts_filtered = {}
contol_insts_filtered = {}
for comp_id in comp_insts:
    comp_samples = comp_insts[comp_id]
    comp_times = set(list(zip(*comp_samples))[1])
    control_samples = contol_insts[comp_id]
    control_times = set(list(zip(*control_samples))[1])
    if len(comp_times) > 1 and len(control_times) > 1:
        comp_insts_filtered[comp_id] = comp_samples
        contol_insts_filtered[comp_id] = control_samples

pool = mp.Pool(args.cores_num)
jobs = []
list_split = list(split(list(comp_insts_filtered.keys()), args.cores_num))
for comp_list in list_split:
    comp_insts_chunk = {}
    control_insts_chunk = {}
    for comp_id in comp_list:
        comp_insts_chunk[comp_id] = comp_insts_filtered[comp_id]
        control_insts_chunk[comp_id] = contol_insts_filtered[comp_id]

    job = pool.apply_async(save_counts, (counts, comp_insts_chunk, control_insts_chunk, sample_num,))
    jobs.append(job)

for job in jobs: 
    job.get()
