import argparse
import os.path
import pkg_resources
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import multiprocessing as mp

def save_counts(counts, comp_insts_filtered):
    for comp_id in comp_insts_filtered:
        counts_df = pd.DataFrame()
        with open(os.path.join(args.output,comp_id),'w') as f:
            samples = comp_insts_filtered[comp_id]
            for inst_time in samples:
                inst = inst_time[0]
                time = inst_time[1]
                f.write(f"{inst}\t{time}\n")
                counts_df[inst] = counts[inst]
            cell_id = comp_id.split('_')[1]
            for control_comp_id in contol_insts:
                if control_comp_id != cell_id:
                    continue

                for control in contol_insts[control_comp_id]:
                    f.write(f"{control}\t0\n")
                    counts_df[control] = counts[control]
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

comp_insts = {}
contol_insts = {}

# get perturbation times and controls
for inst_id, inst in inst_info.iterrows():

    if inst_id not in counts.columns:
        continue

    pert_type = inst['pert_type']

    # get only trt_cp data
    if pert_type != 'trt_cp' and pert_type != 'ctl_vehicle':
        continue

    pert_id = inst['pert_id']
    cell_id = inst['cell_id']
    if pert_type == 'trt_cp':
        pert_dose = inst['pert_dose']
        comp_id = f'{pert_id}_{cell_id}_{pert_dose}'

        if comp_id not in comp_insts:
            comp_insts[comp_id] = []
        pert_time = inst['pert_time']
        comp_insts[comp_id].append((inst_id, pert_time))

    if pert_type == 'ctl_vehicle':
        if cell_id not in contol_insts:
            contol_insts[cell_id] = []
        contol_insts[cell_id].append(inst_id)

# filter compounds with more than 2 time conditions
comp_insts_filtered = {}
for comp_id in comp_insts:
    samples = comp_insts[comp_id]
    times = set(list(zip(*samples))[1])
    if len(times) > 1:
        comp_insts_filtered[comp_id] = samples

pool = mp.Pool(args.cores_num)
jobs = []
list_split = list(split(list(comp_insts_filtered.keys()), args.cores_num))
for comp_list in list_split:
    comp_insts_chunk = {}
    for comp_id in comp_list:
        comp_insts_chunk[comp_id] = comp_insts_filtered[comp_id]

    job = pool.apply_async(save_counts, (counts, comp_insts_chunk,))
    jobs.append(job)

for job in jobs: 
    job.get()
