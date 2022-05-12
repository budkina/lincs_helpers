import argparse
import os.path
import pkg_resources
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import multiprocessing as mp
import random
from random import sample

def save_counts_parallel(data_folder, cores_num, comp_insts_filtered, control_insts_filtered, sample_num, counts, output_folder):
    pool = mp.Pool(cores_num)
    jobs = []
    list_split = list(split(list(comp_insts_filtered.keys()), cores_num))
    for comp_list in list_split:
        comp_insts_chunk = {}
        control_insts_chunk = {}
        for comp_id in comp_list:
            comp_insts_chunk[comp_id] = comp_insts_filtered[comp_id]
            control_insts_chunk[comp_id] = control_insts_filtered[comp_id]
        job = pool.apply_async(save_counts, (counts, comp_insts_chunk, control_insts_chunk, sample_num, output_folder,))
        jobs.append(job)
    for job in jobs: 
        job.get()

def count_samples_per_time(insts):
    insts_num = {}
    for (inst_id, pert_time) in insts:
        if pert_time not in insts_num:
            insts_num[pert_time] = 0
        insts_num[pert_time] +=1
    return insts_num

def save_counts(counts, comp_insts, control_insts, sample_num, output_folder):
    comp_id_list = []
    for comp_id in comp_insts:
        # get compounds data
        comp_samples = comp_insts[comp_id]
        # count compound sample num for each time:
        comp_sample_num = count_samples_per_time(comp_samples)
        # get random sample for controls of the same length or less
        control_samples = control_insts[comp_id]
        control_sample_num = count_samples_per_time(control_samples)
        control_sample_per_time = {}
        for (inst_id, pert_time) in control_samples:
            if pert_time not in control_sample_per_time:
                control_sample_per_time[pert_time] = []
            control_sample_per_time[pert_time].append((inst_id, pert_time))
        random_samples = []
        for i in range(sample_num):
            random_sample = set()
            for pert_time in control_sample_per_time:
                if control_sample_num[pert_time] >= comp_sample_num[pert_time]:
                    sample = random.sample(control_sample_per_time[pert_time], comp_sample_num[pert_time])
                    random_sample.update(sample)
                else:
                    random_sample.update(control_sample_per_time[pert_time])
            random_samples.append(list(random_sample))
        # write data
        counts_df = pd.DataFrame()
        with open(os.path.join(output_folder,f"{comp_id}"),'w') as f:
            for i,sample in enumerate(random_samples):
                # write compound to design matrix
                for inst_time in comp_samples:
                    inst = inst_time[0]
                    counts_df[inst] = counts[inst]
                    inst = '.'.join(inst.split(':'))
                    time = inst_time[1]
                    f.write(f"{inst}\tcp\t{time}\t{i}\n")
                # write control to design matrix
                for inst_time in sample:
                    inst = inst_time[0]
                    counts_df[inst] = counts[inst]
                    inst = '.'.join(inst.split(':'))
                    time = inst_time[1]
                    f.write(f"{inst}\tcontrol\t{time}\t{i}\n")
        # write counts
        counts_df['index'] = counts.index
        counts_df = counts_df.set_index('index')
        counts_df.to_csv(os.path.join(output_folder, f"{comp_id}_counts.csv"), sep = '\t', header = True)
        comp_id_list.append(comp_id)
    comp_id_df = pd.DataFrame({'comp_ids':comp_id_list})
    comp_id_df.to_csv(os.path.join(output_folder, f"comp_ids.csv"), sep = '\t', header = True)


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

def get_compounds_with_dmso(controls, compounds):
    compounds_with_dmso = pd.DataFrame(columns = compounds.columns)
    used_rna_plates = set()
    for inst_id, inst in controls.iterrows():
        rna_plate  = inst['rna_plate']
        if rna_plate in used_rna_plates:
            continue
        compounds_with_dmso = compounds_with_dmso.append(compounds[compounds['rna_plate'] == rna_plate])
        used_rna_plates.add(rna_plate)
    return compounds_with_dmso

def get_compounds_times(controls, compounds, counts):
    control_insts = {}
    comp_insts = {}
    # get perturbation times and controls
    for inst_id, inst in compounds.iterrows():
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
        if comp_id not in control_insts:
            control_insts[comp_id] = set()
        # add compound instances
        comp_insts[comp_id].add((inst_id, pert_time))
        dmso_same_rna_plate = controls[controls['rna_plate'] == rna_plate]
        contol_inst_id = list(dmso_same_rna_plate.index)
        control_pert_time = dmso_same_rna_plate['pert_time']
        control_insts[comp_id].update(list(zip(contol_inst_id, control_pert_time)))
    return control_insts, comp_insts

def filter_compounds(control_insts, comp_insts):
    control_insts_filtered = {}
    comp_insts_filtered = {}
    for comp_id in comp_insts:
        comp_samples = comp_insts[comp_id]
        comp_times = set(list(zip(*comp_samples))[1])
        control_samples = control_insts[comp_id]
        control_times = set(list(zip(*control_samples))[1])
        if len(comp_times) > 1 and len(control_times) > 1:
            comp_insts_filtered[comp_id] = comp_samples
            control_insts_filtered[comp_id] = control_samples
    return control_insts_filtered, comp_insts_filtered

parser = argparse.ArgumentParser()
parser.add_argument('--data', help='LINCS data directory')
parser.add_argument('--output', help='Output directory')
parser.add_argument('--cores_num', type = int, help='Number of processes')
parser.add_argument('--sample_num', type = int, help='Number of random samples for selecting conrols for each compound')

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

controls = inst_info[inst_info['pert_id'] == 'DMSO']
compounds = inst_info[inst_info['pert_type'] == 'trt_cp']

# get compounds with DMSO on the same rna plate
compounds_with_dmso = get_compounds_with_dmso(controls, compounds)

print("Reading GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx")
level2_epsilon = parse(os.path.join(args.data, "GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx"))
counts = pd.DataFrame(data=level2_epsilon.data_df,
    index=level2_epsilon.row_metadata_df.index,
    columns=level2_epsilon.col_metadata_df.index)

control_insts, comp_insts = get_compounds_times(controls, compounds_with_dmso, counts)
print(f"Number of compound instances: {len(comp_insts)}")

# filter compounds with more than 2 time conditions
control_insts_filtered, comp_insts_filtered = filter_compounds(control_insts, comp_insts)

save_counts_parallel(args.data, args.cores_num, comp_insts_filtered, control_insts_filtered, args.sample_num, counts, args.output)
