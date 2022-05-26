import sys
import argparse
import os.path
import pkg_resources
import pandas as pd
from cmapPy.pandasGEXpress.parse import parse
import multiprocessing as mp
import random

def get_compound_items(compounds, counts):
    comp_items = {}
    # get perturbation times
    for inst_id, inst in compounds.iterrows():
        if inst_id not in counts.columns:
            continue
        pert_id = inst['pert_id']
        cell_id = inst['cell_id']
        pert_time = inst['pert_time']
        if (pert_time != 6) and (pert_time !=24):
            continue
        pert_dose = inst['pert_dose']
        comp_item = f'{pert_id}_{cell_id}_{pert_dose}'
        if comp_item not in comp_items:
            comp_items[comp_item] = set()
        # add compound instances
        comp_items[comp_item].add((inst_id, pert_time))
    return comp_items

def count_samples_per_time(insts):
    insts_num = {}
    for (inst_id, pert_time) in insts:
        if pert_time not in insts_num:
            insts_num[pert_time] = 0
        insts_num[pert_time] +=1
    return insts_num

def filter_samples(items, samples_per_time):
    # 3 samples per time point
    items_filtered = {}
    for item in items:
        item_insts = items[item]
        # check if each timepoint has more than samples_per_time instances
        # count instances per time
        insts_num_per_time = count_samples_per_time(item_insts)
        times_to_remove = []
        for time in insts_num_per_time:
            if insts_num_per_time[time] < samples_per_time:
                times_to_remove.append(time)
        # remove times without samples_per_time number of instances
        item_insts_filtered = set()
        for (inst_id, pert_time) in item_insts:
            if pert_time not in times_to_remove:
                item_insts_filtered.add((inst_id, pert_time))
        if not item_insts_filtered:
            continue
        # check if several timepoints present
        times = set(list(zip(*item_insts_filtered))[1])
        if len(times) < 2:
            continue
        items_filtered[item] = item_insts_filtered
    return items_filtered

def get_controls(comp_items, compounds, controls, counts):
    # get all controls from the same rna plate
    control_items = {}
    for comp_item in comp_items:
        for (inst_id, pert_time) in comp_items[comp_item]:
            inst = compounds.loc[inst_id,]
            rna_plate = inst['rna_plate']
            dmso_same_rna_plate = controls[controls['rna_plate'] == rna_plate]
            contol_inst_ids = list(dmso_same_rna_plate.index)
            if comp_item not in control_items:
                control_items[comp_item] = set()
            control_pert_time = dmso_same_rna_plate['pert_time']
            control_items[comp_item].update(list(zip(contol_inst_ids, control_pert_time)))
    return control_items

def get_compounds_with_dmso(comp_items, control_items):
    """ Select compounds with DMSO samples on the same rna plates"""
    comp_items_selected = {}
    control_items_selected = {}
    for comp_item in comp_items:
        if comp_item in control_items:
            comp_item_insts = comp_items[comp_item]
            control_item_insts = control_items[comp_item]
            comp_times = set(list(zip(*comp_item_insts))[1])
            control_times = set(list(zip(*control_item_insts))[1])
            # select compounds where controls and compounds have two and more common timepoints
            if len(comp_times & control_times) >= 2:
                comp_items_selected[comp_item] = comp_items[comp_item]
                control_items_selected[comp_item] = control_items[comp_item]
    return comp_items_selected, control_items_selected

def save_counts_parallel(data_folder,
    cores_num,
    comp_items,
    control_items,
    control_select_type,
    sample_num,
    counts,
    output_folder):
    pool = mp.Pool(cores_num)
    jobs = []
    list_split = list(split(list(comp_items.keys()), cores_num))
    for comp_list in list_split:
        comp_items_chunk = {}
        control_items_chunk = {}
        for item in comp_list:
            comp_items_chunk[item] = comp_items[item]
            control_items_chunk[item] = control_items[item]
        job = pool.apply_async(save_counts, (counts, comp_items_chunk, control_items_chunk, control_select_type, sample_num, output_folder,))
        jobs.append(job)
    for job in jobs: 
        job.get()
    pool.close()
    pool.join()

def arrange_samples_per_time(insts):
    insts_per_time = {}
    for (inst_id, pert_time) in insts:
        if pert_time not in insts_per_time:
            insts_per_time[pert_time] = []
        insts_per_time[pert_time].append((inst_id, pert_time))
    return insts_per_time

def random_sample_controls(control_insts, comp_insts, sample_num):
    control_insts_per_time = arrange_samples_per_time(control_insts)
    # count compound sample num for each time:
    comp_insts_num_per_time = count_samples_per_time(comp_insts)
    random_samples = []
    for i in range(sample_num):
        random_sample = set()
        for pert_time in control_insts_per_time:
            if len(control_insts_per_time[pert_time]) >= comp_insts_num_per_time[pert_time]:
                sample = random.sample(control_insts_per_time[pert_time], comp_insts_num_per_time[pert_time])
                random_sample.update(sample)
            else:
                random_sample.update(control_insts_per_time[pert_time])
        random_samples.append(list(random_sample))
    return random_samples

def write_counts(item, control_insts, comp_insts, counts, output_folder):
    counts_df = pd.DataFrame()
    with open(os.path.join(output_folder,f"{item}"),'w') as f:
        # write compound to design matrix
        for (inst_id, pert_time) in comp_insts:
            inst_id_fixed = '.'.join(inst_id.split(':'))
            f.write(f"{inst_id_fixed}\tcp\t{pert_time}\t0\n")
            counts_df[inst_id_fixed] = counts[inst_id]
        # write control to design matrix
        for (inst_id, pert_time) in control_insts:
            inst_id_fixed = '.'.join(inst_id.split(':'))
            f.write(f"{inst_id_fixed}\tcontrol\t{pert_time}\t0\n")
            counts_df[inst_id_fixed] = counts[inst_id]
        # write counts
        counts_df['index'] = counts.index
        counts_df = counts_df.set_index('index')
        counts_df.to_csv(os.path.join(output_folder, f"{item}_counts.csv"), sep = '\t', header = True)

def write_counts_sampled(item, control_insts_per_time_sampled, comp_insts, counts, output_folder):
    counts_df = pd.DataFrame()
    with open(os.path.join(output_folder,f"{item}"),'w') as f:
        for i,control_insts in enumerate(control_insts_per_time_sampled):
            # write compound to design matrix
            for (inst_id, pert_time) in comp_insts:
                inst_id_fixed = '.'.join(inst_id.split(':'))
                f.write(f"{inst_id_fixed}\tcp\t{pert_time}\t{i}\n")
                counts_df[inst_id_fixed] = counts[inst_id]
            # write control to design matrix
            for (inst_id, pert_time) in control_insts:
                inst_id_fixed = '.'.join(inst_id.split(':'))
                f.write(f"{inst_id_fixed}\tcontrol\t{pert_time}\t{i}\n")
                counts_df[inst_id_fixed] = counts[inst_id]
        # write counts
        counts_df['index'] = counts.index
        counts_df = counts_df.set_index('index')
        counts_df.to_csv(os.path.join(output_folder, f"{item}_counts.csv"), sep = '\t', header = True)

def save_counts(counts, comp_items, control_items, control_select_type, sample_num, output_folder):
    for item in comp_items:
        comp_insts = comp_items[item]
        control_insts = control_items[item]
        # select controls based on type and write design matrices and counts
        if control_select_type == 1:
            write_counts(item, control_insts, comp_insts, counts, output_folder)
        elif control_select_type == 2:
            # get random sample for controls of the same length or less
            control_insts_per_time_sampled = random_sample_controls(control_insts, comp_insts, sample_num)
            write_counts_sampled(item, control_insts_per_time_sampled, comp_insts, counts, output_folder)
        else:
            sys.exit("Unknown control_select_type")

def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))

parser = argparse.ArgumentParser()
parser.add_argument('--data',
    help='LINCS data directory')
parser.add_argument('--output',
    help='Output directory')
parser.add_argument('--cores_num',
    type = int,
    help='Number of processes')
parser.add_argument('--control_select_type',
    type = int,
    help="""DMSO controls selection: 1 - select all DMSO from the same rna plates as compound samples, 
    2 - select the same number of DMSO samples as the number of compound samples at each timepoint from
    the same rna plates, random sampling of controls""",
    default=1)
parser.add_argument('--sample_num',
    type = int,
    help='Number of random samples for selecting conrols for each compound (type 2 control selection)')
parser.add_argument('--samples_per_time',
    type = int,
    help="Minimum number of samples per time point, default: 3",
    default=3)

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
print("Reading GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx")
level2_epsilon = parse(os.path.join(args.data, "GSE92742_Broad_LINCS_Level2_GEX_epsilon_n1269922x978.gctx"))
counts = pd.DataFrame(data=level2_epsilon.data_df,
    index=level2_epsilon.row_metadata_df.index,
    columns=level2_epsilon.col_metadata_df.index)

# item is {pert_id}_{cell_id}_{pert_dose}
comp_items = get_compound_items(compounds, counts)
print(f"Number of compound instances: {len(comp_items)}")

# filter compounds with more than 2 time conditions
comp_items_filtered = filter_samples(comp_items, args.samples_per_time)
control_items = get_controls(comp_items_filtered, compounds, controls, counts)
control_items_filtered = filter_samples(control_items, args.samples_per_time)

# get compounds with DMSO on the same rna plate
comp_items_selected, control_items_selected = get_compounds_with_dmso(comp_items_filtered, control_items_filtered)

save_counts_parallel(args.data,
    args.cores_num,
    comp_items_selected,
    control_items_selected,
    args.control_select_type,
    args.sample_num,
    counts,
    args.output)

comp_items_df = pd.DataFrame(columns = ['comp_ids', 'pert_id', 'pert_iname', 'cell_id', 'pert_dose'])

for item in comp_items_selected:
    pert_id, cell_id, pert_dose = item.split('_')
    pert_name = list(inst_info[inst_info['pert_id'] == pert_id]['pert_iname'])[0]
    comp_items_df = comp_items_df.append(pd.DataFrame([[item, pert_id, pert_iname, cell_id, pert_dose]], columns  = comp_items_df.columns))
comp_items_df.to_csv(os.path.join(args.output, f"comp_ids.csv"), sep = '\t', header = True)