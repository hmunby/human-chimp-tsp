 #!/usr/bin/env python3

### Import libraries
import tskit
import pandas as pd
from treeseq_functions import *
import statistics

### Inputs from snakemake
input_ts = snakemake.input.tree
samples = snakemake.input.samples
candidate_snps = snakemake.input.all_shared
output_tsp_metric = snakemake.output.out
region = snakemake.wildcards.region
start_index = snakemake.params.start_index
end_index = snakemake.params.end_index

# Sept 10th edit. If oldest mutation is on the root, then record the age of the *second* oldest mutation instead.
# April 16th 2025 edits -- add params for start and end index of tree sequences for generalization

# Extract start and end of region
chrom = int(region.split("_")[0])
start = int(region.split("_")[1])
end = int(region.split("_")[2])

# Remove _100.trees 
prefix = input_ts.replace(f"_{start_index}.trees", "")

### Prepare dictionaries 

# Read samples into list
with open(samples) as f:
    sample_list = f.readlines()
    sample_list = [x.strip() for x in sample_list]

# Make a list of haplotype names (2 per sample - {sample_name}_0 and {sample_name}_1)
haplotype_names = []
for sample in sample_list:
    haplotype_names.append(sample + "_0")
    haplotype_names.append(sample + "_1")

# Make sample_name_dict - key is index (0 to 187) and value is haplotype name
sample_name_dict = {}
for i in range(0,len(haplotype_names)):
    sample_name_dict[i] = haplotype_names[i]


# Species dict - all starting with "HG" or "NA" are human
species_dict = {}
for i in range(0,len(haplotype_names)):
    if haplotype_names[i].startswith("HG"):
        species_dict[i] = "Human"
    elif haplotype_names[i].startswith("NA"):
        species_dict[i] = "Human"
    else:
        species_dict[i] = "Chimp"

# Load candidates df and subset positions within region
candidates = pd.read_csv(candidate_snps, header=None, sep="\t")
# Set col names to CHROM and POS
candidates.columns = ['CHROM', 'POS']
# Remove any CHROM = X 
candidates = candidates[candidates['CHROM'] != 'X']
# Set both columns to int
candidates['CHROM'] = candidates['CHROM'].astype(int)
candidates['POS'] = candidates['POS'].astype(int)
shared_pos = candidates[(candidates['CHROM'] == chrom) & (candidates['POS'] >= start) & (candidates['POS'] <= end)].iloc[:,1].tolist()

# Init dict
results_dict = {}
for sample in range(start_index,end_index):
    # Get ts_file (prefix + "_" + sample + ".trees")
    input_ts = prefix + "_" + str(sample) + ".trees"
    # Load tree sequence
    ts = tskit.load(input_ts)
    # Add population assignments 
    ts_new = assign_population(ts, species_dict)
    # Get shared snps that are present in ARG 
    pos_dict = {}
    for site in ts_new.sites():
        pos_dict[int(start + site.position)] = site.position
    shared_pos = [x for x in shared_pos if x in pos_dict] 
    # Initialize a nested dictionary to store results
    for pos in shared_pos:
        pos_idx = shared_pos.index(pos)
        tree = ts_new.at(pos_dict[pos])
        pos_ts = pos_dict[pos]
        mut_age_dict = order_mutations(ts_new, pos_ts)
        oldest_mut_id = list(mut_age_dict.keys())[0]
        oldest_age = mut_age_dict[oldest_mut_id][0]
        # Get number of forward mutations (derived_state == state)
        n_mut = len(mut_age_dict)
        n_mut_fwd = len([x for x in mut_age_dict.values() if x[1] == "fwd"])
        n_mut_rev = n_mut - n_mut_fwd
        # Calc. TSP critera
        coal_crit = descendants_condition(ts_new, tree, oldest_mut_id, species_dict)
        nodes_crit_n = n_nodes_before_mut(ts_new, tree, oldest_mut_id)
        if nodes_crit_n > 2: 
            nodes_crit = False
        else:
            nodes_crit = True
        nonroot_crit = nonroot_condition(ts_new, tree, oldest_mut_id)
        tsp = tsp_topology(ts_new, tree, oldest_mut_id, species_dict)
        if nonroot_crit == False:
            if len(mut_age_dict) > 1:
                # If tsp fails, if there is a check the second oldest mutation
                second_mut_id = list(mut_age_dict.keys())[1]
                tsp_second = tsp_topology(ts_new, tree, second_mut_id, species_dict)
                # Modify oldest_age to be the second oldest age
                oldest_age = mut_age_dict[second_mut_id][0]
        else:
            tsp_second = False
        # If pos not in results_dict, add it
        if pos not in results_dict:
            # Make sub-dictionary for sample
            results_dict[pos] = {}
            results_dict[pos][sample] = {"n_mut": n_mut, "n_mut_fwd": n_mut_fwd, "n_mut_rev": n_mut_rev, "oldest_age": oldest_age, "coal_crit": coal_crit, "nodes_crit_n": nodes_crit_n, "nodes_crit": nodes_crit, "nonroot_crit": nonroot_crit, "tsp": tsp, "tsp_second": tsp_second}
        else:
            results_dict[pos][sample] = {"n_mut": n_mut, "n_mut_fwd": n_mut_fwd, "n_mut_rev": n_mut_rev, "oldest_age": oldest_age, "coal_crit": coal_crit, "nodes_crit_n": nodes_crit_n, "nodes_crit": nodes_crit, "nonroot_crit": nonroot_crit, "tsp": tsp, "tsp_second": tsp_second}

# Convert results_dict to a df - make first column the position
results_df = pd.DataFrame.from_dict({(i,j): results_dict[i][j]
                            for i in results_dict.keys()
                            for j in results_dict[i].keys()},
                          orient='index')

# Make position a column
results_df['position'] = [x[0] for x in results_df.index]

# Groupby position and calculate:
# 1. Mean number of mutations
# 2. Mean number of forward mutations
# 3. Mean number of reverse mutations
# 4. Median oldest mutation age
# 5. Count of coal_crit == True
# 6. Count of nodes_crit == True
# 7. Count of nonroot_crit == True
# 8. Count of tsp == True
# 9. Count of tsp_second == True
results_df_grouped = results_df.groupby('position').agg({'n_mut': 'mean', 'n_mut_fwd': 'mean', 'n_mut_rev': 'mean', 'oldest_age': 'median', 'coal_crit': 'sum', 'nodes_crit': 'sum', 'nonroot_crit': 'sum', 'tsp': 'sum', 'tsp_second': 'sum'}).reset_index()

# Rename columns to be more descriptive of how they were aggregated
results_df_grouped.columns = ['position', 'mean_n_mut', 'mean_n_mut_fwd', 'mean_n_mut_rev', 'median_oldest_age', 'coal_crit_count', 'nodes_crit_count', 'nonroot_crit_count', 'tsp_count', 'tsp_second_count']

# Sum tsp_count and tsp_second_count to get tsp_first_or_second_count
results_df_grouped['tsp_first_or_second_count'] = results_df_grouped['tsp_count'] + results_df_grouped['tsp_second_count']

# Insert chrom column at start
results_df_grouped.insert(0, 'chrom', chrom)

# Write full results df and grouped results df to file
results_df.to_csv(output_tsp_metric, sep="\t", index=False)
results_df_grouped.to_csv(output_tsp_metric.replace(".txt", ".grouped.txt"), sep="\t", index=False)
