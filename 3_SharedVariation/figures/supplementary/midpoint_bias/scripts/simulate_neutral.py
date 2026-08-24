import msprime
import numpy as np
import pickle

goal_n = snakemake.params.goal_n
outfile = snakemake.output[0]
Ne = snakemake.params.human_Ne
frequency_threshold = snakemake.params.frequency_threshold
sequence_length = snakemake.params.sequence_length

mu = 1e-8  # mutation rate per site per generation

mutation_model = msprime.Mutation

# Function to run a single simulation, returning tuples of actual mutation age and midpoint age of the branch on which the mutation occurred
def run_sim(Ne, mu, sequence_length, frequency_threshold):
    ts = msprime.sim_ancestry(
        samples=2000,
        population_size=Ne,
        sequence_length=sequence_length
    )
    # Add mutations -- simple model with same mutation rate across time
    mts = msprime.sim_mutations(ts, rate=mu)
    # For each mutation record the age and the midpoint age of the branch on which it occurred.
    # Parent lookup is tree-specific, so we iterate trees then sites/mutations within each tree.
    mutation_ages = []
    for tree in mts.trees():
        for site in tree.sites():
            for mut in site.mutations:
                node = mut.node
                # Skip mutation if derived allele frequency is below the threshold
                freq = tree.num_samples(node) / mts.num_samples
                if freq < frequency_threshold:
                    continue
                parent = tree.parent(node)
                if parent == msprime.NULL:
                    continue  # skip root mutations (no parent branch)
                node_age = mts.node(node).time
                parent_age = mts.node(parent).time
                midpoint_age = (node_age + parent_age) / 2
                # Round ages to nearest generation (since msprime times are in generations)
                midpoint_age = round(midpoint_age)
                actual_age = round(mut.time)
                # mut.time is the actual time the mutation occurred on the branch
                mutation_ages.append((actual_age, midpoint_age))
    return mutation_ages
    
# Run until we have at least goal_n mutations
all_mutation_ages = []
while len(all_mutation_ages) < goal_n:
    mutation_ages = run_sim(Ne, mu, sequence_length, frequency_threshold)
    all_mutation_ages.extend(mutation_ages)
    
# Save the results to a file
with open(outfile, 'wb') as f:
    pickle.dump(all_mutation_ages, f)