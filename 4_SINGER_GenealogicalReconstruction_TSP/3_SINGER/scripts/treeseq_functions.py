#!/usr/bin/env python3
import matplotlib.pyplot as plt

# Assign population from species_dict
def assign_population(ts, species_dict):
    # Count length of species_dict
    n_haps = len(species_dict)
    # Get number of unique values (species) in species_dict
    n_species = len(set(species_dict.values()))
    # Make a dict of species and their corresponding population number
    species_idx_dict = {}
    for i, species in enumerate(set(species_dict.values())): 
        species_idx_dict[species] = i
        species_idx_dict[i] = species
    # Get tables
    new_tables = ts.dump_tables()
    # Get node table
    node_table = new_tables.nodes
    # Get population array
    populations = node_table.population
    # Get population array, replace first n_haps values with 0 (human) or 1 (chimp) according to species_dict 
    for j in range(0,n_haps):
        populations[j] = species_idx_dict[species_dict[j]]
    # Assign new populations to tables
    new_tables.nodes.population = populations
    # Edit population table
    pop_table = new_tables.populations
    # Add n_species new popluations - 0 for human and 1 for chimp
    for k in range(0,n_species):
        pop_table.add_row(metadata=b'{"name": species_idx_dict[k]}')
    # Update tables
    new_tables.populations.replace_with(pop_table)
    # Make new tree sequence
    ts_new = new_tables.tree_sequence()
    return ts_new


# Conditions:
# 1. At the time of the mutation, there are 3 or fewer total lineages in the tree 
# 2. Below there mutation, there are individuals from multiple species (populations) and the first coalescent event is between ind. from different species
# 3. Mutation is not on the root branch

# Function to determine if mutation's tree has a topology consistent with trans-species polymorphism - returns True or False dep on iif all condidtions at met:
# Conditions:
# 1. At the time of the mutation, there are 3 or fewer total lineages in the tree 
# 2. Below there mutation, there are individuals from multiple species (populations) and the first coalescent event is between ind. from different species
# 3. Mutation is not on the root branch

def tsp_topology(ts, tree, mut_id, species_dict):
    if nonroot_condition(ts, tree, mut_id) == False:
        return False
    # Check if there are 3 or fewer total lineages in the tree before the mutation
    elif (n_nodes_before_mut(ts, tree, mut_id) > 2):
        return False
    # Check if the first coalescence below the mutation is between monophyletic groups of different species
    elif (descendants_condition(ts, tree, mut_id, species_dict) == False):
        return False
    else:
        return True

# New condition - should be true in case of the 'yellow' topology defined in Gao et al. 2014 
# When there are only 3 lineages, is first coalescence is between two lineages that are from diff species, and that the 3rd lineage is also a mono-species lineage
def yellow_topology_condition(ts, tree, species_dict):
    # Get list of node ages in order
    node_ages_dict = {}
    for node in tree.nodes():
        node_ages_dict[node] = tree.time(node)
    # Sort dictionary by values (ages) in descending order
    node_ages_dict_sorted = sorted(node_ages_dict.items(), key=lambda x: x[1], reverse=True)
    # Get the ID of the second oldest node
    second_oldest_node = node_ages_dict_sorted[1][0]
    # Get the children of the second oldest node
    child_1, child_2 = tree.children(second_oldest_node)
    # Check whether the first coalescence below this node is between monophyletic groups of different species
    # Get the species of the terminal nodes under child 1
    species_child_1 = []
    for child in tree.leaves(child_1):
        species_child_1.append(species_dict[child])
    # Get unique only
    species_child_1 = list(set(species_child_1))
    # Get the species of the terminal nodes under child 2
    species_child_2 = []
    for child in tree.leaves(child_2):
        species_child_2.append(species_dict[child])
    # Get unique only
    species_child_2 = list(set(species_child_2))
    # Get the species under the 3rd lineage (other child of oldest node)
    oldest_node_children = tree.children(node_ages_dict_sorted[0][0])
    # Get the id of node that is not the second oldest 
    if oldest_node_children[0] == child_1:
        third_child = oldest_node_children[1]
    else:
        third_child = oldest_node_children[0]
    # Get the species of the terminal nodes under the third child
    species_child_3 = []
    for child in tree.leaves(third_child):
        species_child_3.append(species_dict[child])
    # Get unique only
    species_child_3 = list(set(species_child_3))
    # Check if the first coalescence is between monophyletic groups of different species
    if len(species_child_1) > 1 or len(species_child_2) > 1:
        return False
    elif species_child_1[0] == species_child_2[0]:
        return False
    else:
        return True

def nonroot_condition(ts, tree, mut_id):
    # Get mutation
    mut = ts.mutation(mut_id)
    # Check if mutation is on root branch (mutation node is first in node list)
    node = ts.node(mut.node)
    node_list = list(tree.nodes())
    if node_list[0] == node.id:
        return False
    else:
        return True

def n_nodes_before_mut(ts, tree, mut_id):
    ### Calculate the number of nodes in the tree before the mutation (easier to extract than number of lineages)
    mut = ts.mutation(mut_id)
    top = tree.time(ts.edge(mut.edge).parent)
    bottom = tree.time(ts.edge(mut.edge).child)
    age_est = (top + bottom)/2
    # Get list of node ages in order
    node_ages = []
    for node in tree.nodes():
        node_ages.append(tree.time(node))
    # Order - oldest to youngest
    node_ages.sort(reverse=True)
    # Get number of nodes older than mutation
    nodes_above_mut = len([x for x in node_ages if x > age_est])
    return nodes_above_mut

def descendants_condition(ts, tree, mut_id, species_dict):
    ### Check whether the first coalescence below the mutation is between monophyletic groups of different species
    # Get mutation
    mut = ts.mutation(mut_id)
    # Get node
    node = ts.node(mut.node)
    # Get the descendant nodes of the mutation node
    children = tree.children(node.id)
    # A mutation sitting on a leaf (or any non-bifurcating node) has no first
    # coalescence below it, so the between-species coalescence condition cannot
    # be met -> not a TSP by this criterion.
    if len(children) != 2:
        return False
    child_1, child_2 = children
    # Get the species of the terminal nodes under child 1
    species_child_1 = []
    for child in tree.leaves(child_1):
        species_child_1.append(species_dict[child])
    # Get unique only
    species_child_1 = list(set(species_child_1))
    # Get the species of the terminal nodes under child 2
    species_child_2 = []
    for child in tree.leaves(child_2):
        species_child_2.append(species_dict[child])
    # Get unique only
    species_child_2 = list(set(species_child_2))
    # Check if the first coalescence is between monophyletic groups of different species
    if len(species_child_1) > 1 or len(species_child_2) > 1:
        return False
    elif species_child_1[0] == species_child_2[0]:
        return False
    else:
        return True

def mutation_age(ts, mut_id):
    # Calculate mutation age as midpoint between parent and child nodes
    mut = ts.mutation(mut_id)
    top = ts.node(mut.node).time
    bottom = ts.node(ts.edge(mut.edge).child).time
    age_est = (top + bottom)/2
    return age_est

# Function to order mutations at a site by age and determine which are fwd vs. rev mutations
# Returns a dictionary (ordered by age) of mutation ids as keys, and value as tuple of age and state (fwd or rev)
def order_mutations(ts, pos):
    # Retrieve list of mutations at the site
    mut_list = ts.site(position = pos).mutations
    # Initialize dictionary to store mutation ages
    mut_age_dict = {}
    # For each mutation, calculate age and add to dictionary
    for mut in mut_list:
        mut_age_dict[mut.id] = [mutation_age(ts, mut.id), ts.mutation(mut.id).derived_state]
    # Order mutations by age (oldest to youngest)
    mut_age_dict = dict(sorted(mut_age_dict.items(), key = lambda x: x[1][0], reverse = True))
    oldest_mut_id = max(mut_age_dict, key = lambda x: mut_age_dict[x][0])
    oldest_state = ts.mutation(oldest_mut_id).derived_state
    # Re-assign derived states in dictionary as fwd or rev depending on if they match the state of the oldest mutation
    for mut_id in mut_age_dict.keys():
        if mut_age_dict[mut_id][1] == oldest_state:
            mut_age_dict[mut_id] = (mut_age_dict[mut_id][0], "fwd")
        else:
            mut_age_dict[mut_id] = (mut_age_dict[mut_id][0], "rev")
    return mut_age_dict

def plot_tree_svg_mut(ts, pos, pos_other, start, species_dict):
    # Get tree at pos
    tree = ts.at(pos)
    svg_size = (1000, 600)
    tree_length = tree.interval[1] - tree.interval[0]
    # Initialize CSS string
    css_string_shared = ""
    # Get mut_age_dict for focal position
    mut_dict_focal = order_mutations(ts, pos)
    # Loop over keys in mut_age_dict
    for mut_id in mut_dict_focal.keys():
        # If state of mutation is fwd --> color branches below magenta
        if mut_dict_focal[mut_id][1] == "fwd":
            css_string_shared = css_string_shared + ".mut.m" + str(mut_id) + " .sym, .m" + str(mut_id) + ">line, .m" + str(mut_id) + ">.node .edge{stroke:magenta} .mut.m" + str(mut_id) + " .lab{fill:magenta}"
        # Else revert colors to black, color the mutation itself green
        else:
            css_string_shared = css_string_shared + ".mut.m" + str(mut_id) + " .sym, .m" + str(mut_id) + ">line, .m" + str(mut_id) + ">.node .edge{stroke:black} .mut.m" + str(mut_id) + " .lab{fill:green}"
    # Check if mut_id_other entry 0 is "NA"
    if pos_other != "NA":
        for pos in pos_other:
            mut_dict_pos = order_mutations(ts, pos)
            for mut_id in mut_dict_pos.keys():
                # Color forward mutation nodes (but not edges below) magenta
                if mut_dict_pos[mut_id][1] == "fwd":
                    css_string_shared = css_string_shared + ".mut.m" + str(mut_id) + " .sym, .lab{fill:magenta}"
                # Color reverse mutation nodes (but not edges below) green
                else:
                    css_string_shared = css_string_shared + ".mut.m" + str(mut_id) + " .sym, .lab{fill:green}"
    # Title string
    # Get real position (pos + start)
    real_pos = pos + start
    title_string = "Tree at position " + str(int(real_pos)) + " - " + str(int(tree_length)) + "bp"
    # Make list of oldest 4 nodes in tree for y-axis ticks
    y_ticks = []
    for node in tree.nodes():
        y_ticks.append(tree.time(node))
    y_ticks.sort(reverse=True)
    # Get only the 4 oldest nodes
    y_ticks = y_ticks[:4]
    # Other formatting
    css_string = (
        ".node.p0 > .sym {fill: red}"
        ".node.p1 > .sym {fill: blue}"
        # Hide internal node labels & symbols
        ".node:not(.leaf) > .sym, .node:not(.leaf) > .lab {display: none}"
        # Display xlabel at the bottom, large
        ".xlabel {font-size: 20px; text-anchor: middle}"
    )
    # Add background colouring to distinguish between tree based on whether they pass or fail the TSP test
    if tsp_topology(ts, tree, max(mut_dict_focal), species_dict) == True:
        css_string = css_string + "svg {background-color: white}"
    else:
        css_string = css_string + "svg {background-color: lightgrey}"
    # Combine css strings
    css_string = css_string + css_string_shared
    # Draw tree
    svg = tree.draw_svg(
        size = svg_size,
        node_labels = {},
        style = css_string,
        time_scale = "time",
        x_label = title_string,
        y_axis = True,
        y_label = "Time (generations)",
        y_ticks = y_ticks
    )
    return svg
    
def plot_tmrca_dict(tmrca_dict, interval_dict, title, start, end):
    # Create a figure with two subplots
    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True, figsize=(10,5), gridspec_kw={'height_ratios': [5, 1]})
    # Unpack dictionary and convert intervals to be plotted with stairs function 
    tmrcas = [0,1,2] # Human, Chimp, All
    tmrca_labels = ['Human', 'Chimp', 'All']
    for tmrca in tmrcas:
        x = []
        y= []
        for interval in tmrca_dict.keys():
            x.append(interval[0])
            y.append(tmrca_dict[interval][tmrca])
        ax1.plot(x, y, label=str(tmrca_labels[tmrca])) 
    # Set labels and title
    ax1.set_xlabel('Genomic position')
    ax1.set_xlim(start, end)
    ax1.set_ylabel('TMRCA (generations)')
    ax1.set_title(title)
    ax1.legend()
    # Plot segment topology states
    ax2.set_ylabel('Topology state')
    # No y ticks
    ax2.set_yticks([])
    # Set x limits
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_ylim(0,10)
    # Plot the topology states
    for interval, state in interval_dict.items():
        ax2.fill_betweenx([0, 10], interval[0], interval[1], color=state, edgecolor='None')
    return fig, ax1, ax2

def plot_states_main(master_dict, start, end, tsp):
    # Get keys and make int list
    keys = list(master_dict.keys())
    min_keys = min(keys)
    max_keys = max(keys)
    # Get number of samples
    num_samples = len(master_dict)
    # Create a figure and axis
    fig, ax = plt.subplots(figsize=(10, 5))
    # Set axis labels
    ax.set_xlabel("Position")
    ax.set_ylabel("Tree sequence sample")
    # Set title
    ax.set_title("Topology states along the genome")
    # Set x-axis limits
    ax.set_xlim(start, end)
    # Set y-axis limits
    ax.set_ylim(min_keys-5, max_keys+15)
    # Hide y-ticks
    ax.set_yticks([])  
    # Plot each sample
    for i in master_dict.keys():
        # Value is interval dictionary
        intervals = master_dict[i]
        # Plot the intervals
        plot_sample(intervals, y_axis=int(i), height=5)  # Adjust y_axis and height as needed
    for pos in tsp:
        # Plot each position as an X, above sample bars
        plt.plot(pos, max_keys+10, marker=11, color='black', markersize=8)
    # Return the figure and axis
    return fig, ax


# Function to plot a line for a single sample of the tree sequence, witb each segment colored according to its topology state
def plot_sample(interval_dict, y_axis, height):
    # For each interval, fill between interval start and end with the corresponding color
    for i in interval_dict:
        plt.fill_betweenx([y_axis, y_axis + height], i[0], i[1], color=interval_dict[i], edgecolor='None')
