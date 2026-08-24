#!/usr/bin/env python3
"""Convert raw SINGER ARG output to tskit .trees, patched for strict tskit.

SINGER's bundled `convert_to_tskit` builds the mutation table without sorting it
or computing mutation parents, so `tables.tree_sequence()` fails on newer tskit
with TSK_ERR_MUTATION_PARENT_AFTER_CHILD. This version adds
sort() + build_index() + compute_mutation_parents() before materialising the
tree sequence. ARG-building logic (read_ts / read_mutation) is otherwise copied
verbatim from SINGER so the conversion is identical.

Operates only on existing raw SINGER files (nodes/branches/muts) -- it does NOT
run singer_master, so no resampling is needed.
"""
import sys
import argparse
import numpy as np
import tskit


def read_ts(node_file, edge_file):
    node_time = np.loadtxt(node_file)
    edge_span = np.loadtxt(edge_file)
    edge_span = edge_span[edge_span[:, 2] >= 0, :]
    length = max(edge_span[:, 1])
    tables = tskit.TableCollection(sequence_length=length)
    node_table = tables.nodes
    edge_table = tables.edges
    prev_time = -1
    for t in node_time:
        if (t == 0):
            node_table.add_row(flags=tskit.NODE_IS_SAMPLE)
        else:
            t = max(prev_time + 1e-4, t)
            node_table.add_row(time=t)
            prev_time = t
    parent_indices = np.array(edge_span[:, 2], dtype=np.int32)
    child_indices = np.array(edge_span[:, 3], dtype=np.int32)
    edge_table.set_columns(left=edge_span[:, 0], right=edge_span[:, 1],
                           parent=parent_indices, child=child_indices)
    tables.sort()
    return tables


def read_mutation(tables, mutation_file):
    mutations = np.loadtxt(mutation_file)
    n = mutations.shape[0]
    mut_pos = 0
    for i in range(n):
        if mutations[i, 0] != mut_pos:
            tables.sites.add_row(position=mutations[i, 0], ancestral_state='0')
            mut_pos = mutations[i, 0]
        site_id = tables.sites.num_rows - 1
        tables.mutations.add_row(site=site_id, node=int(mutations[i, 1]),
                                 derived_state=str(int(mutations[i, 3])))
    return


def read_ARG(node_file, branch_file, mutation_file):
    tables = read_ts(node_file, branch_file)
    read_mutation(tables, mutation_file)
    # --- patch: make the mutation table valid for strict tskit ---
    tables.sort()
    tables.build_index()
    tables.compute_mutation_parents()
    ts = tables.tree_sequence()
    return ts


def write_trees(input_prefix, output_prefix, start, end, step, fast=False):
    tag = "_fast" if fast else ""
    for i in range(start, end, step):
        trees_file = f"{output_prefix}_{i}.trees"
        node_file = f"{input_prefix}{tag}_nodes_{i}.txt"
        branch_file = f"{input_prefix}{tag}_branches_{i}.txt"
        mutation_file = f"{input_prefix}{tag}_muts_{i}.txt"
        ts = read_ARG(node_file, branch_file, mutation_file)
        ts.dump(trees_file)


def main():
    p = argparse.ArgumentParser(description='Convert SINGER ARG to tskit (patched)')
    p.add_argument('-input', required=True, help='Prefix of raw ARG files.')
    p.add_argument('-output', required=True, help='Prefix of output .trees files.')
    p.add_argument('-start', type=int, required=True)
    p.add_argument('-end', type=int, required=True)
    p.add_argument('-step', type=int, default=1)
    p.add_argument('-fast', action='store_true')
    if len(sys.argv) == 1:
        p.print_help(sys.stderr)
        sys.exit(1)
    a = p.parse_args()
    write_trees(a.input, a.output, a.start, a.end, a.step, a.fast)


if __name__ == '__main__':
    main()
