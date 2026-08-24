#!/usr/bin/env python3
"""Figure 5A: the local genealogy at two of the MUC7 TSPs.

Draws one posterior local tree with the trans-species clade picked out, ready to place in the
figure layout: the y-axis in years rather than generations, mutation labels in genomic rather
than region-local coordinates, and the two legends.

The reported panel uses seed 3, sample 134.

Usage:
    python plot_panel_a.py --seed 3 --sample 134 --polarity epo --out fig5A.svg
    python plot_panel_a.py --seed 3 --sample 134 --fs 1.3        # larger text
"""
import argparse
import os
import re
import sys

import tskit

HERE = os.path.dirname(os.path.abspath(__file__))

# ---- inputs (edit) -------------------------------------------------------------------------
ANALYSES = "resources"
ARG_RUNS = {s: f"{ANALYSES}/singer_seed{s}" for s in (1, 2, 3)}
TREES = [os.path.join(HERE, "work", "seed{seed}", "tskit",
                      "4_70860866_71636175_{sample}.trees"),
         f"{ANALYSES}/figure5_muc7/trees/seed{{seed}}_sample{{sample}}.trees"]
# --------------------------------------------------------------------------------------------

REGION, START = "4_70860866_71636175", 70860866
TSPS = [71210321, 71210336]
GEN = 29

MAGENTA, GREEN, RED, BLUE = "magenta", "green", "red", "blue"
FONT = "Arial, Helvetica, 'Nimbus Sans', sans-serif"
FS_AXIS, FS_TICK, FS_MUT, FS_LEGEND = 30, 24, 22, 25

#              REF, ALT, EPO ancestral
ALLELES = {71210321: ("A", "T", "T"), 71210336: ("T", "C", "C")}


def trees_path(seed, sample):
    for pat in TREES:
        f = pat.format(seed=seed, sample=sample)
        if os.path.exists(f):
            return f
    raise FileNotFoundError(
        f"no .trees for seed {seed} sample {sample}; convert that posterior ARG sample to tskit "
        f"(stage 4, 3_SINGER rule convert_to_tskit) and put it at one of: " + ", ".join(TREES))


def species_dict(seed):
    for sub in ("merged_regions_downsampled", "merged_regions"):
        f = f"{ARG_RUNS[seed]}/{sub}/{REGION}.samples"
        if os.path.exists(f):
            hap = [f"{s.strip()}_{i}" for s in open(f) for i in (0, 1)]
            return {i: ("Human" if h.startswith(("HG", "NA")) else "Chimp")
                    for i, h in enumerate(hap)}
    raise FileNotFoundError(f"no {REGION}.samples under {ARG_RUNS[seed]}")


def set_populations(ts, sd):
    """Chimp -> p0 (red), Human -> p1 (blue). Fixed, not inferred from set() order."""
    t = ts.dump_tables()
    pop = t.nodes.population
    for i, sp in sd.items():
        pop[i] = 0 if sp == "Chimp" else 1
    t.nodes.population = pop
    while len(t.populations) < 2:
        t.populations.add_row()
    return t.tree_sequence()


def fmt_My(t_gen):
    return f"{t_gen * GEN / 1e6:.1f}M"


def build(seed, sample, size=(1250, 780), n_ticks=2, polarity="epo", fs=1.0):
    sys.path.insert(0, f"{ARG_RUNS[seed]}/scripts")
    import treeseq_functions as tf                          # noqa: E402

    ts = set_populations(tskit.load(trees_path(seed, sample)), species_dict(seed))
    tree = ts.at(TSPS[0] - START)

    focal, labels = {}, {}
    for gp in TSPS:
        mid = max(tf.order_mutations(ts, gp - START))
        m = ts.mutation(mid)
        ref, alt, epo = ALLELES[gp]
        nt = {"0": ref, "1": alt}
        if polarity == "epo":
            a, b = epo, (alt if epo == ref else ref)
        else:
            a, b = nt[ts.site(m.site).ancestral_state], nt[m.derived_state]
        focal[gp] = mid
        labels[mid] = f"chr4:{gp}  {a}→{b}"

    on_tree = [m for m in ts.mutations()
               if tree.interval[0] <= ts.site(m.site).position < tree.interval[1]]
    for m in on_tree:                       # every other mutation is drawn but unlabelled
        labels.setdefault(m.id, "")

    css = [f"svg {{font-family: {FONT}}}",
           f".node.p0 > .sym {{fill: {RED}}}",
           f".node.p1 > .sym {{fill: {BLUE}}}",
           ".node:not(.leaf) > .sym, .node:not(.leaf) > .lab {display: none}",
           f".y-axis .ticks text {{font-size: {FS_TICK * fs:.0f}px}}",
           f".mut .lab {{font-size: {FS_MUT * fs:.0f}px; font-style: italic; font-weight: bold}}",
           "svg {background-color: white}"]
    for mid in focal.values():              # the trans-species clade
        css.append(f".mut.m{mid} .sym, .m{mid}>line, .m{mid}>.node .edge {{stroke:{MAGENTA}}}")
        css.append(f".mut.m{mid} .lab {{fill:{MAGENTA}}}")
    for m in on_tree:
        if m.id not in focal.values():
            css.append(f".mut.m{m.id} .sym {{fill:{GREEN}; stroke:{GREEN}}}")

    node = ts.mutation(max(focal.values())).node
    ticks = {tree.time(tree.root): fmt_My(tree.time(tree.root)),
             tree.time(node): fmt_My(tree.time(node))}
    if n_ticks > 2:
        ticks = {t: fmt_My(t)
                 for t in sorted({tree.time(n) for n in tree.nodes()}, reverse=True)[:n_ticks]}

    # left pad from the widest tick label (right-anchored 5 px inside the axis) plus the
    # rotated title; right pad keeps the legend off the tree
    tick_w = max(len(t) for t in ticks.values()) * 0.62 * FS_TICK * fs
    pad_l, pad_b, pad_r = int(tick_w + 52 * fs), int(20 * fs), int(300 * fs)
    lo, hi = int(tree.interval[0] + START), int(tree.interval[1] + START)
    svg = tree.draw_svg(size=size, node_labels={}, mutation_labels=labels,
                        style="".join(css), time_scale="time", y_axis=True, y_label="",
                        y_ticks=ticks, canvas_size=(size[0] + pad_l + pad_r, size[1] + pad_b))
    W, H = size[0] + pad_l + pad_r, size[1] + pad_b
    svg = re.sub(r"(<svg\b[^>]*>)", rf"\1<g transform='translate({pad_l},0)'>", svg, count=1)
    svg = svg.replace("</svg>", "</g></svg>")

    yx = int(24 * fs)
    axis_labels = (
        f'<g font-family="{FONT}">'
        f'<text x="{yx}" y="{size[1] / 2:.0f}" font-size="{FS_AXIS * fs:.0f}" '
        f'text-anchor="middle" transform="rotate(-90,{yx},{size[1] / 2:.0f})">'
        f'Time (years)</text>'
        f'<text x="{pad_l + size[0] / 2:.0f}" y="{H - 6 * fs:.0f}" '
        f'font-size="{FS_AXIS * fs:.0f}" text-anchor="middle">'
        f'Sample of local genealogy at chr4: {lo:,} - {hi:,}</text></g>')

    rows = [("Mutations:", "black", True), ("x TSPs", MAGENTA, False)]
    if len(on_tree) > len(focal):
        rows.append(("x Other", GREEN, False))
    rows += [(None, None, None), ("Haplotypes:", "black", True),
             ("Chimpanzee", RED, False), ("Human", BLUE, False)]
    step, gap, y, items = 32 * fs, 52 * fs, 0.0, []
    for text, colour, header in rows:
        if text is None:
            y += gap
            continue
        weight = ' font-weight="bold"' if header or colour in (RED, BLUE) else ""
        deco = ' text-decoration="underline"' if header else ""
        items.append(f'<text y="{y:.0f}" fill="{colour}"{weight}{deco}>{text}</text>')
        y += step
    legend = (f'<g transform="translate({pad_l + size[0] + 26 * fs:.0f},{40 * fs:.0f})" '
              f'font-family="{FONT}" font-size="{FS_LEGEND * fs:.0f}">'
              + "".join(items) + "</g>")

    return svg.replace("</svg>", axis_labels + legend + "</svg>"), tree, ts, focal


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--seed", type=int, default=3)
    ap.add_argument("--sample", type=int, default=134)
    ap.add_argument("--ticks", type=int, default=2, help="y-axis tickmarks (default root + clade)")
    ap.add_argument("--polarity", choices=("epo", "singer"), default="epo")
    ap.add_argument("--fs", type=float, default=1.0, help="font/padding scale")
    ap.add_argument("--out")
    a = ap.parse_args()

    svg, tree, ts, focal = build(a.seed, a.sample, n_ticks=a.ticks,
                                 polarity=a.polarity, fs=a.fs)
    out = a.out or os.path.join(HERE, f"panelA_seed{a.seed}_sample{a.sample}.svg")
    open(out, "w").write(svg)
    print(f"  wrote {out}")
    print(f"  tree chr4:{tree.interval[0] + START:,.0f}-{tree.interval[1] + START:,.0f} "
          f"({tree.interval[1] - tree.interval[0]:.0f} bp)")
    print(f"  TMRCA {fmt_My(tree.time(tree.root))}, "
          f"clade below the mutations {fmt_My(tree.time(ts.mutation(max(focal.values())).node))}")
