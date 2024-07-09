import sys
import tsinfer
import tskit
import argparse
from ete3 import Tree

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if len(sys.argv) < 2:
    print("Usage: python evaluate_tests.py <seed> <sequence_length>")
    sys.exit()

parser = argparse.ArgumentParser(
    prog='generate_tests.py',
    description='Generate test cases for a tskit / libairs test suite')

parser.add_argument('seed', type=int, help='test case seed')
parser.add_argument('sequence_length', type=int, help='sequence length')

args = parser.parse_args()

seed = args.seed
sequence_length = args.sequence_length
simulation_dir = f"simulation-{seed}"
simulation_filename = f"{simulation_dir}/sim{seed}"

tsinfer_ts = tskit.load(f"{simulation_filename }.ancestors.trees")

nodes = open(f"{simulation_dir}/nodes.tsv")
edges = open(f"{simulation_dir}/edges.tsv")
mutations = open(f"{simulation_dir}/mutations.tsv")
sites = open(f"{simulation_dir}/sites.tsv")
airs_ts = tskit.load_text(nodes, edges, mutations=mutations, sites=sites, sequence_length=sequence_length)

tree_sequence_length = airs_ts.num_trees

print("airs generated tree sequence with", tree_sequence_length, "trees")

# Check the ancestors are equal
tsinfer_ancestors = tsinfer.load(f"{simulation_filename}.ancestors")
airs_ancestor_file = open(f"{simulation_dir}/ancestors.tsv", 'r')
lines = [line.rstrip() for line in airs_ancestor_file]

# Check that the number of ancestors is the same (minus the header for lines, and minus the tsinfer virtual root ancestor)
if len(lines) - 1 != tsinfer_ancestors.num_ancestors - 1:
    eprint("airs and tsinfer disagree on number of ancestors: (airs:", len(lines) - 1, ", tsi: ",
           tsinfer_ancestors.num_ancestors - 1, ")")
    sys.exit(1)

# Check that the ancestors are the same
for line in lines[1:]:
    # extract ancestor definition
    start, end, age, focal_sites, genotypes = line.split('\t')
    # convert focal sites in int array
    focal_sites = list(map(lambda s: int(s), filter(lambda s: s != "", focal_sites[1:-1].split(', '))))

    # convert genotypes in int array
    genotypes = list(map(lambda s: int(s), list(genotypes)))

    # find ancestor in tsinfer ancestors using focal site
    ts_anc = next(a for a in tsinfer_ancestors.ancestors() if
                  len(a.focal_sites) == len(focal_sites) and all(x == y for x, y in zip(a.focal_sites, focal_sites)))

    # check if ancestors are equal
    if len(genotypes) == len(ts_anc.haplotype) and all(x == y for x, y in zip(genotypes, ts_anc.haplotype)):
        continue
    else:
        eprint("airs and tsinfer disagree on ancestor: (airs: ", focal_sites, ": ", genotypes, ", tsi: ", ts_anc, ")")
        sys.exit(1)


# Check that the tree sequences have the same number of trees
if tree_sequence_length != tsinfer_ts.num_trees:
    eprint("airs and tsinfer disagree on number of trees: (airs:", tree_sequence_length, ", tsi: ",
           tsinfer_ts.num_trees, ")")
    sys.exit(1)

# Check that the tree sequences have the same number of nodes (ancestors)
if airs_ts.num_nodes != tsinfer_ts.num_nodes:
    eprint("airs and tsinfer disagree on number of nodes: (airs:", airs_ts.num_nodes, ", tsi: ", tsinfer_ts.num_nodes,
           ")")
    sys.exit(1)

tsinfer_trees = tsinfer_ts.trees(sample_lists=True)

# Check that the tree sequences have the same trees by comparing the kc_distance of corresponding trees
rfs = list()
for airs_tree in airs_ts.trees(sample_lists=True):
    tsinfer_tree = next(tsinfer_trees)

    if airs_tree.interval != tsinfer_tree.interval:
        eprint(
            f"airs and tsinfer disagree on tree {airs_tree.index}: intervals (airs: {airs_tree.interval}, tsi: {tsinfer_tree.interval})")
        sys.exit(1)

    # rename sample nodes (leafs) from 0 to num_samples-1
    num_samples = airs_ts.num_samples
    airs_samples_dict = dict()
    tsinfer_samples_dict = dict()
    for i, node in enumerate(airs_tree.nodes()):
        if node >= num_samples:
            airs_samples_dict[node] = node - num_samples
        else:
            airs_samples_dict[node] = "A" + str(node)

    for i, node in enumerate(tsinfer_tree.nodes()):
        if node >= num_samples:
            tsinfer_samples_dict[node] = node - num_samples
        else:
            tsinfer_samples_dict[node] = "T" + str(node)

    airs_newick = [airs_tree.as_newick(root=root, node_labels=airs_samples_dict, include_branch_lengths=False) for root in airs_tree.roots]
    tsinfer_newick = [tsinfer_tree.as_newick(root=root, node_labels=tsinfer_samples_dict, include_branch_lengths=False) for root in tsinfer_tree.roots]
    for (aT, tT) in zip(airs_newick, tsinfer_newick):
        aT = Tree(aT, format=8)
        tT = Tree(tT, format=8)
        rf, max_rf, _, _, _, _, _ = aT.robinson_foulds(tT)
        rfs.append(rf / max_rf)

if len(rfs) > 0:
    mean_rf = sum(rfs) / len(rfs)
    if mean_rf > 0:
        eprint(f"airs and tsinfer disagree on tree sequences: mean robinson foulds distance: {mean_rf}")
        sys.exit(1)