import sys
import tskit
import argparse

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
airs_ts = tskit.load_text(nodes, edges, sequence_length=sequence_length)

tree_sequence_length = airs_ts.num_trees

print("airs generated tree sequence with", tree_sequence_length, "trees")

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

# Check that the tree sequences have the same trees by comparing the kc_distance of corresponding trees
for tree in range(tree_sequence_length):
    airs_tree = airs_ts.at(tree)
    tsinfer_tree = tsinfer_ts.at(tree)
    if airs_tree.interval != tsinfer_tree.interval:
        eprint(f"airs and tsinfer disagree on {tree}. tree interval (airs: {airs_tree.interval}, tsi: {tsinfer_tree.interval})")
        sys.exit(1)

    kc = airs_tree.kc_distance(tsinfer_tree)
    if kc != 0:
        eprint(f"airs and tsinfer disagree on {tree}. tree (kc_distance: {kc})")
        sys.exit(1)