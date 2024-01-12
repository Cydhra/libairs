import sys
import tsinfer
import tskit
import argparse

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if len(sys.argv) < 2:
    print("Usage: python evaluate_tests.py <seed>")
    sys.exit()

parser = argparse.ArgumentParser(
    prog='generate_tests.py',
    description='Generate test cases for a tskit / libairs test suite')

parser.add_argument('seed', type=int, help='test case seed')

args = parser.parse_args()

seed = args.seed
simulation_dir = f"simulation-{seed}"
simulation_filename = f"{simulation_dir}/sim{seed}"

tsinfer_ts = tskit.load(f"{simulation_filename }.ancestors.trees")

nodes = open(f"{simulation_dir}/nodes.tsv")
edges = open(f"{simulation_dir}/edges.tsv")
airs_ts = tskit.load_text(nodes, edges)

tree_sequence_length = len(airs_ts.trees())

print("airs generated tree sequence with", tree_sequence_length, "trees")

if tree_sequence_length != len(tsinfer_ts.trees()):
    eprint("airs and tsinfer disagree on number of trees")
    sys.exit(1)

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