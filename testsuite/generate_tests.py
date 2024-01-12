import sys
import os
import msprime
import tsinfer
import argparse

if len(sys.argv) < 4:
    print("Usage: python generate_tests.py <seed> <pop_size> <seq_length>")
    sys.exit()

parser = argparse.ArgumentParser(
    prog='generate_tests.py',
    description='Generate test cases for a tskit / libairs test suite')

parser.add_argument('seed', type=int, help='test case seed')
parser.add_argument('-p', '--pop_size', type=int, help='population size')
parser.add_argument('-s', '--seq_length', type=int, help='sequence length')
args = parser.parse_args()

seed = args.seed
pop_size = args.pop_size
seq_length = args.seq_length

os.mkdir(f"simulation-{seed}")

sweep_model = msprime.SweepGenicSelection(
    position=seq_length / 2, start_frequency=0.0001, end_frequency=0.9999, s=0.25, dt=1e-6)

ts = msprime.sim_ancestry(
    4,
    model=[sweep_model, msprime.StandardCoalescent()],
    population_size=pop_size,
    sequence_length=seq_length,
    recombination_rate=1e-8,
    random_seed=seed,
)
# Optionally add finite-site mutations to the ts using the Jukes & Cantor model, creating SNPs
ts = msprime.sim_mutations(ts, rate=1e-8, random_seed=seed)

print(ts)

output_file = open(f"simulation-{seed}/sim{seed}.vcf", "w")
ts.write_vcf(output_file)

with tsinfer.SampleData(path=f"simulation-{seed}/sim{seed}.samples", sequence_length=seq_length) as sample_data:
    for variant in ts.variants():
        sample_data.add_site(variant.position, variant.genotypes, variant.alleles)

    sample_data.finalise()
    ancestor_data = tsinfer.generate_ancestors(
        sample_data,
        path=f"simulation-{seed}/sim{seed}.ancestors"
    )
    ancestor_tree = tsinfer.match_ancestors(sample_data, ancestor_data)
    ancestor_tree.dump(f"simulation-{seed}/sim{seed}.ancestors.trees")
