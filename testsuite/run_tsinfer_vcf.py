import sys
import tsinfer
import argparse
from cyvcf2 import VCF

if len(sys.argv) < 3:
    print("Usage: python run_tsinfer_vcf.py <seed> <seq_length>")
    sys.exit()

parser = argparse.ArgumentParser(
    prog='generate_tests.py',
    description='Manually run tsinfer from a VCF file')

parser.add_argument('seed', type=str, help='identifier for the test suite')
parser.add_argument('seq_length', type=int, help='sequence length', default=1e5)

args = parser.parse_args()
seed = args.seed
seq_length = args.seq_length
test_dir = f"testdata/simulation-{seed}"
vcf_file = f"{test_dir}/sim{seed}.vcf"
vcf = VCF(vcf_file)

with tsinfer.SampleData(path=f"{test_dir}/sim{seed}.samples", sequence_length=seq_length,
                        max_file_size=2 ** 20) as sample_data:
    for variant in vcf:
        genotypes = list()
        basemap = dict()
        counter = 1
        basemap[variant.REF] = 0
        for bases in variant.gt_bases:
            left, right = bases.split('|')
            if left not in basemap:
                basemap[left] = counter
                counter += 1
            if right not in basemap:
                basemap[right] = counter
                counter += 1

            genotypes.append(basemap[left])
            genotypes.append(basemap[right])

        sample_data.add_site(variant.POS, genotypes, list(basemap.keys()))

    sample_data.finalise()
    ancestor_data = tsinfer.generate_ancestors(
        sample_data,
        path=f"{test_dir}/sim{seed}.ancestors",
        max_file_size=2 ** 20,
    )

    ancestor_tree = tsinfer.match_ancestors(sample_data, ancestor_data, path_compression=False, precision=20)
    ancestor_tree.dump(f"{test_dir}/sim{seed}.ancestors.trees")
    print("tsinfer generated tree sequence with", len(ancestor_tree.trees()), "trees")
