import os
import pickle
import subprocess
import sys
import tsinfer


def generate_test_inputs(file, compressed, filter):
    # Generates inputs for airs and tsinfer from a VCF file such that tsinfer and airs can be guaranteed to get the
    # same input.
    args = ["cargo", "run", "--release", "--example", "convert_inputs", "--", "-i", file];

    if compressed:
        args.append("--compressed")

    for contig in filter:
        args.append(contig)

    environment = os.environ.copy()
    environment["RUSTFLAGS"] = "-C target-cpu=native"

    process = subprocess.Popen(args, env=environment)
    process.wait()

    # for each contig.variants file in the directory of file, generate a sample_data file
    input_dir = os.path.dirname(file)
    for contig in filter:
        variants = pickle.load(open(f"{input_dir}/{contig}.variants", "rb"))
        with tsinfer.SampleData(path=f"{input_dir}/{contig}.samples", sequence_length=variants["sequence_length"],
                                max_file_size=2 ** 31) as sample_data:
            for site in variants["sites"]:
                sample_data.add_site(site["position"], site["genotypes"],
                                     [site["ancestral_state"], site["derived_state"]])
            sample_data.finalise()


if __name__ == "__main__":
    vcf_file = sys.argv[1]
    filter = sys.argv[2:]

    generate_test_inputs(vcf_file, False, filter)
