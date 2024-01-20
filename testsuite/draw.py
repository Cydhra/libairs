import sys
import tskit
import argparse


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


if len(sys.argv) < 2:
    print("Usage: python draw.py <seed>")
    sys.exit()

parser = argparse.ArgumentParser(
    prog='generate_tests.py',
    description='Generate test cases for a tskit / libairs test suite')

parser.add_argument('seed', type=str, help='test case seed')
parser.add_argument("--width", type=int, default=2048)
parser.add_argument("--height", type=int, default=800)

args = parser.parse_args()

seed = args.seed
width = args.width
height = args.height

simulation_dir = f"testdata/simulation-{seed}"
simulation_filename = f"{simulation_dir}/sim{seed}"

tsinfer_ts = tskit.load(f"{simulation_filename}.ancestors.trees")

nodes = open(f"{simulation_dir}/nodes.tsv")
edges = open(f"{simulation_dir}/edges.tsv")
mutations = open(f"{simulation_dir}/mutations.tsv")
sites = open(f"{simulation_dir}/sites.tsv")
airs_ts = tskit.load_text(nodes, edges, mutations=mutations, sites=sites)

output = open(f"{simulation_dir}/tsinfer.svg", "w+")
svg_size = (width, height)
svg_string = tsinfer_ts.draw_svg(
    size=svg_size,
    y_axis=True,
    y_label=" ",
    time_scale="rank",
    x_scale="treewise",  # Force same axis settings as the text view
)
output.write(svg_string)

output = open(f"{simulation_dir}/airs.svg", "w+")
svg_size = (width, height)
svg_string = airs_ts.draw_svg(
    size=svg_size,
    y_axis=True,
    y_label=" ",
    time_scale="rank",
    x_scale="treewise",  # Force same axis settings as the text view
)
output.write(svg_string)
