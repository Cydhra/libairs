if [ ! -d "testdata" ]; then
    echo "Please run this script from the repository root directory"
    exit 1
fi

# setup cargo
export RUSTFLAGS="-C target-cpu=native"

pushd testdata

# delete old test data
if [ -d "simulation-10000" ]; then
    rm -rf simulation-10000
fi

py ../testsuite/generate_tests.py 10000 -p 10000 -s 400000
cargo run --release --example match_ancestors -- ./simulation-10000/sim10000.vcf

popd
