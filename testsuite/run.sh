if [ ! -d "testdata" ]; then
    echo "Please run this script from the repository root directory"
    exit 1
fi

run_test() {
  target_dir="simulation-$1"
  target_file="$target_dir/sim$1.vcf"
  echo "Running test $target_dir"

  # delete old test data
  if [ -d "$target_dir" ]; then
      rm -rf "$target_dir"
  fi

  py ../testsuite/generate_tests.py $1 -p $2 -s $3 -i $4 && \
    cargo run --release --example match_ancestors -- $target_file
}

# setup cargo
export RUSTFLAGS="-C target-cpu=native"

pushd testdata > /dev/null

run_test 10000 10000 400000 4
run_test 20000 10000 800000 16

popd > /dev/null
