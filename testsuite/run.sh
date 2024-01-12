if [ ! -d "testdata" ]; then
    printf "Please run this script from the repository root directory\n"
    exit 1
fi

# exit function for error reporting
exit_on_error() {
  printf "$1\n"
  exit 1
}

number_tests=0
failed_tests=0

# arguments: $1 = seed, $2 = population size, $3 = sequence length, $4 = number of samples
run_test() {
  target_dir="simulation-$1"
  target_file="$target_dir/sim$1.vcf"
  printf "Running test $target_dir\n"

  # delete old test data
  if [ -d "$target_dir" ]; then
      rm -rf "$target_dir"
  fi

  py ../testsuite/generate_tests.py $1 -p $2 -s $3 -i $4 && \
    RUSTFLAGS="-C target-cpu=native" cargo run --release --example match_ancestors -- $target_file 2> /dev/null && \
    py ../testsuite/evaluate_tests.py $1 $3 || ((failed_tests++))

  ((number_tests++))

  printf "\n"
}

# check if cargo is installed and python venv is available
cargo --version &> /dev/null || exit_on_error "cargo not installed"
py --version &> /dev/null || exit_on_error "no python virtual env available"

# run tests
pushd testdata > /dev/null

run_test 10000 10000 400000 4
run_test 20000 10000 800000 16

printf "Ran $number_tests tests, $failed_tests failed\n"
popd > /dev/null
