//! tests for a regression where airs would not correctly count the number of ancestors during the viterbi algorithm.
//! It used the amount of active ancestor at any given site, but the Markov Chain requires it to use the full amount
//! of ancestors theoretically available, even if the ancestors state is unknown at that site.

mod common;

#[test]
fn test_markov_ancestor_count2() {
    let sites = [
        [0, 1, 1, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 0, 0, 1, 1],
        [0, 0, 1, 1, 0, 0, 1, 0],
        [1, 0, 1, 1, 1, 0, 1, 0],
        [1, 0, 1, 1, 1, 0, 1, 0],
        [1, 1, 0, 0, 1, 0, 0, 1],
        [1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 0, 0, 1, 1, 0, 1],
    ];

    let ag = common::create_ancestor_generator(9, &sites);
    let ancestors = common::generate_ancestors(ag);
    let ts = common::match_ancestors(ancestors);

    // when the ancestor count is incorrect, the algorithm will recombine the last ancestor too early, which will
    // increase the number of trees in the sequence

    assert_eq!(ts[7].edges()[0].end, 4);
}
