//! tests for a regression where airs would not correctly count the number of ancestors during the viterbi algorithm.
//! It used the amount of active ancestor at any given site, but the Markov Chain requires it to use the full amount
//! of ancestors theoretically available, even if the ancestors state is unknown at that site.

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::TreeSequenceGenerator;

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

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let len = SequencePosition::from_usize(9);
    let ancestors = ag.generate_ancestors(len);
    let ancestor_matcher = TreeSequenceGenerator::new(
        ancestors,
        len,
        1e-2,
        1e-20,
        SequencePosition::from_vec(vec![1, 2, 3, 4, 5, 6, 7, 8]),
    );
    let ts = ancestor_matcher.generate_tree_sequence().nodes;

    // when the ancestor count is incorrect, the algorithm will recombine the last ancestor too early, which will
    // increase the number of trees in the sequence

    assert_eq!(ts[7].edges()[0].end, 4);
}
