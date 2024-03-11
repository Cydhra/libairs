//! tests a regression in which the generator would end an edge in the tree sequence too early, because it considered
//! the right edge end position to be inclusive, while it is exclusive.

use std::ops::Deref;

mod common;

#[test]
fn test_incomplete_nodes() {
    let sites = [
        [0, 0, 1, 0, 0, 1],
        [1, 1, 0, 1, 1, 1],
        [1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 0],
    ];

    let ag = common::create_ancestor_generator(6, &sites);
    let ancestors = common::generate_ancestors(ag);

    assert_eq!(ancestors.len(), 4);
    assert_eq!(ancestors.deref()[3].haplotype(), &vec![1, 1, 1, 1]);

    let ts = common::match_ancestors(ancestors);

    assert_eq!(ts.len(), 4);
    assert_eq!(ts[3].edges().len(), 2);
    assert_eq!(ts[3].edges()[0].start, 0);
    assert_eq!(ts[3].edges()[0].end, ts[3].edges()[1].start);
    assert_eq!(ts[3].edges()[1].end, 5); // test that the end is placed correctly, i.e. exclusive and less than the sequence length
}
