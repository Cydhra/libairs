//! Tests a regression wherein airs used sequence-positions as indices into the variant-site-vector instead of the actual
//! vector indices and thus didn't realize that an ancestor should have been removed from the tableau during the viterbi
//! algorithm. The test case tests that the incomplete ancestor stops copying from any ancestor, and that the last ancestor,
//! which copies from the incomplete ancestor, stops copying from it before the incomplete ancestor ends

use std::ops::Deref;

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::TreeSequenceGenerator;

#[test]
fn test_incomplete_inner_nodes() {
    let sites = [
        [0, 0, 1, 1, 0, 1, 1, 0],
        [0, 1, 0, 0, 1, 0, 0, 1],
        [1, 1, 0, 0, 1, 0, 0, 1],
        [0, 1, 1, 1, 1, 1, 1, 1],
        [0, 1, 1, 1, 0, 0, 1, 1],
        [1, 0, 1, 1, 1, 0, 1, 0],
        [1, 0, 1, 1, 1, 0, 1, 0],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let ancestors = ag.generate_ancestors();
    assert_eq!(ancestors.deref()[5].len(), 6); // only 6 sites, instead of 7

    let ancestor_matcher = TreeSequenceGenerator::new(
        ancestors,
        SequencePosition::from_usize(8),
        1e-2,
        1e-20,
        SequencePosition::from_vec(vec![1, 2, 3, 4, 5, 6, 7]),
    );
    let ts = ancestor_matcher.generate_tree_sequence().0;

    assert_eq!(ts.len(), 7);

    // test that the incomplete ancestor copies from another one and stops at the site where it has no more state
    assert_eq!(ts[5].node_intervals.len(), 1);
    assert_eq!(ts[5].node_intervals[0].end, SequencePosition::from_usize(7)); // not equal to the sequence-length 8, because the ancestor doesnt have state for the last site

    assert_eq!(ts[6].node_intervals.len(), 2);
    assert_eq!(ts[6].node_intervals[0].parent, 5); // check that it copies from the incomplete ancestor
    assert_eq!(ts[6].node_intervals[0].end, SequencePosition::from_usize(7));

    // check that it copies from somewhere else
    assert_eq!(
        ts[6].node_intervals[1].start,
        SequencePosition::from_usize(7)
    );
    assert_eq!(ts[6].node_intervals[1].end, SequencePosition::from_usize(8));
}
