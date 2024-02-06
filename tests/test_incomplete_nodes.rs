//! tests a regression in which the generator would end an edge in the tree sequence too early, because it considered
//! the right edge end position to be inclusive, while it is exclusive.

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::ViterbiMatcher;
use std::ops::Deref;

#[test]
fn test_incomplete_nodes() {
    let sites = [
        [0, 0, 1, 0, 0, 1],
        [1, 1, 0, 1, 1, 1],
        [1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 0],
        [1, 1, 1, 1, 1, 0],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let len = SequencePosition::from_usize(6);
    let ancestors = ag.generate_ancestors(len);

    assert_eq!(ancestors.len(), 4);
    assert_eq!(ancestors.deref()[3].haplotype(), &vec![1, 1, 1, 1]);

    let mut ancestor_matcher = ViterbiMatcher::new(ancestors, 1e-2, 1e-20);
    ancestor_matcher.match_ancestors();
    let ts = ancestor_matcher.get_tree_sequence().nodes;

    assert_eq!(ts.len(), 4);
    assert_eq!(ts[3].edges().len(), 2);
    assert_eq!(ts[3].edges()[0].start, 0);
    assert_eq!(ts[3].edges()[0].end, ts[3].edges()[1].start);
    assert_eq!(ts[3].edges()[1].end, 5); // test that the end is placed correctly, i.e. exclusive and less than the sequence length
}
