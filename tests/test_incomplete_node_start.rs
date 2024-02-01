//! Tests a regression where airs set the initial likelihood of ancestors that did not start at the first position to
//! 1. This meant that if during traceback the incomplete ancestor was selected, it would not insert another
//! recombination event at the end of the ancestor, so the non-existing sites of the ancestor were inserted as parents
//! of younger nodes.

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::TreeSequenceGenerator;
use std::ops::Deref;

#[test]
fn test_incomplete_node_start() {
    let sites = [
        [0, 0, 0, 1, 0, 1, 1, 0],
        [1, 1, 1, 0, 1, 0, 0, 1],
        [1, 1, 1, 0, 1, 0, 0, 1],
        [0, 1, 0, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 0, 0, 0, 1],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let ancestors = ag.generate_ancestors();

    // fixme this is broken because start() publicly exposes variant site
    // assert_ne!(ancestors.deref()[2].start(), 0); // the third ancestor is incomplete and doesn't start at position 0.

    let ancestor_matcher = TreeSequenceGenerator::new(
        ancestors,
        SequencePosition::from_usize(6),
        1e-2,
        1e-20,
        SequencePosition::from_vec(vec![1, 2, 3, 4, 5]),
    );
    let ts = ancestor_matcher.generate_tree_sequence().0;

    // assert that no ancestor has the third ancestor as its parent in the first variant site, because it doesnt exist
    // for this position
    ts[3..].iter().for_each(|node| {
        assert!(
            node.node_intervals[0].parent != 2
                || node.node_intervals[0].start > SequencePosition::from_usize(0)
        );
    });
}
