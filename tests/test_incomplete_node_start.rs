//! Tests a regression where airs set the initial likelihood of ancestors that did not start at the first position to
//! 1. This meant that if during traceback the incomplete ancestor was selected, it would not insert another
//! recombination event at the end of the ancestor, so the non-existing sites of the ancestor were inserted as parents
//! of younger nodes.

use libairs::ts::ViterbiMatcher;

mod common;

#[test]
fn test_incomplete_node_start() {
    let sites = [
        [0, 0, 0, 1, 0, 1, 1, 0],
        [1, 1, 1, 0, 1, 0, 0, 1],
        [1, 1, 1, 0, 1, 0, 0, 1],
        [0, 1, 0, 1, 0, 0, 1, 1],
        [0, 1, 0, 0, 0, 0, 0, 1],
    ];

    let ag = common::create_ancestor_generator(8, &sites);
    let ancestors = ag.generate_ancestors();
    let mut ancestor_matcher = ViterbiMatcher::new(ancestors, 1e-2, 1e-20, false, 1);
    ancestor_matcher.match_ancestors();
    let ts = ancestor_matcher.get_tree_sequence().nodes;

    // assert that no ancestor has the third ancestor as its parent in the first variant site, because it doesnt exist
    // for this position
    ts[3..].iter().for_each(|node| {
        assert!(node.edges()[0].parent != 2 || node.edges()[0].start > 0);
    });
}
