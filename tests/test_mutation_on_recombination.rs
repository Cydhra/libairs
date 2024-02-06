//! This tests a regression where airs infers more mutations than tsinfer when mutations happened on sites of recombination
//! This happened because airs checked for a mutation on the ancestor we were recombing from,
//! but not on the ancestor we were recombing to at that site

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::ViterbiMatcher;

#[test]
fn test_mutation_on_recombination_site() {
    let sites = [[0, 1, 1, 0, 1, 1], [0, 1, 1, 1, 1, 1], [1, 1, 1, 1, 0, 1]];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let len = SequencePosition::from_usize(4);
    let ancestors = ag.generate_ancestors(len);

    let mut ancestor_matcher = ViterbiMatcher::new(ancestors, 1e-2, 1e-20);
    ancestor_matcher.match_ancestors();
    let ts = ancestor_matcher.get_tree_sequence().nodes;

    assert_eq!(ts.len(), 4); // root node + 3 ancestor nodes

    assert_eq!(ts[0].edges().len(), 1); // root node has one interval

    // two nodes connect to the root node
    assert_eq!(ts[1].edges().len(), 1);
    assert_eq!(ts[1].edges()[0].parent, 0);
    assert_eq!(ts[2].edges().len(), 1);
    assert_eq!(ts[2].edges()[0].parent, 0);

    // the third node connects to the second and the third node
    assert_eq!(ts[3].edges().len(), 2);

    let left_parent = ts[3].edges()[0].parent;
    let right_parent = ts[3].edges()[1].parent;

    assert!(left_parent == 1 || left_parent == 2);
    assert!(right_parent == 1 || right_parent == 2);
    assert_ne!(left_parent, right_parent);

    // there is a mutation between left_parent and the third node, but not between right_parent and the third node
    assert_eq!(ts[3].mutations().len(), 1);
    assert_eq!(ts[3].mutations()[0], 0); // first site is mutated
}
