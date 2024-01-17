//! tests whether airs mimics tsinfer's behavior of always choosing the oldest ancestor when multiple ancestors are equally
//! likely to recombine. This test will likely get obsolete in the future, when airs intentionally diverges from tsinfer's behavior.

use libairs::ancestors::AncestorGenerator;
use libairs::dna::VariantSite;
use libairs::ts::TreeSequenceGenerator;

#[test]
fn test_incomplete_nodes() {
    let sites = [
        [0, 0, 1, 0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1, 0, 1, 1],
        [1, 1, 0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let ancestors = ag.generate_ancestors();
    let ancestor_matcher =
        TreeSequenceGenerator::new(ancestors, 6, 1e-2, 1e-20, vec![1, 2, 3, 4, 5]);
    let ts = ancestor_matcher.generate_tree_sequence().0;

    // if tsinfer behavior is mimicked, nodes 1 and 2 are connected to root (0), and nodes 3 and 4 are connected to each of
    // them (1, 2) in two different trees. If built incorrectly, one of nodes (3, 4) will connect to the other (4, 3) one respectively

    assert_eq!(ts[1].node_intervals.len(), 1);
    assert_eq!(ts[1].node_intervals[0].parent, 0);
    assert_eq!(ts[2].node_intervals.len(), 1);
    assert_eq!(ts[2].node_intervals[0].parent, 0);

    assert_eq!(ts[3].node_intervals.len(), 2);
    assert!(ts[3].node_intervals.iter().any(|ni| ni.parent == 1));
    assert!(ts[3].node_intervals.iter().any(|ni| ni.parent == 2));
    assert_eq!(ts[4].node_intervals.len(), 2);
    assert!(ts[4].node_intervals.iter().any(|ni| ni.parent == 1));
    assert!(ts[4].node_intervals.iter().any(|ni| ni.parent == 2));
}
