//! tests whether airs mimics tsinfer's behavior of always choosing the oldest ancestor when multiple ancestors are equally
//! likely to recombine. This test will likely get obsolete in the future, when airs intentionally diverges from tsinfer's behavior.

mod common;

#[test]
fn test_incomplete_nodes() {
    let sites = [
        [0, 0, 1, 0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1, 0, 1, 1],
        [1, 1, 0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1],
        [1, 1, 1, 1, 1, 0, 1, 1],
    ];

    let ag = common::create_ancestor_generator(6, &sites);
    let ancestors = common::generate_ancestors(ag);
    let ts = common::match_ancestors(ancestors);

    // if tsinfer behavior is mimicked, nodes 1 and 2 are connected to root (0), and nodes 3 and 4 are connected to each of
    // them (1, 2) in two different trees. If built incorrectly, one of nodes (3, 4) will connect to the other (4, 3) one respectively

    assert_eq!(ts[1].edges().len(), 1);
    assert_eq!(ts[1].edges()[0].parent, 0);
    assert_eq!(ts[2].edges().len(), 1);
    assert_eq!(ts[2].edges()[0].parent, 0);

    assert_eq!(ts[3].edges().len(), 2);
    assert!(ts[3].edges().iter().any(|ni| ni.parent == 1));
    assert!(ts[3].edges().iter().any(|ni| ni.parent == 2));
    assert_eq!(ts[4].edges().len(), 2);
    assert!(ts[4].edges().iter().any(|ni| ni.parent == 1));
    assert!(ts[4].edges().iter().any(|ni| ni.parent == 2));
}
