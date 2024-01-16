//! This is a test for an issue where airs generates a sequence with less recombinations than
//! tsinfer

use libairs::ancestors::AncestorGenerator;
use libairs::dna::VariantSite;
use libairs::ts::TreeSequenceGenerator;

#[test]
fn test_unknown() {
    let sites = [
        [0, 1, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 1],
        [0, 1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 0],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let ancestors = ag.generate_ancestors();

    for (i, a) in ancestors.iter().enumerate() {
        println!("{}; {:?}", i + 1, a);
    }

    let ancestor_matcher = TreeSequenceGenerator::new(ancestors, 5, 1e-2, 1e-20, vec![1, 2, 3, 4]);
    let ts = ancestor_matcher.generate_tree_sequence().0;

    // there have to be three trees.
    // in the first tree, three nodes are linked to the root, and one of them has a mutation, and that one has the last
    // node as a child
    // in the second, only two nodes are linked to the root, and one of them has a mutation, and that one has still the
    // same child, but that also has the third previously root-linked node as a child
    // in the third, the same two nodes are linked to the root, but the children switch to the other one

    // identify the three nodes linked to root
    let root_linked_nodes;
    let non_root_linked_node;

    // check that exactly three nodes are linked to root throughout the sequence
    assert_eq!(
        ts[1..4]
            .iter()
            .filter(|t| t.node_intervals.iter().any(|ni| ni.parent == 0))
            .count(),
        3
    );

    if !ts[1].node_intervals.iter().any(|ni| ni.parent == 0) {
        non_root_linked_node = &ts[1];
        root_linked_nodes = [&ts[2], &ts[3], &ts[4]];
    } else if !ts[2].node_intervals.iter().any(|ni| ni.parent == 0) {
        non_root_linked_node = &ts[2];
        root_linked_nodes = [&ts[1], &ts[3], &ts[4]];
    } else if !ts[3].node_intervals.iter().any(|ni| ni.parent == 0) {
        non_root_linked_node = &ts[3];
        root_linked_nodes = [&ts[1], &ts[2], &ts[4]];
    } else {
        non_root_linked_node = &ts[4];
        root_linked_nodes = [&ts[1], &ts[2], &ts[3]];
    }

    // identify the node that moves and the two that are always linked to the root
    let always_root_linked_nodes;
    let moving_node;

    if root_linked_nodes[1].node_intervals.len() > 1 {
        moving_node = root_linked_nodes[0];
        always_root_linked_nodes = [root_linked_nodes[1], root_linked_nodes[2]];
    } else if root_linked_nodes[2].node_intervals.len() > 1 {
        moving_node = root_linked_nodes[1];
        always_root_linked_nodes = [root_linked_nodes[0], root_linked_nodes[2]];
    } else {
        moving_node = root_linked_nodes[2];
        always_root_linked_nodes = [root_linked_nodes[0], root_linked_nodes[1]];
    }

    // check which of the linked nodes has the first mutation, and set that to be the left one (to mirror how the
    // tsinfer sequence looks visually)
    let (left_root_linked, right_root_linked) =
        if always_root_linked_nodes[0].mutations[0] < always_root_linked_nodes[1].mutations[0] {
            (always_root_linked_nodes[0], always_root_linked_nodes[1])
        } else {
            (always_root_linked_nodes[1], always_root_linked_nodes[0])
        };

    // check that the non-root-linked node is connected to the right root-linked node, then to the left one
    assert_eq!(non_root_linked_node.node_intervals.len(), 2);
    assert_eq!(
        non_root_linked_node.node_intervals[0].parent,
        right_root_linked.ancestor_index
    );
    assert_eq!(
        non_root_linked_node.node_intervals[1].parent,
        left_root_linked.ancestor_index
    );

    // check that the moving node is connected to the root, then to the non-root-linked node
    assert_eq!(moving_node.node_intervals.len(), 2);
    assert_eq!(moving_node.node_intervals[0].parent, 0);
    assert_eq!(
        moving_node.node_intervals[1].parent,
        non_root_linked_node.ancestor_index
    );
}
