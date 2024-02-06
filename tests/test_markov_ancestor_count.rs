//! tests for a regression where airs would not correctly count the number of ancestors during the viterbi algorithm.
//! It only counted ancestors from which recombination events were possible, but to obtain correct probabilies, it is necessary
//! to also count the current ancestor (i.e. count the possibility of recombining to itself, because this is how the
//! Markov Chain in tsinfer is modelled)

use libairs::ancestors::AncestorGenerator;
use libairs::dna::{SequencePosition, VariantSite};
use libairs::ts::TreeSequenceGenerator;

#[test]
fn test_markov_ancestor_count() {
    let sites = [
        [1, 1, 1, 1, 1, 0, 1, 1, 1, 1],
        [1, 1, 1, 1, 1, 1, 1, 1, 0, 1],
        [0, 1, 0, 1, 1, 1, 1, 0, 0, 1],
        [0, 1, 0, 1, 1, 1, 1, 0, 1, 1],
        [1, 0, 1, 0, 0, 0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1, 1, 1, 0, 1, 1],
    ];

    let ag = AncestorGenerator::from_iter(
        sites
            .iter()
            .enumerate()
            .map(|(i, site)| VariantSite::new(site.to_vec(), i + 1)),
    );

    let len = SequencePosition::from_usize(7);
    let ancestors = ag.generate_ancestors(len);
    let ancestor_matcher = TreeSequenceGenerator::new(
        ancestors,
        len,
        1e-2,
        1e-20,
        SequencePosition::from_vec(vec![1, 2, 3, 4, 5, 6]),
    );
    let ts = ancestor_matcher.generate_tree_sequence().0;

    // when the ancestor count is incorrect, the algorithm will recombine at the wrong spots. So we test the correct
    // layout of the tree sequence

    // the first nodes (1, 2) after the root are only connected to the root
    assert_eq!(ts[1].edges().len(), 1);
    assert_eq!(ts[1].edges()[0].parent, 0);

    assert_eq!(ts[2].edges().len(), 1);
    assert_eq!(ts[2].edges()[0].parent, 0);

    // the next node is connected to one of the previous nodes, and then to the other one
    assert_eq!(ts[3].edges().len(), 2);
    assert!(ts[3].edges().iter().any(|ni| ni.parent == 1));
    assert!(ts[3].edges().iter().any(|ni| ni.parent == 2));

    // the next two nodes build a path from the previous node downwards
    assert_eq!(ts[4].edges().len(), 1);
    assert_eq!(ts[4].edges()[0].parent, 3);

    assert_eq!(ts[5].edges().len(), 1);
    assert_eq!(ts[5].edges()[0].parent, 4);

    // the last node is the most interesting one, because it is the node that got recombined at the wrong side when
    // airs counted ancestors incorrectly. It is supposed to connect to (3), and then at site 5 it will reconnect to the
    // root. airs did that at site 6 when the error was present
    assert_eq!(ts[6].edges().len(), 2);
    assert_eq!(ts[6].edges()[0].parent, 3);
    assert_eq!(ts[6].edges()[0].end, 5); // exclusive
    assert_eq!(ts[6].edges()[1].parent, 0);
    assert_eq!(ts[6].edges()[1].start, 5);
}
