//! tests a regression in which the generator would end an edge in the tree sequence too early, because it considered
//! the right edge end position to be inclusive, while it is exclusive.

use libairs::ancestors::AncestorGenerator;
use libairs::dna::VariantSite;
use libairs::ts::TreeSequenceGenerator;

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

    let ancestors = ag.generate_ancestors();

    assert_eq!(ancestors.len(), 4);
    assert_eq!(ancestors[3].haplotype(), &vec![1, 1, 1, 1]);

    let ancestor_matcher =
        TreeSequenceGenerator::new(ancestors, 6, 1e-2, 1e-20, vec![1, 2, 3, 4, 5]);
    let ts = ancestor_matcher.generate_tree_sequence().0;

    assert_eq!(ts.len(), 4);
    assert_eq!(ts[3].node_intervals.len(), 2);
    assert_eq!(ts[3].node_intervals[0].start, 0);
    assert_eq!(ts[3].node_intervals[0].end, ts[3].node_intervals[1].start);
    assert_eq!(ts[3].node_intervals[1].end, 5); // test that the end is placed correctly, i.e. exclusive and less than the sequence length
}
