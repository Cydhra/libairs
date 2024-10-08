use libairs::ancestors::{AncestorArray, AncestorGenerator};
use libairs::ts::{TreeSequenceNode, ViterbiMatcher};
use libairs::variants::VariantDataBuilder;

#[allow(dead_code)]
pub fn create_ancestor_generator<const S: usize>(
    sequence_len: usize,
    variant_sites: &[[u8; S]],
) -> AncestorGenerator {
    let variant_data = VariantDataBuilder::from_iter(
        sequence_len,
        variant_sites
            .iter()
            .enumerate()
            .map(|(i, site)| (site.to_vec(), i + 1, /* non-sense data */ 'T', 'A')),
    )
    .finalize();

    AncestorGenerator::from_variant_data(variant_data)
}

#[allow(dead_code)]
pub fn generate_ancestors(ag: AncestorGenerator) -> AncestorArray {
    ag.generate_ancestors()
}

#[allow(dead_code)]
pub fn match_ancestors(ancestors: AncestorArray) -> Vec<TreeSequenceNode> {
    let mut ancestor_matcher =
        ViterbiMatcher::with_parallelism(ancestors.clone(), 1e-2, 1e-20, false, 1, 1, 1);
    ancestor_matcher.match_ancestors();
    let sequential_sequence = ancestor_matcher.get_tree_sequence();

    let mut ancestor_matcher = ViterbiMatcher::new(ancestors, 1e-2, 1e-20, false, 1);
    ancestor_matcher.match_ancestors();
    let parallel_sequence = ancestor_matcher.get_tree_sequence();

    assert_eq!(sequential_sequence.nodes, parallel_sequence.nodes);
    sequential_sequence.nodes
}
