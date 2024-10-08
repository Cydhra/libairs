mod ancestor_index;
mod matcher;
mod partial_sequence;
mod tree_sequence;

pub use matcher::ViterbiMatcher;
pub use partial_sequence::PartialTreeSequence;
pub use tree_sequence::*;

#[cfg(test)]
mod tests {
    use crate::ancestors::AncestorGenerator;
    use crate::ts::matcher::ViterbiMatcher;
    use crate::variants::VariantDataBuilder;
    use std::ops::Deref;

    #[test]
    fn trivial_tree_test() {
        let site1 = vec![0, 0, 0, 1, 1, 1];
        let site2 = vec![0, 1, 1, 0, 0, 0];
        let site3 = vec![0, 1, 1, 0, 0, 0];
        let site4 = vec![0, 0, 0, 1, 1, 1];
        let site5 = vec![0, 1, 0, 0, 0, 1];

        let variant_data = VariantDataBuilder::from_iter(
            6,
            vec![
                (site1, 1, 'C', 'A'),
                (site2, 2, 'C', 'T'),
                (site3, 3, 'C', 'T'),
                (site4, 4, 'C', 'A'),
                (site5, 5, 'C', 'G'),
            ],
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::with_parallelism(ancestors, 1e-2, 1e-20, false, 1, 1, 1);
        ancestor_matcher.match_ancestors();
        let ts = ancestor_matcher.get_tree_sequence().nodes;

        assert_eq!(ts.len(), 4);
        assert_eq!(ts[0].ancestor(), 0);
        assert_eq!(ts[0].edges().len(), 1);
        assert_eq!(ts[0].edges()[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 0, 0, 1, 0])
            .unwrap();
        assert_eq!(seq.edges().len(), 1);
        assert_eq!(seq.edges()[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 0, 0, 1, 1])
            .unwrap();
        assert_eq!(seq.edges().len(), 1);
        assert_eq!(seq.edges()[0].parent, 1);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![0, 1, 1, 0, 0])
            .unwrap();
        assert_eq!(seq.edges().len(), 1);
        assert_eq!(seq.edges()[0].parent, 0);
    }

    #[test]
    fn trivial_recombinant_test() {
        let site1 = vec![0, 0, 0, 1, 1, 1];
        let site2 = vec![0, 0, 0, 1, 1, 1];
        let site4 = vec![0, 1, 0, 0, 0, 1];
        let site5 = vec![0, 0, 0, 1, 1, 0];
        let site6 = vec![0, 1, 1, 0, 0, 1];
        let site7 = vec![0, 1, 1, 0, 0, 1];

        let variant_data = VariantDataBuilder::from_iter(
            7,
            vec![
                (site1, 1, 'T', 'A'),
                (site2, 2, 'T', 'A'),
                (site4, 4, 'T', 'A'),
                (site5, 5, 'T', 'A'),
                (site6, 6, 'T', 'A'),
                (site7, 7, 'T', 'A'),
            ]
            .into_iter(),
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::with_parallelism(ancestors, 1e-2, 1e-20, false, 1, 1, 1);
        ancestor_matcher.match_ancestors();
        let ts = ancestor_matcher.get_tree_sequence().nodes;

        assert_eq!(ts.len(), 5);

        assert_eq!(ts[0].ancestor(), 0);
        assert_eq!(ts[0].edges().len(), 1);

        let seq1 = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 1, 0, 0, 0, 0])
            .unwrap();
        assert_eq!(seq1.edges().len(), 1);
        assert_eq!(seq1.edges()[0].parent, 0);

        let seq2 = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![0, 0, 0, 0, 1, 1])
            .unwrap();
        assert_eq!(seq2.edges().len(), 1);
        assert_eq!(seq2.edges()[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 1, 1, 0, 1, 1])
            .unwrap();
        assert_eq!(seq.edges().len(), 2);
        assert_eq!(seq.edges()[0].parent as usize, seq1.ancestor());
        assert_eq!(seq.edges()[1].parent as usize, seq2.ancestor());

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 1, 0, 1, 0, 0])
            .unwrap();
        assert_eq!(seq.edges().len(), 1);
        assert_eq!(seq.edges()[0].parent as usize, seq1.ancestor());
    }
}
