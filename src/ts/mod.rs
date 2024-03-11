mod ancestor_iterator;
mod matcher;
mod partial_sequence;
mod tree_sequence;

pub use matcher::ViterbiMatcher;

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
            vec![(site1, 1), (site2, 2), (site3, 3), (site4, 4), (site5, 5)],
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::new(ancestors, 1e-2, 1e-20, ag.variant_data.len());
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
                (site1, 1),
                (site2, 2),
                (site4, 4),
                (site5, 5),
                (site6, 6),
                (site7, 7),
            ]
            .into_iter(),
        )
        .finalize();

        let ag = AncestorGenerator::from_variant_data(variant_data);
        let ancestors = ag.generate_ancestors();
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::new(ancestors, 1e-2, 1e-20, ag.variant_data.len());
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
        assert_eq!(seq.edges()[0].parent, seq1.ancestor());
        assert_eq!(seq.edges()[1].parent, seq2.ancestor());

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor()].haplotype() == vec![1, 1, 0, 1, 0, 0])
            .unwrap();
        assert_eq!(seq.edges().len(), 1);
        assert_eq!(seq.edges()[0].parent, seq1.ancestor());
    }
}
