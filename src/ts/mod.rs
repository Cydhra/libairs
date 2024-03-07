mod ancestor_iterator;
mod matcher;
mod partial_sequence;
mod tree_sequence;

pub use matcher::ViterbiMatcher;

#[cfg(test)]
mod tests {
    use crate::ancestors::AncestorGenerator;
    use crate::ts::matcher::ViterbiMatcher;
    use crate::variants::{SequencePosition, VariantSite};
    use std::ops::Deref;

    #[test]
    fn trivial_tree_test() {
        let site1 = vec![0, 0, 0, 1, 1, 1];
        let site2 = vec![0, 1, 1, 0, 0, 0];
        let site3 = vec![0, 1, 1, 0, 0, 0];
        let site4 = vec![0, 0, 0, 1, 1, 1];
        let site5 = vec![0, 1, 0, 0, 0, 1];

        let ag = AncestorGenerator::from_iter(
            vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site3, 3),
                VariantSite::new(site4, 4),
                VariantSite::new(site5, 5),
            ]
            .into_iter(),
        );
        let len = SequencePosition::from_usize(6);

        let ancestors = ag.generate_ancestors(len);
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::new(ancestors, 1e-2, 1e-20, ag.variant_positions().len());
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

        let ag = AncestorGenerator::from_iter(
            vec![
                VariantSite::new(site1, 1),
                VariantSite::new(site2, 2),
                VariantSite::new(site4, 4),
                VariantSite::new(site5, 5),
                VariantSite::new(site6, 6),
                VariantSite::new(site7, 7),
            ]
            .into_iter(),
        );
        let len = SequencePosition::from_usize(7);

        let ancestors = ag.generate_ancestors(len);
        let ancestors_copy = ancestors.deref().clone();
        let mut ancestor_matcher =
            ViterbiMatcher::new(ancestors, 1e-2, 1e-20, ag.variant_positions().len());
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
