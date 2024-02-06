mod ancestor_iterator;
mod matcher;
mod partial_sequence;
mod tree_sequence;

use crate::ancestors::AncestorArray;
use crate::dna::SequencePosition;
use crate::ts::matcher::ViterbiMatcher;
use crate::ts::tree_sequence::TreeSequence;
use std::fs::File;
use std::io;
use std::io::Write;
use std::ops::Deref;
use std::path::Path;

// TODO remove this once the ViterbiMatcher can be used conveniently
pub struct TreeSequenceGenerator {
    pub ancestor_sequences: AncestorArray,
    matcher: ViterbiMatcher,
    variant_positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
}

impl TreeSequenceGenerator {
    pub fn new(
        ancestor_sequences: AncestorArray,
        sequence_length: SequencePosition,
        recombination_rate: f64,
        mismatch_rate: f64,
        variant_positions: Vec<SequencePosition>,
    ) -> Self {
        Self {
            ancestor_sequences: ancestor_sequences.clone(),
            matcher: ViterbiMatcher::new(ancestor_sequences, recombination_rate, mismatch_rate),
            variant_positions,
            sequence_length,
        }
    }

    pub fn generate_tree_sequence(mut self) -> TreeSequence {
        let partial_sequence = self.matcher.match_ancestors();
        partial_sequence.as_tree_sequence(self.ancestor_sequences)
    }

    /// Export the tree sequence in a TSV format that can be read by the test suite
    pub fn export_ancestors(&self, path: &Path) -> io::Result<()> {
        let mut node_file = path.to_path_buf();
        node_file.push("ancestors.tsv");
        let mut writer = File::create(node_file)?;

        writer.write_fmt(format_args!("start\tend\tage\tfocal_sites\tstate\n"))?;
        for ancestor in self.ancestor_sequences.deref() {
            ancestor.export(&mut writer)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::AncestorGenerator;
    use crate::dna::{SequencePosition, VariantSite};
    use crate::ts::TreeSequenceGenerator;
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
        let ancestor_matcher = TreeSequenceGenerator::new(
            ancestors,
            SequencePosition::from_usize(5),
            1e-2,
            1e-20,
            SequencePosition::from_vec(vec![1, 2, 3, 4, 5]),
        );
        let ts = ancestor_matcher.generate_tree_sequence().0;

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
        let ancestor_matcher = TreeSequenceGenerator::new(
            ancestors,
            len,
            1e-2,
            1e-20,
            SequencePosition::from_vec(vec![1, 2, 4, 5, 6, 7]),
        );
        let ts = ancestor_matcher.generate_tree_sequence().0;

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
