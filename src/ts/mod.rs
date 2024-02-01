mod ancestor_iterator;
mod matcher;
mod partial_sequence;
mod tree_sequence;

use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence};
use crate::dna::SequencePosition;
use crate::ts::matcher::ViterbiMatcher;
use crate::ts::tree_sequence::{TreeSequence, TreeSequenceInterval, TreeSequenceNode};
use radix_heap::RadixHeapMap;
use std::cmp::{Ordering, Reverse};
use std::fs::File;
use std::io;
use std::io::Write;
use std::ops::Deref;
use std::path::Path;

pub struct TreeSequenceGenerator {
    pub ancestor_sequences: AncestorArray,
    matcher: ViterbiMatcher,
    partial_tree_sequence: Vec<TreeSequenceNode>,
    variant_positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
    recombination_prob: f64,
    mismatch_prob: f64,
}

impl TreeSequenceGenerator {
    pub fn new(
        ancestor_sequences: AncestorArray,
        sequence_length: SequencePosition,
        recombination_rate: f64,
        mismatch_rate: f64,
        variant_positions: Vec<SequencePosition>,
    ) -> Self {
        let num_ancestors = ancestor_sequences.len();

        Self {
            ancestor_sequences: ancestor_sequences.clone(),
            matcher: ViterbiMatcher::new(ancestor_sequences, recombination_rate, mismatch_rate),
            partial_tree_sequence: (0..num_ancestors)
                .map(|i| TreeSequenceNode::empty(i))
                .collect(),
            variant_positions,
            sequence_length,
            recombination_prob: recombination_rate,
            mismatch_prob: mismatch_rate,
        }
    }

    /// For a given [`AncestralSequence`] and a set of ancestor sequences, calculate the most likely
    /// copying path within the LS model using the Viterbi algorithm.
    fn find_hidden_path(
        &mut self,
        node: Ancestor,
        candidate: &AncestralSequence,
        mut sweep_line_queue: RadixHeapMap<Reverse<usize>, SweepEvent>,
        num_ancestors: usize,
    ) -> (Vec<TreeSequenceInterval>, Vec<usize>) {
        let (edges, mutations) = self.matcher.find_copy_path(candidate, num_ancestors);
        self.matcher
            .insert_edges(node, edges.clone(), mutations.clone());

        let mut nodes = Vec::new();
        for edge in edges {
            nodes.push(TreeSequenceInterval {
                parent: edge.parent().0,
                start: if edge.start().0 == 0 {
                    SequencePosition::from_usize(0)
                } else {
                    self.variant_positions[edge.start().0]
                },
                end: if edge.end().0 == self.variant_positions.len() {
                    self.sequence_length
                } else {
                    self.variant_positions[edge.end().0]
                },
            })
        }
        let mut mutations = mutations.iter().map(|m| m.0).collect::<Vec<_>>();

        (nodes, mutations)
    }

    pub fn generate_tree_sequence(mut self) -> TreeSequence {
        // the first ancestor is the ancestral state and doesnt need to be processed
        // current_age_set.push(0);
        self.partial_tree_sequence[0]
            .node_intervals
            .push(TreeSequenceInterval::new(
                0,
                SequencePosition::from_usize(0),
                self.sequence_length,
            ));

        let ancestor_sequences = self.ancestor_sequences.clone();
        for (ancestor_index, ancestor) in ancestor_sequences.iter().enumerate().skip(1) {
            let mut sweep_line_queue = RadixHeapMap::<Reverse<usize>, SweepEvent>::new();
            let mut num_ancestors = 0;

            for (old_ancestor_index, old_ancestor) in self
                .ancestor_sequences
                .iter()
                .enumerate()
                .take(ancestor_index)
            {
                if old_ancestor.relative_age() > ancestor.relative_age() {
                    // TODO we can perform an overlap check here
                    num_ancestors += 1;
                }
            }

            let (intervals, mutations) = self.find_hidden_path(
                Ancestor(ancestor_index),
                ancestor,
                sweep_line_queue.clone(),
                num_ancestors,
            );
            self.partial_tree_sequence[ancestor_index].node_intervals = intervals;
            self.partial_tree_sequence[ancestor_index].mutations = mutations;
        }

        TreeSequence(self.partial_tree_sequence, self.ancestor_sequences)
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

/// A single event in the sweep line algorithm.
#[derive(Debug, Eq, PartialEq, Clone)]
struct SweepEvent {
    kind: SweepEventKind,
    position: usize,
    ancestor_index: Ancestor,
}

#[derive(Debug, Eq, PartialEq, Clone)]
enum SweepEventKind {
    Start,
    End,
}

impl PartialOrd<Self> for SweepEvent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SweepEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .position
            .cmp(&self.position)
            .then(other.ancestor_index.cmp(&self.ancestor_index))
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

        let ancestors = ag.generate_ancestors();
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
        assert_eq!(ts[0].ancestor_index, 0);
        assert_eq!(ts[0].node_intervals.len(), 1);
        assert_eq!(ts[0].node_intervals[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![1, 0, 0, 1, 0])
            .unwrap();
        assert_eq!(seq.node_intervals.len(), 1);
        assert_eq!(seq.node_intervals[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![1, 0, 0, 1, 1])
            .unwrap();
        assert_eq!(seq.node_intervals.len(), 1);
        assert_eq!(seq.node_intervals[0].parent, 1);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![0, 1, 1, 0, 0])
            .unwrap();
        assert_eq!(seq.node_intervals.len(), 1);
        assert_eq!(seq.node_intervals[0].parent, 0);
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

        let ancestors = ag.generate_ancestors();
        let ancestors_copy = ancestors.deref().clone();
        let ancestor_matcher = TreeSequenceGenerator::new(
            ancestors,
            SequencePosition::from_usize(7),
            1e-2,
            1e-20,
            SequencePosition::from_vec(vec![1, 2, 4, 5, 6, 7]),
        );
        let ts = ancestor_matcher.generate_tree_sequence().0;

        assert_eq!(ts.len(), 5);

        assert_eq!(ts[0].ancestor_index, 0);
        assert_eq!(ts[0].node_intervals.len(), 1);

        let seq1 = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![1, 1, 0, 0, 0, 0])
            .unwrap();
        assert_eq!(seq1.node_intervals.len(), 1);
        assert_eq!(seq1.node_intervals[0].parent, 0);

        let seq2 = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![0, 0, 0, 0, 1, 1])
            .unwrap();
        assert_eq!(seq2.node_intervals.len(), 1);
        assert_eq!(seq2.node_intervals[0].parent, 0);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![1, 1, 1, 0, 1, 1])
            .unwrap();
        assert_eq!(seq.node_intervals.len(), 2);
        assert_eq!(seq.node_intervals[0].parent, seq1.ancestor_index);
        assert_eq!(seq.node_intervals[1].parent, seq2.ancestor_index);

        let seq = ts
            .iter()
            .find(|n| ancestors_copy[n.ancestor_index].haplotype() == vec![1, 1, 0, 1, 0, 0])
            .unwrap();
        assert_eq!(seq.node_intervals.len(), 1);
        assert_eq!(seq.node_intervals[0].parent, seq1.ancestor_index);
    }
}
