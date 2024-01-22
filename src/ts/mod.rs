mod tree_sequence;

use crate::ancestors::AncestralSequence;
use crate::dna::SequencePosition;
use crate::ts::tree_sequence::{TreeSequence, TreeSequenceInterval, TreeSequenceNode};
use crate::ts::SweepEventKind::Start;
use radix_heap::RadixHeapMap;
use std::cmp::{Ordering, Reverse};
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;

pub struct TreeSequenceGenerator {
    pub ancestor_sequences: Vec<AncestralSequence>,
    partial_tree_sequence: Vec<TreeSequenceNode>,
    variant_positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
    recombination_prob: f64,
    mismatch_prob: f64,
}

impl TreeSequenceGenerator {
    pub fn new(
        ancestor_sequences: Vec<AncestralSequence>,
        sequence_length: SequencePosition,
        recombination_rate: f64,
        mismatch_rate: f64,
        variant_positions: Vec<SequencePosition>,
    ) -> Self {
        let num_ancestors = ancestor_sequences.len();

        Self {
            ancestor_sequences,
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
        &self,
        candidate: &AncestralSequence,
        mut sweep_line_queue: RadixHeapMap<Reverse<usize>, SweepEvent>,
        num_ancestors: usize,
    ) -> (Vec<TreeSequenceInterval>, Vec<usize>) {
        let mut active_ancestors = Vec::with_capacity(num_ancestors);
        let mut next_event_position = sweep_line_queue.peek_key().unwrap().0;

        let mut likelihoods = vec![1f64; num_ancestors];
        let mut recombination_points = vec![vec![false; num_ancestors]; candidate.len()];
        let mut mutation_points = vec![vec![false; num_ancestors]; candidate.len()];

        let mut max_likelihoods = vec![0; candidate.len()];
        let candidate_start = candidate.start();

        let rho = self.recombination_prob;
        let mu = self.mismatch_prob;

        for (site, &state) in candidate.site_iter() {
            // update active ancestors
            while next_event_position <= site {
                let (_, event) = sweep_line_queue.pop().unwrap();
                if event.kind == Start {
                    active_ancestors.push(event.ancestor_index);

                    // if the ancestor is not available from the start, the initial likelihood must be zero,
                    // so we recombine away from it on tracback
                    if event.position > candidate_start {
                        likelihoods[event.ancestor_index] = 0f64;
                    }
                    let end_event_pos = self.ancestor_sequences[event.ancestor_index].end();

                    let ancestor_end = SweepEvent {
                        kind: SweepEventKind::End,
                        position: end_event_pos,
                        ancestor_index: event.ancestor_index,
                    };
                    sweep_line_queue.push(Reverse(end_event_pos), ancestor_end)
                } else {
                    active_ancestors.retain(|&ancestor| ancestor != event.ancestor_index);
                }
                next_event_position = sweep_line_queue.peek_key().map_or(usize::MAX, |e| e.0);
            }

            let mut max_site_likelihood = -1f64;
            let mut max_site_likelihood_ancestor: Option<usize> = None;

            let k = (num_ancestors + 1) as f64; // number of ancestors in tableau plus the virtual root
            // probability that any one specific ancestor recombines to the current ancestors
            let prob_recomb = rho / k;
            // probability that none of the k-1 active ancestors recombines to the current ancestor
            let prob_no_recomb = 1f64 - rho + rho / k;

            let num_alleles = 2f64; // TODO we might not want to hard-code this
            let rev_mu = 1f64 - (num_alleles - 1f64) * mu;

            debug_assert!(
                active_ancestors.len() > 0,
                "no active ancestors at {}",
                site
            );
            for &ancestor_id in active_ancestors.iter() {
                let ancestral_sequence = &self.ancestor_sequences[ancestor_id];
                let prob_no_recomb = likelihoods[ancestor_id] * prob_no_recomb;

                let pt = if prob_no_recomb > prob_recomb {
                    prob_no_recomb
                } else {
                    recombination_points[site - candidate_start][ancestor_id] = true;
                    prob_recomb
                };

                let pe = if state == ancestral_sequence[site] {
                    rev_mu
                } else {
                    mutation_points[site - candidate_start][ancestor_id] = true;
                    mu
                };

                likelihoods[ancestor_id] = pt * pe;

                if likelihoods[ancestor_id] > max_site_likelihood {
                    max_site_likelihood = likelihoods[ancestor_id];
                    max_site_likelihood_ancestor = Some(ancestor_id);
                } else if likelihoods[ancestor_id] == max_site_likelihood
                    && self.ancestor_sequences[ancestor_id].relative_age()
                    > self.ancestor_sequences[max_site_likelihood_ancestor.unwrap()]
                    .relative_age()
                {
                    // TODO this is a hack to make sure that the oldest ancestor is chosen in case of a tie,
                    //  because this is what tsinfer does implicitly does by not calculating the likelihoods
                    //  of ancestors with equal sequences.
                    max_site_likelihood = likelihoods[ancestor_id];
                    max_site_likelihood_ancestor = Some(ancestor_id);
                }
            }

            // Apparently a measure to maintain numerical stability
            for &ancestor in active_ancestors.iter() {
                likelihoods[ancestor] /= max_site_likelihood;
            }

            max_likelihoods[site - candidate_start] =
                max_site_likelihood_ancestor.expect("no max likelihood calculated");
        }
        let mut nodes = Vec::new();
        let mut mutations = Vec::new();

        let mut ancestor_index = max_likelihoods[candidate.end() - 1 - candidate_start];
        let mut ancestor_coverage_end = if candidate.end() == self.variant_positions.len() {
            self.sequence_length
        } else {
            // convert exclusive site index to exclusive site position
            self.variant_positions[candidate.end()]
        };

        for (site, _) in candidate.site_iter().rev() {
            if mutation_points[site - candidate_start][ancestor_index] {
                mutations.push(site);
            }

            if recombination_points[site - candidate_start][ancestor_index]
                && max_likelihoods[site - 1 - candidate_start] != ancestor_index
            {
                nodes.push(TreeSequenceInterval::new(
                    ancestor_index,
                    self.variant_positions[site],
                    ancestor_coverage_end,
                ));

                ancestor_index = max_likelihoods[site - 1 - candidate_start];
                ancestor_coverage_end = self.variant_positions[site];
            }
        }
        nodes.push(TreeSequenceInterval::new(
            ancestor_index,
            if candidate.start() == 0 {
                SequencePosition::from_usize(0)
            } else {
                self.variant_positions[candidate.start()]
            },
            ancestor_coverage_end,
        ));

        nodes.reverse();
        mutations.reverse();
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

        for (ancestor_index, ancestor) in self.ancestor_sequences.iter().enumerate().skip(1) {
            let mut sweep_line_queue = RadixHeapMap::<Reverse<usize>, SweepEvent>::new();
            let mut num_ancestors = 0;

            for (old_ancestor_index, old_ancestor) in self.ancestor_sequences.iter().enumerate().take(ancestor_index) {
                if old_ancestor.relative_age() > ancestor.relative_age() { // TODO we can perform an overlap check here
                    num_ancestors += 1;
                    sweep_line_queue.push(
                        Reverse(old_ancestor.start()),
                        SweepEvent {
                            kind: Start,
                            position: old_ancestor.start(),
                            ancestor_index: old_ancestor_index,
                        },
                    );
                }
            }

            let (intervals, mutations) = self.find_hidden_path(ancestor, sweep_line_queue.clone(), num_ancestors);
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
        for ancestor in &self.ancestor_sequences {
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
    ancestor_index: usize,
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
        let ancestors_copy = ancestors.clone();
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
        let ancestors_copy = ancestors.clone();
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
