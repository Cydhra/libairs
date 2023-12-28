use crate::ancestors::AncestralSequence;
use crate::ts::SweepEventKind::Start;
use radix_heap::RadixHeapMap;
use std::cmp::{Ordering, Reverse};
use std::io::Write;

/// An interval in an ancestor that is covered by a parent node in the tree sequence.
/// The interval is defined by the start and (exclusive) end position of the interval and the index of the
/// parent node.
#[derive(Debug, Clone)]
pub struct TreeSequenceInterval {
    parent: usize,
    start: usize,
    end: usize,
}

impl TreeSequenceInterval {
    pub fn new(parent: usize, start: usize, end: usize) -> Self {
        Self { parent, start, end }
    }
}

/// A node in the tree sequence. The node is defined by the index of the ancestor sequence it
/// represents and a list of intervals that define what parent nodes cover the ancestor sequence.
#[derive(Debug, Clone)]
pub struct TreeSequenceNode {
    ancestor_index: usize,
    node_intervals: Vec<TreeSequenceInterval>,
}

impl TreeSequenceNode {
    pub fn new(ancestor_index: usize, node_intervals: Vec<TreeSequenceInterval>) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals,
        }
    }

    pub fn empty(ancestor_index: usize) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals: Vec::new(),
        }
    }
}

pub struct TreeSequenceGenerator {
    ancestor_sequences: Vec<AncestralSequence>,
    partial_tree_sequence: Vec<TreeSequenceNode>,
    recombination_probabilities: Vec<f64>,
    mismatch_probabilities: Vec<f64>,
}

impl TreeSequenceGenerator {
    pub fn new(
        mut ancestor_sequences: Vec<AncestralSequence>,
        recombination_rate: f64,
        mismatch_rate: f64,
        variant_positions: Vec<usize>,
    ) -> Self {
        let num_ancestors = ancestor_sequences.len();
        // sort ancestors by age, oldest first
        ancestor_sequences
            .sort_unstable_by(|a, b| b.relative_age()
                .partial_cmp(&a.relative_age()).unwrap()
            );

        // TODO those values are only tsinfer's defaults, the actual calculation of the values
        //  works differently

        let mut recombination_probabilities = vec![recombination_rate; variant_positions.len()];
        recombination_probabilities[0] = 0f64;

        let mismatch_probabilities = vec![mismatch_rate; variant_positions.len()];

        Self {
            ancestor_sequences,
            partial_tree_sequence: (0..num_ancestors)
                .map(|i| TreeSequenceNode::empty(i))
                .collect(),
            recombination_probabilities,
            mismatch_probabilities,
        }
    }

    /// For a given [`AncestralSequence`] and a set of ancestor sequences, calculate the most likely
    /// copying path within the LS model using the Viterbi algorithm.
    fn find_hidden_path(
        &self,
        candidate: &AncestralSequence,
        mut sweep_line_queue: RadixHeapMap<Reverse<usize>, SweepEvent>,
    ) -> Vec<TreeSequenceInterval> {
        let num_ancestors = sweep_line_queue.len();
        let mut active_ancestors = Vec::with_capacity(num_ancestors);
        let mut next_event_position = sweep_line_queue.peek_key().unwrap().0;

        let mut likelihoods = vec![1f64; num_ancestors];
        let mut recombination_points = vec![vec![false; num_ancestors]; candidate.len()];
        let mut max_likelihoods = vec![0; candidate.len()];
        let candidate_start = candidate.start();

        for (site, &state) in candidate.site_iter() {
            // update active ancestors
            while next_event_position <= site {
                let (_, event) = sweep_line_queue.pop().unwrap();
                if event.kind == Start {
                    active_ancestors.push(event.ancestor_index);

                    let ancestor_end = SweepEvent {
                        kind: SweepEventKind::End {
                            next_interval_index: 1,
                        },
                        position: self.partial_tree_sequence[event.ancestor_index].node_intervals[0].end,
                        ancestor_index: event.ancestor_index,
                    };
                    sweep_line_queue.push(Reverse(self.partial_tree_sequence[event.ancestor_index].node_intervals[0].end), ancestor_end)
                } else if let _end_event @ &SweepEventKind::End { next_interval_index } = &event.kind {
                    let node = &self.partial_tree_sequence[event.ancestor_index];
                    if node.node_intervals.len() > next_interval_index {
                        let next_interval_end_event = SweepEvent {
                            kind: SweepEventKind::End {
                                next_interval_index: next_interval_index + 1,
                            },
                            position: node.node_intervals[next_interval_index].end,
                            ancestor_index: event.ancestor_index,
                        };
                        sweep_line_queue.push(
                            Reverse(node.node_intervals[next_interval_index].end),
                            next_interval_end_event,
                        );

                        // todo
                        //  clear L-cache
                        //  change tree topology
                    } else {
                        active_ancestors.retain(|&ancestor| ancestor != event.ancestor_index);
                    }
                }
                next_event_position = sweep_line_queue.peek_key().map_or(usize::MAX, |e| e.0)
            }

            let mut max_site_likelihood = -1f64;
            let mut max_site_likelihood_ancestor: Option<usize> = None;

            let rho = self.recombination_probabilities[site];
            let mu = self.mismatch_probabilities[site];
            let k = active_ancestors.len() as f64;
            let num_alleles = 2f64; // TODO we might not want to hard-code this

            debug_assert!(
                active_ancestors.len() > 0,
                "no active ancestors at {}",
                site
            );
            for &ancestor_id in active_ancestors.iter() {
                let ancestral_sequence = &self.ancestor_sequences[ancestor_id];
                let prob_no_recomb = likelihoods[ancestor_id] * (1f64 - rho - rho / k);
                let prob_recomb = rho / k;

                let pt = if prob_no_recomb > prob_recomb {
                    prob_no_recomb
                } else {
                    recombination_points[site - candidate_start][ancestor_id] = true;
                    prob_recomb
                };

                let pe = if state == ancestral_sequence[site] {
                    1f64 - (num_alleles - 1f64) * mu
                } else {
                    mu
                };

                likelihoods[ancestor_id] = pt * pe;

                if likelihoods[ancestor_id] > max_site_likelihood {
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
        let mut ancestor_index = max_likelihoods[candidate.end() - 1 - candidate_start];
        let mut ancestor_coverage_end = candidate.end();

        for (site, _) in candidate.site_iter().rev() {
            if recombination_points[site - candidate_start][ancestor_index] {
                nodes.push(TreeSequenceInterval::new(
                    ancestor_index,
                    site,
                    ancestor_coverage_end,
                ));
                // assert_ne!(
                //     ancestor_index,
                //     max_likelihoods[site - 1 - candidate_start],
                //     "recombination point {} is not a recombination point",
                //     site
                // );
                ancestor_index = max_likelihoods[site - 1 - candidate_start];
                ancestor_coverage_end = site;
            }
        }
        nodes.push(TreeSequenceInterval::new(
            ancestor_index,
            0,
            ancestor_coverage_end,
        ));

        nodes.reverse();
        nodes
    }

    pub fn generate_tree_sequence(&mut self) -> Vec<TreeSequenceNode> {
        let mut sweep_line_queue = RadixHeapMap::<Reverse<usize>, SweepEvent>::new();

        let mut current_age = f64::INFINITY;
        let mut current_age_set = Vec::new();

        // the first ancestor is the ancestral state and doesnt need to be processed
        current_age_set.push(0);
        self.partial_tree_sequence[0]
            .node_intervals
            .push(TreeSequenceInterval::new(
                0,
                0,
                self.ancestor_sequences[0].len(),
            ));

        for (ancestor_index, ancestor) in self.ancestor_sequences.iter().enumerate().skip(1) {
            if ancestor.relative_age() < current_age {
                current_age_set.iter().for_each(|&index| {
                    sweep_line_queue.push(
                        Reverse(self.ancestor_sequences[index].start()),
                        SweepEvent {
                            kind: Start,
                            position: self.ancestor_sequences[index].start(),
                            ancestor_index: index,
                        },
                    );
                });
                current_age_set.clear();
                current_age = ancestor.relative_age();
            }

            let intervals = self.find_hidden_path(ancestor, sweep_line_queue.clone());
            self.partial_tree_sequence[ancestor_index]
                .node_intervals
                .extend(intervals);

            current_age_set.push(ancestor_index);
        }

        // TODO dont need to clone here if we consume the generator
        self.partial_tree_sequence.clone()
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
    End { next_interval_index: usize },
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
    use crate::ancestors::{AncestorGenerator, AncestralSequence};
    use crate::dna::VariantSite;
    use crate::ts::TreeSequenceGenerator;
    use std::fs::File;
    use std::hint::black_box;
    use std::io::{BufRead, BufReader};

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

        let mut ancestors = ag.generate_ancestors();
        ancestors.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());

        let ancestor_length = ancestors[0].len();
        let mut ancestral_state = AncestralSequence::from_ancestral_state(ancestor_length, 1.0);
        ancestral_state.end = ancestor_length;
        ancestors.insert(0, ancestral_state);
        let ancestors_copy = ancestors.clone();
        let mut ancestor_matcher =
            TreeSequenceGenerator::new(ancestors, 1e-2, 1e-20, vec![1, 2, 3, 4, 5, 6]);
        let ts = ancestor_matcher.generate_tree_sequence();

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

        let mut ancestors = ag.generate_ancestors();
        ancestors.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());

        let ancestor_length = ancestors[0].len();
        let mut ancestral_state = AncestralSequence::from_ancestral_state(ancestor_length, 1.0);
        ancestral_state.end = ancestor_length;
        ancestors.insert(0, ancestral_state);
        let ancestors_copy = ancestors.clone();
        let mut ancestor_matcher =
            TreeSequenceGenerator::new(ancestors, 1e-2, 1e-20, vec![1, 2, 4, 5, 6, 7]);
        let ts = ancestor_matcher.generate_tree_sequence();

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

    /// Reads a specially formatted text file that contains data about variant sites intended for
    /// unit testing. The data was generated by dumping it from tsinfer.
    fn read_variant_dump(path: &str) -> impl Iterator<Item = VariantSite> {
        let input = File::open(path).expect("could not find test data");
        let reader = BufReader::new(input);
        reader.lines().enumerate().map(|(pos, line)| {
            VariantSite::new(
                line.expect("unexpected io error")
                    .trim()
                    .split(" ")
                    .map(|s| s.parse().expect("corrupt input data"))
                    .collect::<Vec<_>>(),
                pos,
            )
        })
    }

    #[test]
    // #[ignore]
    fn compute_chr20_10k_variants() {
        let variant_sites = read_variant_dump("testdata/chr20_10k_variants.txt");
        let ag = AncestorGenerator::from_iter(variant_sites);
        let ancestors = ag.generate_ancestors();

        let mut ancestor_matcher = TreeSequenceGenerator::new(
            ancestors,
            1e-2,
            1e-20,
            ag.sites.iter().map(|s| s.position).collect(),
        );
        let ts = ancestor_matcher.generate_tree_sequence();
        black_box(ts);
    }
}
