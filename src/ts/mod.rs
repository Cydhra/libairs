use std::cmp::Ordering;
use std::collections::BinaryHeap;
use crate::ancestors::AncestralSequence;

pub struct TreeSequenceGenerator {
    ancestor_sequences: Vec<AncestralSequence>,
    recombination_probabilities: Vec<f64>,
    mismatch_probabilities: Vec<f64>,
}

impl TreeSequenceGenerator {
    pub fn new(mut ancestor_sequences: Vec<AncestralSequence>, recombination_rate: f64, mismatch_rate: f64, variant_positions: Vec<usize>) -> Self {
        // sort ancestors by age, oldest first
        ancestor_sequences.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());

        // TODO those values are only tsinfer's defaults, the actual calculation of the values
        //  works differently

        let mut recombination_probabilities = vec![recombination_rate; variant_positions.len()];
        recombination_probabilities[0] = 0f64;

        let mismatch_probabilities = vec![mismatch_rate; variant_positions.len()];

        Self {
            ancestor_sequences,
            recombination_probabilities,
            mismatch_probabilities,
        }
    }

    /// For a given [`AncestralSequence`] and a set of ancestor sequences, calculate the most likely
    /// copying path within the LS model using the Viterbi algorithm.
    // TODO the ancestors are probably required to be packed sequences instead of ancestral sequences
    // TODO ancestors are required to be TSNodeSequences, but those dont work yet
    fn find_hidden_path(
        &self,
        candidate: &AncestralSequence,
        mut sweep_line_queue: BinaryHeap<SweepEvent>,
    ) -> Vec<(usize, usize, usize)> {
        let num_ancestors = sweep_line_queue.len();
        let mut active_ancestors = Vec::with_capacity(num_ancestors);
        let mut next_event_position = sweep_line_queue.peek().unwrap().position;

        let mut likelihoods = vec![1f64; num_ancestors];
        let mut recombination_points = vec![vec![false; num_ancestors]; candidate.len()];
        let mut max_likelihoods = vec![0; candidate.len()];
        let candidate_start = candidate.start();

        for (site, &state) in candidate.site_iter() {
            // update active ancestors
            while next_event_position <= site {
                let event = sweep_line_queue.pop().unwrap();
                if event.is_start {
                    active_ancestors.push(event.ancestor_index);

                    let ancestor_end = SweepEvent { position: self.ancestor_sequences[event.ancestor_index].end, ancestor_index: event.ancestor_index, is_start: false };
                    sweep_line_queue.push(ancestor_end)
                } else {
                    active_ancestors.retain(|&ancestor| ancestor != event.ancestor_index);
                }
                next_event_position = sweep_line_queue.peek().map_or(usize::MAX, |e| e.position)
            }

            let mut max_site_likelihood = -1f64;
            let mut max_site_likelihood_ancestor: Option<usize> = None;

            let rho = self.recombination_probabilities[site];
            let mu = self.mismatch_probabilities[site];
            let k = active_ancestors.len() as f64;
            let num_alleles = 2f64; // TODO we might not want to hard-code this

            debug_assert!(active_ancestors.len() > 0, "no active ancestors at {}", site);
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

        for (site, _) in candidate.site_iter().rev().skip(1) {
            if recombination_points[site - candidate_start][ancestor_index] {
                nodes.push((ancestor_index, site, ancestor_coverage_end));
                assert_ne!(ancestor_index, max_likelihoods[site - 1 - candidate_start], "recombination point {} is not a recombination point", site);
                ancestor_index = max_likelihoods[site - 1 - candidate_start];
                ancestor_coverage_end = site;
            }
        }
        nodes.push((ancestor_index, 0, ancestor_coverage_end));

        nodes.reverse();
        nodes
    }

    pub fn generate_tree_sequence(&self) -> Vec<(Vec<u8>, Vec<(usize, usize, usize)>)> {
        let mut sweep_line_queue = BinaryHeap::new();

        let mut current_age = f64::INFINITY;
        let mut current_age_set = Vec::new();
        let mut tree: Vec<(Vec<u8>, Vec<(usize, usize, usize)>)> = Vec::new();

        // the first ancestor is the ancestral state and doesnt need to be processed
        current_age_set.push(0);
        tree.push((
            Vec::from(self.ancestor_sequences[0].haplotype()),
            vec![(0, 0, self.ancestor_sequences[0].len())],
        ));

        for (ancestor_index, ancestor) in self.ancestor_sequences.iter().enumerate().skip(1) {
            if ancestor.relative_age() < current_age {
                current_age_set.iter().for_each(|&index| {
                    sweep_line_queue.push(SweepEvent { position: self.ancestor_sequences[index].start(), ancestor_index: index, is_start: true });
                });
                current_age_set.clear();
                current_age = ancestor.relative_age();
            }

            tree.push((
                Vec::from(ancestor.haplotype()),
                self.find_hidden_path(&ancestor, sweep_line_queue.clone()),
            ));

            current_age_set.push(ancestor_index);
        }

        tree
    }
}

/// A single event in the sweep line algorithm.
#[derive(Debug, Eq, PartialEq, Clone)]
struct SweepEvent {
    position: usize,
    ancestor_index: usize,
    is_start: bool,
}

impl PartialOrd<Self> for SweepEvent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SweepEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        other.position.cmp(&self.position).then(other.ancestor_index.cmp(&self.ancestor_index))
    }
}

#[cfg(test)]
mod tests {
    use crate::ancestors::{AncestorGenerator, AncestralSequence};
    use crate::dna::VariantSite;
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

        let mut ancestors = ag.generate_ancestors();
        ancestors.sort_unstable_by(|a, b| b.relative_age().partial_cmp(&a.relative_age()).unwrap());

        let ancestor_length = ancestors[0].len();
        let mut ancestral_state = AncestralSequence::from_ancestral_state(ancestor_length, 1.0);
        ancestral_state.end = ancestor_length;
        ancestors.insert(0, ancestral_state);
        let ancestor_matcher = TreeSequenceGenerator::new(ancestors, 1e-2, 1e-20, vec![1, 2, 3, 4, 5, 6]);
        let ts = ancestor_matcher.generate_tree_sequence();

        // expected tree
        let expected = vec![
            (vec![0, 0, 0, 0, 0], vec![(0, 0, 5)]),
            (vec![1, 0, 0, 1, 0], vec![(0, 0, 5)]),
            (vec![0, 1, 1, 0, 0], vec![(0, 0, 5)]),
            (vec![1, 0, 0, 1, 1], vec![(1, 0, 5)]),
        ];

        assert_eq!(ts.len(), expected.len());
        for (expected_ancestor, expected_path) in expected {
            let (_, path) = ts.iter().find(|(ancestor, _)| ancestor == &expected_ancestor).unwrap();
            assert_eq!(path, &expected_path);
        }
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
        let ancestor_matcher = TreeSequenceGenerator::new(ancestors, 1e-2, 1e-20, vec![1, 2, 4, 5, 6, 7]);
        let ts = ancestor_matcher.generate_tree_sequence();

        // expected tree
        let expected = vec![
            (vec![0, 0, 0, 0, 0, 0], vec![(0, 0, 6)]),
            (vec![0, 0, 0, 0, 1, 1], vec![(0, 0, 6)]),
            (vec![1, 1, 0, 0, 0, 0], vec![(0, 0, 6)]),
            (vec![1, 1, 1, 0, 1, 1], vec![(2, 0, 4), (1, 4, 6)]),
            (vec![1, 1, 0, 1, 0, 0], vec![(2, 0, 6)]),
        ];

        assert_eq!(ts.len(), expected.len());
        for (expected_ancestor, expected_path) in expected {
            let (_, path) = ts.iter().find(|(ancestor, _)| ancestor == &expected_ancestor).unwrap();
            let mut path = path.clone();
            path.sort_unstable_by(|(_, a_start, _), (_, b_start, _)| b_start.cmp(a_start));
            assert_eq!(path, expected_path);
        }
    }
}
