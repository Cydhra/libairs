use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence};
use crate::ts::ancestor_iterator::AncestorIndex;
use crate::ts::partial_sequence::{PartialSequenceEdge, PartialTreeSequence};
use crate::ts::tree_sequence::TreeSequence;
use crate::variants::VariantIndex;
use rayon::iter::IndexedParallelIterator;
use rayon::iter::ParallelIterator;
use rayon::iter::{IntoParallelRefIterator, IntoParallelRefMutIterator};
use std::vec;

/// A matcher runs the viterbi algorithm for a set of sequences.
/// It will generate a tree sequence from an array of ancestral sequences and can then match
/// samples against the generated tree sequence.
///
/// It exposes two main methods for matching sequences:
/// - [`ViterbiMatcher::match_ancestors`]: Match the given ancestors into a partial tree sequence
/// - [`ViterbiMatcher::match_samples`]: Match a set of samples against the partial tree sequence
pub struct ViterbiMatcher {
    ancestors: AncestorArray,
    partial_tree_sequence: PartialTreeSequence,
    recombination_prob: f64,
    mutation_prob: f64,
    use_recompression_threshold: bool,
    inverse_recompression_threshold: usize,
}

impl ViterbiMatcher {
    /// Create a new matcher for the given ancestral sequences
    pub fn new(
        ancestors: AncestorArray,
        recombination_prob: f64,
        mutation_prob: f64,
        use_recompression_threshold: bool,
        inverse_recompression_threshold: usize,
    ) -> Self {
        let ancestor_count = ancestors.len();
        Self {
            ancestors,
            partial_tree_sequence: PartialTreeSequence::with_capacity(ancestor_count),
            recombination_prob,
            mutation_prob,
            use_recompression_threshold,
            inverse_recompression_threshold,
        }
    }

    /// Find a copying path for the given sequence given the ancestor iterator as a partial
    /// tree sequence.
    #[must_use]
    fn find_copy_path(
        &self,
        ancestor_iterator: &mut AncestorIndex,
        candidate: &AncestralSequence,
        limit_nodes: usize,
    ) -> (Vec<PartialSequenceEdge>, Vec<VariantIndex>) {
        let mut candidate_iter = candidate.site_iter();

        // index relative to candidate start to index the recombination and mutation arrays
        let mut candidate_site_index = 0;

        let mut max_likelihoods = vec![None; candidate.len()];

        let rho: f64 = self.recombination_prob;
        let mu: f64 = self.mutation_prob;

        let mut sites = ancestor_iterator.sites(candidate.start(), candidate.end(), limit_nodes);

        sites.for_each(|(site, marginal_tree)| {
            let (_, &state) = candidate_iter.next().unwrap();

            let mut max_site_likelihood = -1f64;
            let mut max_site_likelihood_ancestor: Option<Ancestor> = None;

            let k = (marginal_tree.num_nodes() + 1) as f64; // number of ancestors in tableau plus the virtual root
                                                            // probability that any one specific ancestor recombines to the current ancestors
            let prob_recomb = rho / k;
            // probability that none of the k-1 active ancestors recombines to the current ancestor
            let prob_no_recomb = 1f64 - rho + rho / k;

            let num_alleles = 2f64; // TODO we might not want to hard-code this
            let rev_mu = 1f64 - (num_alleles - 1f64) * mu;

            for ancestor_id in marginal_tree.nodes() {
                let ancestral_sequence = &self.ancestors[ancestor_id];
                let prob_no_recomb = *marginal_tree.likelihood(ancestor_id) * prob_no_recomb;

                let pt = if prob_no_recomb > prob_recomb {
                    prob_no_recomb
                } else {
                    *marginal_tree.recombination_site(ancestor_id, candidate_site_index) = true;
                    prob_recomb
                };

                let pe = if state == ancestral_sequence[site] {
                    rev_mu
                } else {
                    *marginal_tree.mutation_site(ancestor_id, candidate_site_index) = true;
                    mu
                };

                let likelihood = pt * pe;
                *marginal_tree.likelihood(ancestor_id) = likelihood;

                if likelihood > max_site_likelihood {
                    max_site_likelihood = likelihood;
                    max_site_likelihood_ancestor = Some(ancestor_id);
                } else if likelihood == max_site_likelihood
                    && self.ancestors[ancestor_id].relative_age()
                        > self.ancestors[max_site_likelihood_ancestor.unwrap()].relative_age()
                {
                    // always select the older one to mimic tsinfer behavior
                    max_site_likelihood_ancestor = Some(ancestor_id);
                }
            }

            // Apparently a measure to maintain numerical stability
            for ancestor in marginal_tree.nodes() {
                *marginal_tree.likelihood(ancestor) /= max_site_likelihood;
            }

            max_likelihoods[candidate_site_index] = max_site_likelihood_ancestor;

            candidate_site_index += 1;
        });

        let mut edges: Vec<PartialSequenceEdge> = Vec::new();
        let mut mutations: Vec<VariantIndex> = Vec::new();

        let mut current_ancestor = max_likelihoods[candidate_site_index - 1]
            .expect("no max likelihood ancestor found at last site");
        let mut ancestor_coverage_end = candidate.end();

        for (site, _) in candidate.site_iter().rev() {
            candidate_site_index -= 1;

            if sites.mutation_site(current_ancestor, candidate_site_index) {
                mutations.push(site);
            }

            if sites.recombination_site(current_ancestor, candidate_site_index)
                && max_likelihoods[candidate_site_index - 1]
                    .expect("no max likelihood ancestor found")
                    != current_ancestor
            {
                debug_assert!(self.ancestors[current_ancestor].start() <= site);
                debug_assert!(self.ancestors[current_ancestor].end() >= site);

                edges.push(PartialSequenceEdge::new(
                    site,
                    ancestor_coverage_end,
                    current_ancestor,
                ));

                current_ancestor = max_likelihoods[candidate_site_index - 1]
                    .expect("no max likelihood ancestor found");
                ancestor_coverage_end = site;
            }
        }
        edges.push(PartialSequenceEdge::new(
            candidate.start(),
            ancestor_coverage_end,
            current_ancestor,
        ));

        edges.reverse();
        mutations.reverse();

        (edges, mutations)
    }

    /// Generate a tree sequence from the given ancestor array.
    /// This will modify the internal state to represent the tree sequence for all ancestors
    /// within the array.
    pub fn match_ancestors(&mut self) {
        // TODO make this configurable
        let num_threads = std::thread::available_parallelism()
            .expect("could not determine number of available threads");

        let mut ancestor_iterators = vec![
            AncestorIndex::new(
                self.ancestors.len(),
                self.ancestors.get_num_variants(),
                self.use_recompression_threshold,
                self.inverse_recompression_threshold,
            );
            num_threads.into()
        ];

        // edges for root node
        let root = Ancestor(0);
        self.partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                self.ancestors[root].start(),
                self.ancestors[root].end(),
                root,
            )],
            vec![],
        );

        let mut current_ancestor_index = 1;
        for chunk in self.ancestors.as_slice()[current_ancestor_index..].chunks(num_threads.into())
        {
            debug_assert!(chunk.len() <= num_threads.into());

            let ancestor_range = current_ancestor_index..current_ancestor_index + chunk.len();

            // prepare the iterators for the chunk by adding the start and end of unmatched ancestors
            // to the iterator, so that the marginal trees are initialized correctly with free nodes
            for ancestor_id in ancestor_range.clone() {
                let ancestor_id = Ancestor(ancestor_id);

                let (start, end) = (
                    self.ancestors[ancestor_id].start(),
                    self.ancestors[ancestor_id].end(),
                );
                ancestor_iterators.iter_mut().for_each(|iterator| {
                    iterator.insert_free_node(ancestor_id, start, end);
                });
            }

            let mut results = chunk
                .par_iter()
                .zip(ancestor_range)
                .zip(ancestor_iterators.par_iter_mut())
                .map(|((ancestor, ancestor_index), ancestor_iterator)| {
                    let mut num_ancestors = 0;

                    // TODO precalculate this
                    for (_, old_ancestor) in self.ancestors.iter().take(ancestor_index) {
                        if old_ancestor.relative_age() > ancestor.relative_age() {
                            // TODO we can perform an overlap check here
                            num_ancestors += 1;
                        }
                    }

                    (
                        ancestor_index,
                        self.find_copy_path(ancestor_iterator, ancestor, num_ancestors),
                    )
                })
                .collect::<Vec<_>>();

            // increase the index start for the next chunk
            debug_assert!(results.len() == chunk.len());
            current_ancestor_index += results.len();

            // insert the results into the partial tree sequence and the iterators
            results.sort_by(|(a, _), (b, _)| a.cmp(b));
            for (ancestor_index, (edges, mutations)) in results {
                ancestor_iterators.iter_mut().for_each(|iterator| {
                    iterator.insert_sequence_node(Ancestor(ancestor_index), &edges, &mutations);
                });
                self.partial_tree_sequence.push(edges, mutations);
            }
        }
    }

    /// Insert a set of samples into an ancestral tree sequence. The samples will be matched
    /// against the existing sequence, but not against each other.
    pub fn match_samples(&mut self, samples: &[AncestralSequence]) {
        // TODO change this method to take the iterators as a parameter, so we dont need to
        //  reallocate them when we match ancestors and samples in one go anyway

        // TODO the iterators are empty, matching samples will fail!
        let num_threads = 1usize;

        let mut ancestor_iterators = vec![
            AncestorIndex::new(
                self.ancestors.len(),
                self.ancestors.get_num_variants(),
                self.use_recompression_threshold,
                self.inverse_recompression_threshold,
            );
            num_threads.into()
        ];

        for sample in samples {
            let (edges, mutations) =
                self.find_copy_path(&mut ancestor_iterators[0], sample, self.ancestors.len());
            self.partial_tree_sequence.push(edges, mutations);
        }
    }

    /// Finalize the tree sequence.
    pub fn get_tree_sequence(&self) -> TreeSequence {
        self.partial_tree_sequence.as_tree_sequence(&self.ancestors)
    }
}
