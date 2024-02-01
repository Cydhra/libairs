use crate::ancestors::{Ancestor, AncestorArray};
use crate::ancestors::{AncestralSequence, VariantIndex};
use crate::ts::ancestor_iterator::AncestorIndex;
use crate::ts::partial_sequence::PartialSequenceEdge;

/// A matcher runs the viterbi algorithm for a set of sequences.
/// It will generate a tree sequence from an array of ancestral sequences and can then match
/// samples against the generated tree sequence.
///
/// It exposes two main methods for matching sequences:
/// - [`ViterbiMatcher::match_ancestors`]: Match the given ancestors into a partial tree sequence
/// - [`ViterbiMatcher::match_samples`]: Match a set of samples against the partial tree sequence
pub struct ViterbiMatcher {
    ancestors: AncestorArray,
    ancestor_iterator: AncestorIndex,
    recombination_prob: f64,
    mutation_prob: f64,
}

impl ViterbiMatcher {
    /// Create a new matcher for the given ancestral sequences
    pub fn new(ancestors: AncestorArray, recombination_prob: f64, mutation_prob: f64) -> Self {
        let ancestor_iterator = AncestorIndex::new();
        Self {
            ancestors,
            ancestor_iterator,
            recombination_prob,
            mutation_prob,
        }
    }

    /// Find a copying path for the given sequence given the ancestor iterator as a partial
    /// tree sequence.
    pub fn find_copy_path(
        &self,
        candidate: &AncestralSequence,
        limit_nodes: usize,
    ) -> (Vec<PartialSequenceEdge>, Vec<VariantIndex>) {
        let mut candidate_iter = candidate.site_iter();

        // index relative to candidate start to index the recombination and mutation arrays
        let mut candidate_site_index = 0;

        let mut max_likelihoods = vec![None; candidate.len()];

        let rho: f64 = self.recombination_prob;
        let mu: f64 = self.mutation_prob;

        let mut sites =
            self.ancestor_iterator
                .sites(candidate.start(), candidate.end(), limit_nodes);

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

    // TODO remove this method
    pub(crate) fn insert_edges(
        &mut self,
        node: Ancestor,
        edges: Vec<PartialSequenceEdge>,
        mutations: Vec<VariantIndex>,
    ) {
        self.ancestor_iterator
            .insert_sequence_node(node, edges, mutations)
    }

    /// Generate a tree sequence from the given ancestor array.
    /// This will modify the internal state to represent the tree sequence for all ancestors
    /// within the array.
    pub fn match_ancestors(&mut self) {
        todo!()
    }

    /// Insert a set of samples into an ancestral tree sequence. The samples will be matched
    /// against the existing sequence, but not against each other.
    pub fn match_samples(&mut self, samples: &[AncestralSequence]) {
        todo!()
    }

    /// Finalize the tree sequence.
    pub fn get_tree_sequence(&self) -> ! {
        todo!()
    }
}
