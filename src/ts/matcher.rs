use std::cmp::min;
use std::sync::mpsc::channel;
use std::vec;

use rayon::iter::ParallelBridge;
use rayon::iter::ParallelIterator;
use scoped_thread_pool::Pool;

use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence};
use crate::ts::ancestor_index::EdgeSequence;
use crate::ts::ancestor_index::{ViterbiEventKind, ViterbiIterator};
use crate::ts::partial_sequence::{PartialSequenceEdge, PartialTreeSequence};
use crate::ts::tree_sequence::TreeSequence;
use crate::variants::VariantIndex;

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
    edge_sequence: EdgeSequence,
    recombination_prob: f64,
    mutation_prob: f64,
    use_recompression_threshold: bool,
    inverse_recompression_threshold: u16,
    num_threads: u16,
    per_thread: u16,
}

impl ViterbiMatcher {
    /// Create a new matcher for the given ancestral sequences
    pub fn new(
        ancestors: AncestorArray,
        recombination_prob: f64,
        mutation_prob: f64,
        use_recompression_threshold: bool,
        inverse_recompression_threshold: u16,
    ) -> Self {
        let ancestor_count = ancestors.len();
        Self {
            ancestors,
            partial_tree_sequence: PartialTreeSequence::with_capacity(ancestor_count),
            edge_sequence: EdgeSequence::new(),
            recombination_prob,
            mutation_prob,
            use_recompression_threshold,
            inverse_recompression_threshold,
            num_threads: 4, // fixme
            per_thread: 4,
        }
    }

    /// Find a copying path for the given sequence given the ancestor iterator as a partial
    /// tree sequence.
    #[must_use]
    fn find_copy_path(
        &self,
        ancestor_iterator: &mut ViterbiIterator,
        candidate: &AncestralSequence,
        limit_nodes: usize,
    ) -> (Vec<PartialSequenceEdge>, Vec<VariantIndex>) {
        let mut candidate_iter = candidate.site_iter();

        // index relative to candidate start to index the recombination and mutation arrays
        let mut candidate_site_index = 0;

        let mut max_likelihoods = vec![None; candidate.len()];

        let rho: f64 = self.recombination_prob;
        let mu: f64 = self.mutation_prob;

        let mut sites = ancestor_iterator.iter_sites(
            &self.edge_sequence,
            candidate.start(),
            candidate.end(),
            limit_nodes,
        );

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
                    marginal_tree.insert_recombination_event(ancestor_id, site);
                    prob_recomb
                };

                let pe = if state == ancestral_sequence[site] {
                    rev_mu
                } else {
                    marginal_tree.insert_mutation_event(ancestor_id, site);
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
        let mut last_site = candidate.end();

        // iterate through the viterbi events of the current ancestor until we find a recombination
        // event that coincides with a new most likely ancestor and then switch.
        let mut viterbi_event_iter = sites.traceback(&self.partial_tree_sequence, current_ancestor);

        while let Some(event) = viterbi_event_iter.next() {
            candidate_site_index -= last_site.get_variant_distance(event.site);
            last_site = event.site;

            match event.kind {
                ViterbiEventKind::Mutation => mutations.push(event.site),
                ViterbiEventKind::Recombination => {
                    if max_likelihoods[candidate_site_index - 1]
                        .expect("no max likelihood ancestor found")
                        != current_ancestor
                    {
                        edges.push(PartialSequenceEdge::new(
                            event.site,
                            ancestor_coverage_end,
                            current_ancestor,
                        ));

                        current_ancestor = max_likelihoods[candidate_site_index - 1]
                            .expect("no max likelihood ancestor found");
                        ancestor_coverage_end = event.site;

                        viterbi_event_iter.switch_ancestor(current_ancestor);
                    }
                }
                _ => unreachable!("unexpected viterbi event kind in find_copy_path"),
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
        let mut ancestor_iterators = vec![
            ViterbiIterator::new(
                self.ancestors.len(),
                self.use_recompression_threshold,
                self.inverse_recompression_threshold,
            );
            self.num_threads.into()
        ];

        self.do_match_ancestors(&mut ancestor_iterators)
    }

    fn do_match_ancestors(&mut self, ancestor_iterators: &mut [ViterbiIterator]) {
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
        let pool = Pool::new(self.num_threads.into());

        while current_ancestor_index < self.ancestors.len() {
            let chunk_size = min(
                self.num_threads as usize * self.per_thread as usize,
                self.ancestors.len() - current_ancestor_index,
            );

            // insert the next chunk of ancestors into the edge sequence as free nodes
            for ancestor_index in current_ancestor_index..(current_ancestor_index + chunk_size - 1)
            {
                let ancestor_index = Ancestor(ancestor_index);
                self.edge_sequence.insert_free_node(
                    ancestor_index,
                    self.ancestors[ancestor_index].start(),
                    self.ancestors[ancestor_index].end(),
                )
            }

            let (sender, receiver) = channel();

            pool.scoped(|scope| {
                ancestor_iterators.iter_mut().enumerate().for_each(
                    |(ancestor_offset, iterator)| {
                        let ancestors = &self.ancestors;
                        let matcher = &self;
                        let sender = sender.clone();
                        let per_thread = self.per_thread as usize;

                        scope.execute(move || {
                            for i in 0..per_thread {
                                let ancestor_index = Ancestor(
                                    current_ancestor_index + per_thread * ancestor_offset + i,
                                );

                                if ancestor_index.0 >= ancestors.len() {
                                    break;
                                }

                                let ancestor = &ancestors[ancestor_index];

                                let mut num_ancestors = 0;
                                // TODO precalculate this
                                for (_, old_ancestor) in ancestors.iter().take(ancestor_index.0) {
                                    if old_ancestor.relative_age() > ancestor.relative_age() {
                                        // TODO we can perform an overlap check here
                                        num_ancestors += 1;
                                    }
                                }

                                sender
                                    .send((
                                        ancestor_index,
                                        matcher.find_copy_path(iterator, ancestor, num_ancestors),
                                    ))
                                    .expect("failed to send result");
                            }
                        });
                    },
                );
            });

            let mut results = receiver.iter().take(chunk_size).collect::<Vec<_>>();
            results.sort_by(|(a, _), (b, _)| a.cmp(b));
            current_ancestor_index += results.len();

            for (ancestor_index, (edges, mutations)) in results {
                self.edge_sequence
                    .insert_sequence_node(ancestor_index, &edges, &mutations);
                self.partial_tree_sequence.push(edges, mutations);
            }
        }
    }

    /// Insert a set of samples into an ancestral tree sequence. The samples will be matched
    /// against the existing sequence, but not against each other.
    pub fn match_samples(&mut self) {
        let mut ancestor_iterators = vec![
            ViterbiIterator::new(
                self.ancestors.len(),
                self.use_recompression_threshold,
                self.inverse_recompression_threshold,
            );
            self.num_threads as usize
        ];

        self.do_match_samples(&mut ancestor_iterators)
    }

    fn do_match_samples(&mut self, ancestor_iterators: &mut [ViterbiIterator]) {
        let samples = self.ancestors.generate_sample_data();
        // assign samples to chunk in a way that all samples get assigned and one chunk may get less
        let samples_per_chunk =
            (samples.len() + (ancestor_iterators.len() - 1)) / ancestor_iterators.len();
        let chunks = samples.as_slice().chunks(samples_per_chunk);

        debug_assert!(ancestor_iterators.len() == chunks.len());

        let results = ancestor_iterators
            .into_iter()
            .zip(chunks)
            .par_bridge()
            .map(|(ancestor_index, samples)| {
                samples
                    .into_iter()
                    .map(|sample| self.find_copy_path(ancestor_index, sample, self.ancestors.len()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        // insert results into the partial tree sequence and the iterators
        for results in results {
            for (edges, mutations) in results {
                self.partial_tree_sequence.push(edges, mutations);
            }
        }
    }

    /// Calls [`do_match_ancestors`] and [`do_match_samples`] in sequence, without requiring
    /// additional allocations in between.
    /// This is more efficient than calling both individually.
    pub fn infer_tree_sequence(&mut self) {
        // TODO make this configurable
        let num_threads = 4usize;

        let mut ancestor_iterators = vec![
            ViterbiIterator::new(
                self.ancestors.len(),
                self.use_recompression_threshold,
                self.inverse_recompression_threshold,
            );
            num_threads.into()
        ];

        self.do_match_ancestors(&mut ancestor_iterators);
        self.do_match_samples(&mut ancestor_iterators);
    }

    /// Finalize the tree sequence.
    pub fn get_tree_sequence(&self) -> TreeSequence {
        self.partial_tree_sequence.as_tree_sequence(&self.ancestors)
    }
}
