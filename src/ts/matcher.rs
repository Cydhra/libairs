use std::cell::RefCell;
use std::cmp::min;
use std::num::NonZeroUsize;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::thread::available_parallelism;
use std::vec;

use rayon::prelude::ParallelSliceMut;
use rayon::scope;
use thread_local::ThreadLocal;

use crate::ancestors::{Ancestor, AncestorArray, AncestralSequence};
use crate::ts::ancestor_index::EdgeSequence;
use crate::ts::ancestor_index::{ViterbiEventKind, ViterbiIterator};
use crate::ts::partial_sequence::PartialSequenceEdge;
use crate::ts::tree_sequence::TreeSequence;
use crate::ts::PartialTreeSequence;
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
        let num_threads = available_parallelism()
            .unwrap_or_else(|err| {
                eprintln!(
                    "Error getting number of threads: {}. Defaulting to 1 thread.",
                    err
                );
                NonZeroUsize::new(1).unwrap()
            })
            .get() as u16;
        let per_thread = 4; // TODO sensible default

        Self::with_parallelism(
            ancestors,
            recombination_prob,
            mutation_prob,
            use_recompression_threshold,
            inverse_recompression_threshold,
            num_threads,
            per_thread,
        )
    }

    /// Create a new matcher with custom parallelism settings.
    pub fn with_parallelism(
        ancestors: AncestorArray,
        recombination_prob: f64,
        mutation_prob: f64,
        use_recompression_threshold: bool,
        inverse_recompression_threshold: u16,
        num_threads: u16,
        per_thread: u16,
    ) -> Self {
        let ancestor_count = ancestors.len();
        Self {
            ancestors,
            partial_tree_sequence: PartialTreeSequence::with_capacity(ancestor_count as usize),
            edge_sequence: EdgeSequence::new(),
            recombination_prob,
            mutation_prob,
            use_recompression_threshold,
            inverse_recompression_threshold,
            num_threads,
            per_thread,
        }
    }

    /// Update the matcher's internal state with the provided partial tree sequence.
    /// Any data currently in the matcher will be overridden.
    pub fn read_partial_tree_sequence(&mut self, partial_tree_sequence: PartialTreeSequence) {
        self.partial_tree_sequence = partial_tree_sequence;
        self.edge_sequence = EdgeSequence::from_tree_sequence(&self.partial_tree_sequence);
    }

    /// Find a copying path for the given sequence given the ancestor iterator as a partial
    /// tree sequence.
    #[must_use]
    fn find_copy_path(
        &self,
        ancestor_iterator: &mut ViterbiIterator,
        candidate: &AncestralSequence,
        limit_nodes: u32,
    ) -> (Vec<PartialSequenceEdge>, Vec<VariantIndex>) {
        let mut candidate_iter = candidate.site_iter();

        // index relative to candidate start to index the recombination and mutation arrays
        let mut candidate_site_index: u32 = 0;

        let mut max_likelihoods = vec![None; candidate.len() as usize];

        let mu: f64 = self.mutation_prob;

        let mut sites = ancestor_iterator.iter_sites(
            &self.edge_sequence,
            candidate.start(),
            candidate.end(),
            limit_nodes,
        );

        sites.for_each(|(site, marginal_tree, nodes)| {
            // no recombination at the first genome site (I don't know whether results change if this special case is removed, we keep it for fidelity)
            let rho: f64 = if site == VariantIndex::from_usize(0) {
                0.0
            } else {
                self.recombination_prob
            };

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

            for &ancestor_id in nodes {
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
            for &ancestor in nodes {
                *marginal_tree.likelihood(ancestor) /= max_site_likelihood;
            }

            max_likelihoods[candidate_site_index as usize] = max_site_likelihood_ancestor;

            candidate_site_index += 1;
        });

        let mut edges: Vec<PartialSequenceEdge> = Vec::new();
        let mut mutations: Vec<VariantIndex> = Vec::new();

        let mut current_ancestor = max_likelihoods[(candidate_site_index - 1) as usize]
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
                    if candidate_site_index > 0
                        && max_likelihoods[(candidate_site_index - 1) as usize]
                            .expect("no max likelihood ancestor found")
                            != current_ancestor
                    {
                        edges.push(PartialSequenceEdge::new(
                            event.site,
                            ancestor_coverage_end,
                            current_ancestor,
                        ));

                        current_ancestor = max_likelihoods[(candidate_site_index - 1) as usize]
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
        let mut iterator = ViterbiIterator::new(
            self.ancestors.len(),
            self.use_recompression_threshold,
            self.inverse_recompression_threshold,
        );

        if self.num_threads > 1 {
            self.do_match_ancestors(&iterator)
        } else {
            self.do_match_ancestors_sequentially(&mut iterator)
        }
    }

    fn do_match_ancestors(&mut self, ancestor_iterator: &ViterbiIterator) {
        // edges for root node
        let root = Ancestor(0);
        self.partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                self.ancestors[root].start(),
                self.ancestors[root].end(),
                root,
            )],
            vec![],
            true,
        );

        let (queue_sender, queue_receiver) =
            flume::bounded(self.num_threads as usize * self.per_thread as usize);

        let thread_local_iterator = Arc::new(ThreadLocal::with_capacity(self.num_threads.into()));

        let mut current_ancestor_index = 1;

        while current_ancestor_index < self.ancestors.len() {
            let chunk_size = min(
                self.num_threads as u32 * self.per_thread as u32,
                self.ancestors.len() - current_ancestor_index,
            );

            for i in current_ancestor_index as u32..(current_ancestor_index + chunk_size) as u32 {
                queue_sender.send(i).expect("queue size insufficient");
            }

            // insert the next chunk of ancestors into the edge sequence as free nodes
            for ancestor_index in current_ancestor_index..(current_ancestor_index + chunk_size - 1)
            {
                let ancestor_index = Ancestor(ancestor_index as u32);
                self.edge_sequence.insert_free_node(
                    ancestor_index,
                    self.ancestors[ancestor_index].start(),
                    self.ancestors[ancestor_index].end(),
                )
            }

            let (sender, receiver) = channel();

            scope(|scope| {
                (0..self.num_threads).for_each(|_| {
                    let thread_local_iterator = &thread_local_iterator;
                    let original_ancestor_iterator = ancestor_iterator;
                    let ancestors = &self.ancestors;
                    let matcher = &self;
                    let sender = sender.clone();
                    let queue = &queue_receiver;

                    scope.spawn(move |_| {
                        let mut iterator = thread_local_iterator
                            .get_or(|| RefCell::new(original_ancestor_iterator.clone()))
                            .borrow_mut();

                        let mut res = queue.try_recv();
                        while let Ok(next_ancestor_index) = res {
                            let ancestor_index = Ancestor(next_ancestor_index);

                            if ancestor_index.0 >= ancestors.len() as u32 {
                                break;
                            }

                            let ancestor = &ancestors[ancestor_index];

                            let mut num_ancestors = 0;
                            // TODO precalculate this
                            for (_, old_ancestor) in
                                ancestors.iter().take(ancestor_index.0 as usize)
                            {
                                if old_ancestor.relative_age() > ancestor.relative_age() {
                                    // TODO we can perform an overlap check here
                                    num_ancestors += 1;
                                }
                            }

                            sender
                                .send((
                                    ancestor_index,
                                    matcher.find_copy_path(&mut iterator, ancestor, num_ancestors),
                                ))
                                .expect("failed to send result");

                            res = queue.try_recv();
                        }
                    });
                });
            });

            let mut results = receiver
                .iter()
                .take(chunk_size as usize)
                .collect::<Vec<_>>();
            results.sort_by(|(a, _), (b, _)| a.cmp(b));
            current_ancestor_index += results.len() as u32;
            // println!(
            //     "Progress: {}/{}",
            //     current_ancestor_index,
            //     self.ancestors.len()
            // );

            for (ancestor_index, (edges, mutations)) in results {
                self.edge_sequence
                    .insert_sequence_node(ancestor_index, &edges, &mutations);
                self.partial_tree_sequence.push(edges, mutations, true);
            }
        }
    }

    fn do_match_ancestors_sequentially(&mut self, iterator: &mut ViterbiIterator) {
        // edges for root node
        let root = Ancestor(0);
        self.partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                self.ancestors[root].start(),
                self.ancestors[root].end(),
                root,
            )],
            vec![],
            true,
        );

        for (ancestor_index, ancestor) in self.ancestors.iter().skip(1) {
            let mut num_ancestors = 0;

            // TODO precalculate this
            for (_, old_ancestor) in self.ancestors.iter().take(ancestor_index.0 as usize) {
                if old_ancestor.relative_age() > ancestor.relative_age() {
                    // TODO we can perform an overlap check here
                    num_ancestors += 1;
                }
            }

            let (edges, mutations) = self.find_copy_path(iterator, ancestor, num_ancestors);
            self.edge_sequence
                .insert_sequence_node(ancestor_index, &edges, &mutations);
            self.partial_tree_sequence.push(edges, mutations, true);
        }
    }

    /// Insert a set of samples into an ancestral tree sequence. The samples will be matched
    /// against the existing sequence, but not against each other.
    pub fn match_samples(&mut self) {
        let mut iterator = ViterbiIterator::new(
            self.ancestors.len(),
            self.use_recompression_threshold,
            self.inverse_recompression_threshold,
        );

        if self.num_threads > 1 {
            self.do_match_samples(&iterator)
        } else {
            self.do_match_sample_sequentially(&mut iterator)
        }
    }

    fn do_match_samples(&mut self, iterator: &ViterbiIterator) {
        let samples = self.ancestors.generate_sample_data();
        // assign samples to chunk in a way that all samples get assigned and one chunk may get less

        let (queue_sender, queue_receiver) = flume::bounded(samples.len());
        for i in 0..samples.len() {
            queue_sender.send(i).expect("queue size insufficient");
        }

        let (sender, receiver) = flume::bounded(samples.len());

        let thread_local_iterators = Arc::new(ThreadLocal::with_capacity(self.num_threads.into()));

        scope(|s| {
            let recv = &queue_receiver;
            let num_threads = self.num_threads as usize;
            let num_ancestors = self.ancestors.len();
            let thread_local_iterators = &thread_local_iterators;
            let result_sender = &sender;
            let samples = &samples;
            let matcher = &self;

            (0..self.num_threads).for_each(|_| {
                s.spawn(move |_| {
                    let mut iterator = thread_local_iterators
                        .get_or(|| RefCell::new(iterator.clone()))
                        .borrow_mut();
                    let mut results =
                        Vec::with_capacity((samples.len() + (num_threads - 1)) / num_threads);
                    while let Ok(i) = recv.try_recv() {
                        let sample = &samples[i];
                        results.push((
                            i,
                            matcher.find_copy_path(&mut iterator, &sample, num_ancestors),
                        ));
                    }

                    result_sender.send(results).unwrap();
                })
            });
        });
        let mut results = receiver
            .iter()
            .take(self.num_threads as usize)
            .flatten()
            .collect::<Vec<_>>();

        results.par_sort_by(|(a, _), (b, _)| a.cmp(b));

        // insert results into the partial tree sequence
        for (_, (edges, mutations)) in results {
            self.partial_tree_sequence.push(edges, mutations, false);
        }
    }

    fn do_match_sample_sequentially(&mut self, iterator: &mut ViterbiIterator) {
        let samples = self.ancestors.generate_sample_data();
        for sample in samples {
            let (edges, mutations) = self.find_copy_path(iterator, &sample, self.ancestors.len());
            self.partial_tree_sequence.push(edges, mutations, false);
        }
    }

    /// Get the current tree sequence state that was generated by the matcher.
    /// It can be serialized and later used to restore the matcher state.
    pub fn get_partial_tree_sequence(&self) -> &PartialTreeSequence {
        &self.partial_tree_sequence
    }

    /// Finalize the tree sequence.
    pub fn get_tree_sequence(&self) -> TreeSequence {
        self.partial_tree_sequence.as_tree_sequence(&self.ancestors)
    }
}
