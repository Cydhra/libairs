use std::iter::Peekable;
use std::num::NonZeroUsize;

use crate::ancestors::Ancestor;
use crate::ts::ancestor_index::tree::MarginalTree;
use crate::ts::partial_sequence::PartialTreeSequence;
use crate::variants::VariantIndex;
mod tree;

mod sequence;
mod traceback;

pub(in crate::ts) use sequence::{EdgeSequence, SequenceEvent};
pub(in crate::ts) use traceback::TracebackSequenceIterator;

/// A DNA sequence site visited during the viterbi algorithm.
/// Consists of a [`VariantIndex`] and a mutable reference to the marginal tree at that site.
pub(super) type Site<'a, 'o> = (VariantIndex, &'a mut MarginalTree<'o>, &'a Vec<Ancestor>);

/// Events that a single ancestor can experience during the Viterbi algorithm.
/// These events are generated during the forward search and are used to reconstruct
/// the most likely path during backtracing by storing where to look for mutations and recombinations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ViterbiEvent {
    pub(super) kind: ViterbiEventKind,
    pub(super) site: VariantIndex,

    // preceding viterbi event for the same ancestor. Cannot be zero, because this index is reserved
    // for ancestors that are never decompressed
    pub(super) prev: Option<NonZeroUsize>,
}

/// The kind of event that a single ancestor can experience during the Viterbi algorithm.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ViterbiEventKind {
    /// Sentinel event, that is only allowed as the first event in the linked list, and should never
    /// turn up during traceback. This allows several optimizations in the memory layout, since
    /// the value 0 is now an invalid index into the linked lists.
    Sentinel,

    /// Mutation event here
    Mutation,

    /// Recombination event here
    Recombination,

    /// Beginning from here, traceback is impossible from the current ancestor
    /// The compressed_until argument is inclusive
    Compressed,

    Uncompressed,
}

/// A helper structure for the Viterbi algorithm that helps to iterate through the sites, updating
/// the marginal tree at the currently visited site and helping with tree compression.
/// The iterator only produces indices into the ancestor array and sequences, it does not hold the
/// actual ancestors.
///
/// Whenever new tree edges are calculated, they can be inserted into the iterator so the tree
/// compression improves for future iterations.
///
/// To insert nodes that have no edges yet, call [`ViterbiIterator::insert_free_node`]. This will
/// make the ancestor available during Viterbi without compressing it into the tree.
///
/// To upgrade nodes with edges (enabling tree compression during the Viterbi algorithm), call
/// [`ViterbiIterator::insert_sequence_node`].
///
/// The iterator will start at the provided start point and iterate through all sites,
/// updating the marginal tree for each site and then calling a consumer function
/// with the [`VariantIndex`] and the [`MarginalTree`] (see: [`Site`])
#[derive(Clone)]
pub(super) struct ViterbiIterator {
    /// The actual parents of the tree nodes. `None` for free nodes and the root, parents might be
    /// compressed.
    parents: Vec<Option<Ancestor>>,

    /// Node children
    children: Vec<Vec<Ancestor>>,

    /// The most recent uncompressed ancestor node of each node in the marginal tree.
    /// `None` for root and free nodes
    uncompressed_parents: Vec<Option<Ancestor>>,

    /// Whether the node is currently compressed.
    is_compressed: Vec<bool>,

    /// The list of active nodes, i.e. nodes that are not compressed and are used in the Viterbi
    /// algorithm.
    active_nodes: Vec<Ancestor>,

    /// The likelihood of each node in the tree. Likelihood for compressed nodes is wrong, but
    /// they are ignored during the Viterbi algorithm.
    likelihoods: Vec<f64>,

    /// A super-array containing several interleaved linked lists with viterbi events for all
    /// ancestors. Every event is pushed to the back of this list, and connected to the last event
    /// inserted for the same ancestor. That event is retrieved from the `last_event` array.
    linked_viterbi_events: Vec<ViterbiEvent>,

    /// A mapping of ancestors to their last event index in viterbi_event_lists. Used as starting
    /// points for traceback.
    last_event: Vec<usize>,

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: Vec<VariantIndex>,

    /// Whether to use the recompression threshold to avoid recompressing when only a few nodes are
    /// active. If this is false, the algorithm will recompress in every iteration, regardless of
    /// efficacy.
    use_recompression_threshold: bool,

    /// Inverse recompression threshold.
    /// A threshold of 100 means that recompression of the marginal tree is attempted when
    /// more than 1% of the nodes are active,
    /// a threshold of 50 means that recompression is attempted when more than 2% of the nodes are active.
    /// A threshold of 1 means that recompression is attempted when more than 100% of the nodes are active,
    /// i.e. never.
    inv_recompression_threshold: u16,
}

impl ViterbiIterator {
    /// Create a new empty ancestor iterator.
    ///
    /// # Parameters
    /// - `max_nodes` the maximum number of tree nodes supported by the iterator. Memory allocations
    /// will orient themselves around this number.
    /// - `variant_count` the number of variants in the genome sequence
    /// - `use_recompression_threshold` whether to use the recompression threshold to avoid
    /// recompressing when only a few nodes are active. If this is false, the algorithm will
    /// recompress in every iteration, regardless of efficacy.
    /// - `inv_recompression_threshold` the inverse recompression threshold. The iterator will
    /// attempt to recompress the marginal tree when more than `1/inv_recompression_threshold * num_nodes`
    /// nodes of the marginal tree are active. If `use_recompression_threshold` is false, this
    /// parameter is ignored.
    pub(crate) fn new(
        max_nodes: u32,
        use_recompression_threshold: bool,
        inv_recompression_threshold: u16,
    ) -> Self {
        assert!(
            inv_recompression_threshold > 0 || !use_recompression_threshold,
            "Recompression interval must be greater than 0"
        );
        Self {
            parents: vec![None; max_nodes as usize],
            children: vec![Vec::new(); max_nodes as usize],
            uncompressed_parents: vec![None; max_nodes as usize],
            is_compressed: vec![true; max_nodes as usize],
            active_nodes: Vec::new(),
            likelihoods: vec![-1.0f64; max_nodes as usize],
            linked_viterbi_events: Vec::with_capacity(max_nodes as usize), // TODO this will grow strongly, so we should reserve much more
            last_event: vec![0; max_nodes as usize],
            last_compressed: vec![VariantIndex(0); max_nodes as usize],
            use_recompression_threshold,
            inv_recompression_threshold,
        }
    }

    /// Iterate through the variant sites of the ancestral tree sequence.
    /// The iterator will start at the provided start point and iterate through all sites
    /// until the given end index (exclusive).
    /// It provides the [`MarginalTree`] for each site.
    /// If multiple sites have the same marginal tree, it is provided multiple times.
    ///
    /// # Parameters
    /// - `start`: The index of the first site to iterate through
    /// - `end`: The index of the last site to iterate through (exclusive)
    /// - `limit_nodes`: The number of nodes to consider for the marginal tree. All nodes with a
    /// index equal or greater than this value will be ignored.
    pub(crate) fn iter_sites<'a>(
        &'a mut self,
        edge_sequence: &'a EdgeSequence,
        start: VariantIndex,
        end: VariantIndex,
        limit_nodes: u32,
    ) -> TreeSequenceState<impl Iterator<Item = &'a SequenceEvent>> {
        let mut marginal_tree = MarginalTree::new(
            start,
            edge_sequence.num_nodes,
            limit_nodes,
            &mut self.parents[0..limit_nodes as usize],
            &mut self.children[0..limit_nodes as usize],
            &mut self.uncompressed_parents[0..limit_nodes as usize],
            &mut self.is_compressed[0..limit_nodes as usize],
            &mut self.active_nodes,
            &mut self.likelihoods[0..limit_nodes as usize],
            &mut self.linked_viterbi_events,
            &mut self.last_event[0..limit_nodes as usize],
            &mut self.last_compressed[0..limit_nodes as usize],
            self.use_recompression_threshold,
            self.inv_recompression_threshold,
        );

        let site = start;
        let mut queue = edge_sequence.edge_index.iter().peekable();

        marginal_tree.advance_to_site(&mut self.active_nodes, &mut queue, site, true, false);

        TreeSequenceState::new(
            marginal_tree,
            &mut self.active_nodes,
            start,
            site,
            end,
            queue,
        )
    }
}

/// Borrowed iterator through the partial tree sequence held by an [`ViterbiIterator`].
/// Not an actual implementation of [`std::iter::IterMut`], because it borrows the same marginal tree for each site
/// mutably, so we need to make sure the reference cannot escape the iterator.
///
/// This means this struct does not actually implement Iterator traits, but offers a single method
/// [`TreeSequenceState::for_each`], which consumers have to call to process the tree sequence.
pub(super) struct TreeSequenceState<'a, 'o, I: Iterator<Item = &'a SequenceEvent>> {
    pub(super) marginal_tree: MarginalTree<'o>,
    pub(super) active_nodes: &'o mut Vec<Ancestor>,
    pub(super) start: VariantIndex,
    pub(super) site: VariantIndex,
    pub(super) end: VariantIndex,
    pub(super) queue: Peekable<I>,
}

impl<'a, 'o, I: Iterator<Item = &'a SequenceEvent>> TreeSequenceState<'a, 'o, I> {
    pub(super) fn new(
        marginal_tree: MarginalTree<'o>,
        active_nodes: &'o mut Vec<Ancestor>,
        start: VariantIndex,
        site: VariantIndex,
        end: VariantIndex,
        queue: Peekable<I>,
    ) -> Self {
        TreeSequenceState {
            marginal_tree,
            active_nodes,
            start,
            site,
            end,
            queue,
        }
    }

    pub(crate) fn for_each<'s, F: FnMut(Site)>(&'s mut self, mut consumer: F) {
        if self.site == self.end {
            return;
        }

        // first site has special treatment because we don't need to decompress edges that start here
        self.marginal_tree.advance_to_site(
            &mut self.active_nodes,
            &mut self.queue,
            self.site.next(),
            false,
            true,
        );
        consumer((self.site, &mut self.marginal_tree, self.active_nodes));
        self.marginal_tree
            .recompress_tree(&mut self.active_nodes, self.site);
        self.site = self.site.next();

        while self.site < self.end {
            self.marginal_tree.advance_to_site(
                &mut self.active_nodes,
                &mut self.queue,
                self.site.next(),
                false,
                false,
            );

            consumer((self.site, &mut self.marginal_tree, self.active_nodes));
            if !self.marginal_tree.use_recompression_threshold
                || self.active_nodes.len()
                    > (self.marginal_tree.limit_nodes
                        / self.marginal_tree.inv_recompression_threshold as u32)
                        as usize
            {
                self.marginal_tree
                    .recompress_tree(&mut self.active_nodes, self.site);
            }
            self.site = self.site.next();
        }

        // insert a compression event into all compressed ancestors, so traceback can follow the
        // compressed path
        for (node, &is_compressed) in self.marginal_tree.is_compressed.iter().enumerate() {
            if is_compressed
                && self.marginal_tree.last_compressed[node] < self.end
                && self.marginal_tree.parents[node].is_some()
            {
                // inlined version of `self.marginal_tree.insert_compression_event(Ancestor(node), self.end);`
                let last_compressed_begin = self.marginal_tree.last_compressed[node];

                self.marginal_tree.linked_viterbi_events.push(ViterbiEvent {
                    kind: ViterbiEventKind::Uncompressed,
                    site: last_compressed_begin,
                    prev: NonZeroUsize::new(self.marginal_tree.last_event_index[node]),
                });

                self.marginal_tree.linked_viterbi_events.push(ViterbiEvent {
                    kind: ViterbiEventKind::Compressed,
                    site: self.end,
                    prev: NonZeroUsize::new(self.marginal_tree.linked_viterbi_events.len() - 1),
                });
                if self.marginal_tree.last_event_index[node] > 0 {
                    debug_assert!(
                        last_compressed_begin
                            >= self.marginal_tree.linked_viterbi_events
                                [self.marginal_tree.last_event_index[node]]
                                .site,
                        "ancestor contains events despite being compressed: site {} < {:?}",
                        last_compressed_begin,
                        self.marginal_tree.linked_viterbi_events
                            [self.marginal_tree.last_event_index[node]],
                    );
                }

                self.marginal_tree.last_event_index[node] =
                    self.marginal_tree.linked_viterbi_events.len() - 1;
            }
        }
    }

    /// Create a traceback iterator through the Viterbi table generated during a previous call to
    /// `for_each`. The iterator will be initialized with `current_ancestor` for the traceback and
    /// automatically retrieve events from its correct ancestor whenever the actual ancestor is
    /// compressed.
    pub(crate) fn traceback<'s>(
        &'s self,
        partial_tree_sequence: &'s PartialTreeSequence,
        current_ancestor: Ancestor,
    ) -> TracebackSequenceIterator<'s, 'o>
    where
        's: 'o,
    {
        TracebackSequenceIterator::new(
            &self.marginal_tree,
            partial_tree_sequence,
            current_ancestor,
            self.start,
            self.end,
            self.marginal_tree.last_event_index[current_ancestor.0 as usize],
        )
    }
}

#[cfg(test)]
mod tests {
    use std::hint::black_box;

    use super::*;
    use crate::ts::ancestor_index::ViterbiEventKind::{Mutation, Recombination};
    use crate::ts::partial_sequence::PartialSequenceEdge;

    fn find_uncompressed_parent(tree: &MarginalTree, node: Ancestor) -> Option<Ancestor> {
        let mut parent = tree.parents[node.0 as usize];

        while let Some(p) = parent {
            if !tree.is_compressed(p) {
                break;
            }
            parent = tree.parents[p.0 as usize];
        }

        parent
    }

    #[test]
    fn test_empty_range() {
        // test whether iterating over an empty range doesn't crash and calls the closure zero times
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();

        edges.insert_free_node(
            Ancestor(1),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(10),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|_| panic!("closure called for empty range"));
    }

    #[test]
    fn test_ancestor_iteration() {
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        edges.insert_free_node(
            Ancestor(1),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );
        edges.insert_free_node(
            Ancestor(2),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(site, tree, nodes)| {
            assert_eq!(site, VariantIndex::from_usize(counter));
            assert_eq!(tree.num_nodes(), 2);
            assert_eq!(nodes.len(), 2);
            assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
            assert_eq!(find_uncompressed_parent(tree, Ancestor(1)), None);
            counter += 1;
        });
    }

    #[test]
    fn test_simple_tree_compression() {
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // insert edge from first to root node
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(5)],
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(site, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(site, VariantIndex::from_usize(counter));
            assert_eq!(tree.num_nodes(), 2);
            assert_eq!(
                nodes.len(),
                if counter < 5 { 1 } else { 2 },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1;
        });
    }

    #[test]
    fn test_simple_tree_recompression() {
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // insert edge from first to root node
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(site, tree, nodes)| {
            if site == VariantIndex::from_usize(0) {
                // fake likelihoods to prevent early recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }
            } else if site == VariantIndex::from_usize(4) {
                // trigger recompression exactly after site 4, meaning the node should be taken
                // away beginning from site 5
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1;
                }
            }

            assert_eq!(
                nodes.len(),
                if counter < 5 { 2 } else { 1 },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1;
        });
    }

    #[test]
    fn test_simple_tree_divergence() {
        // test whether a compressed node is decompressed when its parent changes to a different node
        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // insert edges for two nodes
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(9)],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            &vec![VariantIndex::from_usize(10)],
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(9),
            3,
        )
        .for_each(|(site, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(site, VariantIndex::from_usize(counter));
            assert_eq!(tree.num_nodes(), 3);
            assert_eq!(
                nodes.len(),
                if counter < 5 { 1 } else { 2 },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1;
        });
    }

    #[test]
    fn test_recompression_divergence() {
        // test whether divergence works after recompression
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            // diverge on site 1, then again on site 5, recompress in between
            &vec![VariantIndex::from_usize(1), VariantIndex::from_usize(5)],
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(_, _, nodes)| {
            // likelihoods stay at 0, so recompression happens immediately
            assert_eq!(
                nodes.len(),
                match counter {
                    1 => 2,
                    5 => 2,
                    _ => 1,
                },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1;
        });
    }

    #[test]
    fn test_mutation_on_recombination() {
        // test whether the correct nodes are decompressed when a mutation happens on a recombination
        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            &vec![VariantIndex::from_usize(5)],
        );

        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        )
        .for_each(|(_, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(
                nodes.len(),
                match counter {
                    0..=4 => 2,
                    5.. => 3,
                },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1usize;
        });
    }

    #[test]
    fn test_limit_nodes() {
        // test whether the correct nodes are in the tree if we limit the number of nodes
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // insert two nodes, one will be ignored
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // check the tree size is always 2
        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(_, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(nodes.len(), 2, "wrong number of nodes at site {}", counter);
            counter += 1usize;
        });
    }

    #[test]
    fn test_subset_iterator() {
        // test whether an incomplete range of sites iterated yields correct trees
        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 2;

        // insert two nodes
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(5)],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(1)],
        );

        // check the tree
        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(2),
            VariantIndex::from_usize(10),
            3,
        )
        .for_each(|(_, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(
                nodes.len(),
                match counter {
                    2..=4 => 1,
                    5.. => 2,
                    _ => unreachable!(),
                },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1usize;
        });
    }

    #[test]
    fn test_incomplete_ancestor() {
        // test whether an ancestor that ends early gets removed from the marginal tree
        let mut ix = ViterbiIterator::new(2, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // insert incomplete ancestor
        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(5),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // check the tree
        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|(_, tree, nodes)| {
            // fake likelihoods to prevent recompression
            for i in 0..tree.num_nodes() {
                tree.likelihoods[i] = 0.1 * i as f64
            }

            assert_eq!(
                nodes.len(),
                match counter {
                    0..=4 => 2,
                    5.. => 1,
                },
                "wrong number of nodes at site {}",
                counter
            );
            counter += 1usize;
        });
    }

    #[test]
    fn test_copying_from_parent_on_mutation() {
        // test whether the iterator copies the recombination and mutation sites from the parent when a child is
        // decompressed

        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // prepare partial tree sequence for traceback
        let mut partial_tree_sequence = PartialTreeSequence::with_capacity(3);
        partial_tree_sequence.push(vec![], vec![], true);
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
            false,
        );
        partial_tree_sequence.push(
            vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(0),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(10),
                    Ancestor(1),
                ),
            ],
            vec![VariantIndex::from_usize(4), VariantIndex::from_usize(9)],
            false,
        );

        // insert a child node that serves as a second node
        edges.insert_sequence_node(
            Ancestor(1),
            &partial_tree_sequence.edges[1],
            &partial_tree_sequence.mutations[1],
        );

        // insert another child node that copies from both the root and the child node
        edges.insert_sequence_node(
            Ancestor(2),
            &partial_tree_sequence.edges[2],
            &partial_tree_sequence.mutations[2],
        );

        // add mutations and recombinations with the root and first child node, and check whether the second child node
        // copies them
        let mut iterator = ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        );
        iterator.for_each(|(_, tree, nodes)| {
            match counter {
                // don't recompress first two nodes
                0 => nodes
                    .iter()
                    .copied()
                    .for_each(|n| tree.likelihoods[n.0 as usize] = n.0 as f64 * 0.1),

                // insert mutations and recombination for nodes 0 and 1, and later check if 2 inherited them
                2 => tree.insert_mutation_event(Ancestor(0), VariantIndex::from_usize(2)),
                3 => tree.insert_recombination_event(Ancestor(0), VariantIndex::from_usize(3)),
                4 => {
                    tree.likelihoods[2] = 0.1; // set likelihood of third node to that of second to trigger recompression
                    tree.insert_mutation_event(Ancestor(0), VariantIndex::from_usize(4))
                    // this shouldnt be copied, because the second node is uncompressed at this site
                }
                6 => tree.insert_mutation_event(Ancestor(1), VariantIndex::from_usize(6)),
                7 => tree.insert_recombination_event(Ancestor(1), VariantIndex::from_usize(7)),
                _ => {}
            }

            counter += 1;
        });

        let traceback_iterator = iterator
            .traceback(&partial_tree_sequence, Ancestor(2))
            .filter(|e| !matches!(e.kind, ViterbiEventKind::Uncompressed));
        assert_eq!(traceback_iterator.clone().count(), 4);

        traceback_iterator.for_each(|state| match state.site.0 {
            2 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            6 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            3 => assert_eq!(
                state.kind,
                ViterbiEventKind::Recombination,
                "expected recombination at {}",
                state.site.0
            ),
            7 => assert_eq!(
                state.kind,
                ViterbiEventKind::Recombination,
                "expected recombination at {}",
                state.site.0
            ),
            _ => assert_ne!(
                state.kind,
                ViterbiEventKind::Recombination,
                "expected no event at site {}",
                state.site.0
            ),
        });
    }

    #[test]
    fn test_copying_from_parent_on_change_parent() {
        // test whether the iterator copies the recombination and mutation sites from the parent when a child is
        // decompressed

        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // prepare partial tree sequence for traceback
        let mut partial_tree_sequence = PartialTreeSequence::with_capacity(3);
        partial_tree_sequence.push(vec![], vec![], true);
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
            false,
        );
        partial_tree_sequence.push(
            vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(0),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(7),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(7),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            vec![],
            false,
        );

        // insert a child node that serves as a second node
        edges.insert_sequence_node(
            Ancestor(1),
            &partial_tree_sequence.edges[1],
            &partial_tree_sequence.mutations[1],
        );

        // insert another child node that copies from both the root and the child node, but has no mutations
        edges.insert_sequence_node(
            Ancestor(2),
            &partial_tree_sequence.edges[2],
            &partial_tree_sequence.mutations[2],
        );

        // add mutations and recombinations with the root, and check whether the second child node copies them
        let mut iterator = ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        );
        iterator.for_each(|(_, tree, nodes)| {
            match counter {
                // don't recompress first two nodes
                0 => nodes
                    .iter()
                    .copied()
                    .for_each(|n| tree.likelihoods[n.0 as usize] = n.0 as f64 * 0.1),

                // insert mutations and recombination for nodes 0 and 1, and later check if 2 inherited them
                2 => tree.insert_mutation_event(Ancestor(0), VariantIndex::from_usize(2)),
                3 => tree.insert_recombination_event(Ancestor(0), VariantIndex::from_usize(3)),
                5 => tree.insert_mutation_event(Ancestor(1), VariantIndex::from_usize(5)), // this should not be copied because the second ancestor is uncompressed at this site

                _ => {}
            }

            counter += 1;
        });

        let traceback_iterator = iterator
            .traceback(&partial_tree_sequence, Ancestor(2))
            .filter(|e| !matches!(e.kind, ViterbiEventKind::Uncompressed));
        assert_eq!(traceback_iterator.clone().count(), 2);

        traceback_iterator.for_each(|state| match state.site.0 {
            2 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            3 => assert_eq!(
                state.kind,
                ViterbiEventKind::Recombination,
                "expected recombination at {}",
                state.site.0
            ),
            _ => assert_ne!(
                state.kind,
                ViterbiEventKind::Recombination,
                "expected no events at site {}",
                state.site.0
            ),
        });
    }

    #[test]
    fn test_traceback_compression() {
        // test whether a traceback iterator correctly returns the events that must be copied from
        // parents, as well as the events that are present in the child

        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        // prepare partial tree sequence for traceback
        let mut partial_tree_sequence = PartialTreeSequence::with_capacity(3);
        partial_tree_sequence.push(vec![], vec![], true);
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0), VariantIndex::from_usize(9)],
            true,
        );

        edges.insert_sequence_node(
            Ancestor(1),
            &partial_tree_sequence.edges[1],
            &partial_tree_sequence.mutations[1],
        );

        let mut iterator = ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            2,
        );
        iterator.for_each(|(site, tree, _)| {
            match counter {
                // insert a mutation on first site and then immediately recompress
                0 => tree.insert_mutation_event(Ancestor(1), site),

                // insert mutations for parent for sites 2, 3 and check if they are inherited
                2 => tree.insert_mutation_event(Ancestor(0), site),
                3 => tree.insert_mutation_event(Ancestor(0), site),

                // then one more for the child for site 9
                9 => tree.insert_mutation_event(Ancestor(1), site),
                _ => {}
            }

            counter += 1;
        });

        // check that we inherit all mutations from the parent and all from the child are also discovered
        let traceback_iterator = iterator
            .traceback(&partial_tree_sequence, Ancestor(1))
            .filter(|e| !matches!(e.kind, ViterbiEventKind::Uncompressed));
        assert_eq!(traceback_iterator.clone().count(), 4);

        traceback_iterator.for_each(|state| match state.site.0 {
            0 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            2 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            3 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            9 => assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            ),
            _ => assert_ne!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected no event at site {}",
                state.site.0
            ),
        });
    }

    #[test]
    fn test_transitive_compression() {
        // test whether the traceback iterator works correctly when copying from a grandparent

        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();

        // prepare partial tree sequence for traceback
        let mut partial_tree_sequence = PartialTreeSequence::with_capacity(3);
        partial_tree_sequence.push(vec![], vec![], true);
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(1)],
            true,
        );
        // the second child is uncompressed in first and last site
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(1),
            )],
            vec![VariantIndex::from_usize(0), VariantIndex::from_usize(9)],
            true,
        );

        edges.insert_sequence_node(
            Ancestor(1),
            &partial_tree_sequence.edges[1],
            &partial_tree_sequence.mutations[1],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &partial_tree_sequence.edges[2],
            &partial_tree_sequence.mutations[2],
        );

        let mut iterator = ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        );
        iterator.for_each(|(site, tree, nodes)| {
            // insert recombination on parent on site 0,
            // to make sure we not only copy mutations from the grandparent, but also stop
            // copying from it once the parent is decompressed
            if site.0 == 1 {
                tree.insert_recombination_event(Ancestor(1), site);
            }

            nodes.iter().copied().for_each(|n| {
                // make sure recompression happens instantly, meaning the third node is only present at first and last site
                *tree.likelihood(n) = 1.0;

                // insert mutations on every site
                tree.insert_mutation_event(n, site);
            });
        });

        // check that we inherit all mutations from the parent
        let traceback_iterator = iterator
            .traceback(&partial_tree_sequence, Ancestor(2))
            .filter(|e| !matches!(e.kind, ViterbiEventKind::Uncompressed));
        assert_eq!(traceback_iterator.clone().count(), 11);

        traceback_iterator.clone().take(9).for_each(|state| {
            assert_eq!(
                state.kind,
                ViterbiEventKind::Mutation,
                "expected mutation at {}",
                state.site.0
            )
        });

        assert_eq!(
            traceback_iterator.clone().skip(9).next().unwrap().kind,
            Recombination
        );
        assert_eq!(traceback_iterator.last().unwrap().kind, Mutation);
    }

    #[test]
    fn test_compression_interval_boundaries() {
        // we test if compressed intervals are accounted for, if their boundary is at the end of the range.
        // specifically, we test a regression where a compressed ancestor was not decompressed at the end of the range,
        // because the interval we tried to decompress only overlapped in that one site.

        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();

        // prepare partial tree sequence for traceback
        let mut partial_tree_sequence = PartialTreeSequence::with_capacity(3);
        partial_tree_sequence.push(vec![], vec![], true);
        // this node will be compressed into the previous node
        partial_tree_sequence.push(
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![],
            true,
        );
        // this will be the ancestor we traceback, and it will compress into the parent at the last site.
        partial_tree_sequence.push(
            vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(1),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(1),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            // we want it to be uncompressed at the last site, so we can traceback it without weird issues
            vec![VariantIndex::from_usize(9)],
            true,
        );

        edges.insert_sequence_node(
            Ancestor(1),
            &partial_tree_sequence.edges[1],
            &partial_tree_sequence.mutations[1],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &partial_tree_sequence.edges[2],
            &partial_tree_sequence.mutations[2],
        );

        // simulate the compression events
        let mut iterator = ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        );

        iterator.for_each(|(site, tree, nodes)| {
            // insert a mutation on the first site for the root node. This is what we want to see
            // in the traceback despite the direct parent being Ancestor(1) and its interval
            // only stretching to site 0 (so if site 0 is not included in the iterator, we have a
            // off-by-1 error
            if site.0 == 0 {
                // only root should be active
                assert_eq!(nodes.iter().copied().collect::<Vec<_>>(), vec![Ancestor(0)]);

                *tree.likelihood(Ancestor(0)) = 0.5;
                tree.insert_mutation_event(Ancestor(0), site);
            }

            if site.0 == 1 {
                // Ancestor(0) and Ancestor(2) should be decompressed here
                assert_eq!(*nodes, vec![Ancestor(0), Ancestor(2)]);

                // make sure we recompress the ancestor
                *tree.likelihood(Ancestor(2)) = *tree.likelihood(Ancestor(0));
            }
        });

        // check that we see all expected events
        let mut traceback_iterator = iterator
            .traceback(&partial_tree_sequence, Ancestor(2))
            .filter(|e| !matches!(e.kind, ViterbiEventKind::Uncompressed));
        assert_eq!(traceback_iterator.clone().count(), 1);
        assert_eq!(traceback_iterator.next().unwrap().kind, Mutation);
    }

    #[test]
    fn test_queue_order() {
        // test whether the queue order is sensible in upholding invariants.
        // specifically it shouldn't remove nodes before their children have changed parents or have been removed,
        // and shouldnt insert nodes before their parents have been inserted

        let mut ix = ViterbiIterator::new(4, false, 1);
        let mut edges = EdgeSequence::new();

        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(5),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );
        // this must change its parent before the above node is removed from the tree
        // furthermore, this must be inserted after the Ancestor(1) node, otherwise it has no parent
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            &vec![VariantIndex::from_usize(0)],
        );
        // this must be removed before Ancestor(1) is removed, otherwise the queue will crash when it searches the parent
        edges.insert_sequence_node(
            Ancestor(3),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(5),
                Ancestor(1),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // if one of the above restrictions is violated, the queue will panic, so we don't need to explicitely test this
        // behavior (and we also can't because the closure is only called after all events for a site have been processed)
        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            4,
        )
        .for_each(|tuple| {
            black_box(tuple);
        });
    }

    #[test]
    fn test_uncompressed_tree_integrity() {
        // test whether the uncompressed tree array always points to the correct parent
        let mut ix = ViterbiIterator::new(3, false, 1);
        let mut edges = EdgeSequence::new();
        let mut counter = 0;

        edges.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(2),
                VariantIndex::from_usize(4),
            ],
        );
        edges.insert_sequence_node(
            Ancestor(2),
            &vec![
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(0),
                    VariantIndex::from_usize(5),
                    Ancestor(1),
                ),
                PartialSequenceEdge::new(
                    VariantIndex::from_usize(5),
                    VariantIndex::from_usize(10),
                    Ancestor(0),
                ),
            ],
            &vec![
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(1),
                VariantIndex::from_usize(2),
                VariantIndex::from_usize(4),
            ],
        );

        // check tree
        ix.iter_sites(
            &edges,
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
            3,
        )
        .for_each(|(_, tree, _)| {
            // always trigger recompression
            match counter {
                0 => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(1)),
                        Some(Ancestor(0))
                    );
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(2)),
                        Some(Ancestor(1))
                    );
                }
                1 => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                    // ancestor 1 is compressed
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(2)),
                        Some(Ancestor(0))
                    );
                }
                2 => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(1)),
                        Some(Ancestor(0))
                    );
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(2)),
                        Some(Ancestor(1))
                    );
                }
                4 => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(1)),
                        Some(Ancestor(0))
                    );
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(2)),
                        Some(Ancestor(1))
                    );

                    // prevent recompression, to test if Ancestor(2) still recombines to Ancestor(0)
                    *tree.likelihood(Ancestor(1)) = 0.1;
                }
                5 => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(1)),
                        Some(Ancestor(0))
                    );
                    assert_eq!(
                        find_uncompressed_parent(tree, Ancestor(2)),
                        Some(Ancestor(0))
                    );

                    // recompress everything
                    *tree.likelihood(Ancestor(1)) = 0.0;
                }
                _ => {
                    assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                }
            }
            counter += 1usize;
        });
    }
}
