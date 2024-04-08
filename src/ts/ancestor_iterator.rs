use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::iter::{Peekable, Rev};
use std::slice::Iter;

use crate::ancestors::Ancestor;
use crate::ts::partial_sequence::{PartialSequenceEdge, PartialTreeSequence};
use crate::variants::VariantIndex;

/// A DNA sequence site visited during the viterbi algorithm.
/// Consists of a [`VariantIndex`] and a mutable reference to the marginal tree at that site.
pub(crate) type Site<'a, 'o> = (VariantIndex, &'a mut MarginalTree<'o>);

/// Events that a single ancestor can experience during the Viterbi algorithm.
/// These events are generated during the forward search and are used to reconstruct
/// the most likely path during backtracing by storing where to look for mutations and recombinations.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) struct ViterbiEvent {
    pub(super) kind: ViterbiEventKind,
    pub(super) site: VariantIndex,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum ViterbiEventKind {
    /// Mutation event here
    Mutation,

    /// Recombination event here
    Recombination,

    /// Beginning from here, traceback is possible from the current ancestor
    Decompress,

    /// Beginning from here, traceback must be continued at the given `Ancestor`
    Copy(Ancestor),
}

/// A helper structure for the Viterbi algorithm that helps to iterate through the sites, updating
/// the marginal tree at the currently visited site and helping with tree compression.
/// The iterator only produces indices into the ancestor array and sequences, it does not hold the
/// actual ancestors.
///
/// Whenever new tree edges are calculated, they can be inserted into the iterator so the tree
/// compression improves for future iterations.
///
/// To insert nodes that have no edges yet, call [`AncestorIndex::insert_free_node`]. This will
/// make the ancestor available during Viterbi without compressing it into the tree.
///
/// To upgrade nodes with edges (enabling tree compression during the Viterbi algorithm), call
/// [`AncestorIndex::insert_sequence_node`].
///
/// The iterator will start at the provided start point and iterate through all sites,
/// updating the marginal tree for each site and then calling a consumer function
/// with the [`VariantIndex`] and the [`MarginalTree`] (see: [`Site`])
#[derive(Clone)]
pub(crate) struct AncestorIndex {
    /// A sorted set containing all tree-sequence events:
    /// - Start of a new tree node
    /// - Change of parent of a tree node
    /// - End of a tree node
    edge_index: BTreeSet<SequenceEvent>,

    /// The number of nodes in the tree sequence.
    /// This number is the maximum amount of nodes that have likelihoods during the Viterbi
    /// algorithm.
    /// The actual iterator can be limited to a smaller number of nodes.
    /// Not all nodes included in this number necessarily have edges in the tree sequence.
    num_nodes: usize,

    // the following are pre-allocated arrays that are used to store the state of the marginal tree
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

    /// Mutations, Recombinations and indirections
    viterbi_events: Vec<Vec<ViterbiEvent>>,

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: Vec<usize>,

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
    inv_recompression_threshold: usize,
}

impl AncestorIndex {
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
        max_nodes: usize,
        variant_count: usize,
        use_recompression_threshold: bool,
        inv_recompression_threshold: usize,
    ) -> Self {
        assert!(
            inv_recompression_threshold > 0,
            "Recompression interval must be greater than 0"
        );
        Self {
            edge_index: BTreeSet::new(),
            num_nodes: 1,
            parents: vec![None; max_nodes],
            children: vec![Vec::new(); max_nodes],
            uncompressed_parents: vec![None; max_nodes],
            is_compressed: vec![true; max_nodes],
            active_nodes: Vec::new(),
            likelihoods: vec![-1.0f64; max_nodes],
            viterbi_events: vec![Vec::new(); max_nodes],
            last_compressed: vec![0; max_nodes],
            use_recompression_threshold,
            inv_recompression_threshold,
        }
    }

    /// Instantiate a new ancestor index and insert up to `max_nodes` nodes from the given tree
    /// sequence.
    ///
    /// # Parameters
    /// - `max_nodes`: The maximum number of nodes to insert into the iterator. Can be lower or higher
    /// than the amount of nodes in the `tree_sequence`.
    /// - `variant_count`: The number of variants in the genome sequence
    /// - `use_recompression_threshold`: Whether to use the recompression threshold to avoid
    /// recompressing when only a few nodes are active. If this is false, the algorithm will
    /// recompress in every iteration, regardless of efficacy.
    /// - `inv_recompression_threshold`: The inverse recompression threshold. The iterator will
    /// attempt to recompress the marginal tree when more than `1/inv_recompression_threshold * num_nodes`
    /// nodes of the marginal tree are active. If `use_recompression_threshold` is false, this
    /// parameter is ignored.
    /// - `tree_sequence`: The tree sequence to insert into the iterator
    pub(crate) fn from_tree_sequence(
        max_nodes: usize,
        variant_count: usize,
        use_recompression_threshold: bool,
        inv_recompression_threshold: usize,
        tree_sequence: &PartialTreeSequence,
    ) -> Self {
        let mut ancestor_index = Self::new(
            max_nodes,
            variant_count,
            use_recompression_threshold,
            inv_recompression_threshold,
        );

        for (node, (edges, mutations)) in tree_sequence
            .edges
            .iter()
            .zip(tree_sequence.mutations.iter())
            .enumerate()
            .skip(1)
        {
            ancestor_index.insert_sequence_node(Ancestor(node), edges, mutations);
        }

        ancestor_index
    }

    /// Insert a free node into the iterator. This node will be available during the Viterbi
    /// algorithm, but it will not be compressed into the tree (since it hasn't been matched against
    /// the tree sequence yet).
    ///
    /// The nodes must be inserted in order of increasing index.
    ///
    /// A free node can be upgraded to a sequence node by calling [`AncestorIndex::insert_sequence_node`].
    /// Upgrading does not need to be done in order.
    ///
    /// # Parameters
    /// - `node`: The tree node that will be inserted into the tree sequence
    /// - `start`: The index of the first variant site of its ancestral sequence. The node will
    /// not participate in the Viterbi algorithm before this site.
    /// - `end`: The index of the last variant site of its ancestral sequence.
    ///
    /// # Panics
    /// Panics if the given tree node index is not the next index in the sequence.
    pub(crate) fn insert_free_node(
        &mut self,
        node: Ancestor,
        start: VariantIndex,
        end: VariantIndex,
    ) {
        assert_eq!(
            node.0, self.num_nodes,
            "Tree nodes must be inserted in order: {} != {}",
            node.0, self.num_nodes
        );
        self.num_nodes += 1;

        // insert start and end of the ancestor
        self.edge_index.insert(SequenceEvent {
            site: start,
            node,
            kind: SequenceEventKind::StartFree,
        });

        self.edge_index.insert(SequenceEvent {
            site: end,
            node,
            kind: SequenceEventKind::End,
        });
    }

    /// Insert a tree sequence node into the partial tree sequence, which will be exploited by the
    /// iterator to compress trees.
    /// Alternatively, if the node already exists and has no edges in the tree sequence
    /// (i.e. is a free node), it will be upgraded into a tree node.
    ///
    /// The nodes must be inserted in order of increasing index. This does not apply to upgrading.
    ///
    /// # Parameters
    /// - `tree_node`: The tree node handle that will be inserted into the tree sequence
    /// - `edges`: The edges that are associated with the tree node
    /// - `mutations`: The mutations that are associated with the tree node
    ///
    /// # Panics
    /// Panics, if the given tree node index is not already a free node and is not the next index
    /// in the sequence.
    pub(crate) fn insert_sequence_node(
        &mut self,
        tree_node: Ancestor,
        edges: &[PartialSequenceEdge],
        mutations: &[VariantIndex],
    ) {
        if tree_node.0 >= self.num_nodes {
            self.num_nodes += 1;
            assert_eq!(
                tree_node.0,
                self.num_nodes - 1,
                "Tree nodes must be inserted in order"
            );

            // insert the end only if this node has not been inserted before as a free node
            self.edge_index.insert(SequenceEvent {
                site: edges.last().unwrap().end(),
                node: tree_node,
                kind: SequenceEventKind::End,
            });
        } else {
            // remove the start event of the inserted node, since we assume it has been inserted before
            let success = self.edge_index.remove(&SequenceEvent {
                site: edges.first().unwrap().start(),
                node: tree_node,
                kind: SequenceEventKind::StartFree,
            });
            assert!(
                success,
                "attempted to upgrade a node that had no free start event associated with it"
            );
        }

        if !edges.is_empty() {
            self.edge_index.insert(SequenceEvent {
                site: edges.first().unwrap().start(),
                node: tree_node,
                kind: SequenceEventKind::Start {
                    parent: edges.first().unwrap().parent(),
                },
            });

            for edge in edges.iter().skip(1) {
                self.edge_index.insert(SequenceEvent {
                    site: edge.start(),
                    node: tree_node,
                    kind: SequenceEventKind::ChangeParent {
                        new_parent: edge.parent(),
                    },
                });
            }
        } else {
            assert_eq!(self.num_nodes, 1, "only the root node can have no edges");
        }

        for &mutation in mutations {
            self.edge_index.insert(SequenceEvent {
                site: mutation,
                node: tree_node,
                kind: SequenceEventKind::Mutation,
            });
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
    pub(crate) fn sites(
        &mut self,
        start: VariantIndex,
        end: VariantIndex,
        limit_nodes: usize,
    ) -> PartialTreeSequenceIterator<impl Iterator<Item = &'_ SequenceEvent>> {
        let current_sequence_length = end.get_variant_distance(start);
        let mut marginal_tree = MarginalTree::new(
            start,
            self.num_nodes,
            limit_nodes,
            current_sequence_length,
            &mut self.parents[0..limit_nodes],
            &mut self.children[0..limit_nodes],
            &mut self.uncompressed_parents[0..limit_nodes],
            &mut self.is_compressed[0..limit_nodes],
            &mut self.active_nodes,
            &mut self.likelihoods[0..limit_nodes],
            &mut self.viterbi_events[0..limit_nodes],
            &mut self.last_compressed[0..limit_nodes],
            self.use_recompression_threshold,
            self.inv_recompression_threshold,
        );

        let site = start;
        let mut queue = self.edge_index.iter().peekable();

        marginal_tree.advance_to_site(&mut queue, site, true, false);

        PartialTreeSequenceIterator {
            marginal_tree,
            start,
            site,
            end,
            queue,
        }
    }
}

/// Borrowed iterator through the partial tree sequence held by an [`AncestorIndex`].
/// Not an actual implementation of [`std::iter::IterMut`], because it borrows the same marginal tree for each site
/// mutably, so we need to make sure the reference cannot escape the iterator.
///
/// This means this struct does not actually implement Iterator traits, but offers a single method
/// [`PartialTreeSequenceIterator::for_each`], which consumers have to call to process the tree sequence.
pub(crate) struct PartialTreeSequenceIterator<'a, 'o, I: Iterator<Item = &'a SequenceEvent>> {
    marginal_tree: MarginalTree<'o>,
    start: VariantIndex,
    site: VariantIndex,
    end: VariantIndex,
    queue: Peekable<I>,
}

impl<'a, 'o, I: Iterator<Item = &'a SequenceEvent>> PartialTreeSequenceIterator<'a, 'o, I> {
    pub(crate) fn for_each<F: FnMut(Site)>(&mut self, mut consumer: F) {
        if self.site == self.end {
            return;
        }

        // first site has special treatment because we don't need to decompress edges that start here
        self.marginal_tree
            .advance_to_site(&mut self.queue, self.site.next(), false, true);
        consumer((self.site, &mut self.marginal_tree));
        self.marginal_tree
            .recompress_tree(self.site.get_variant_distance(self.start));
        self.site = self.site.next();

        while self.site < self.end {
            self.marginal_tree
                .advance_to_site(&mut self.queue, self.site.next(), false, false);

            consumer((self.site, &mut self.marginal_tree));
            if !self.marginal_tree.use_recompression_threshold
                || self.marginal_tree.active_nodes.len()
                    > self.marginal_tree.limit_nodes
                        / self.marginal_tree.inv_recompression_threshold
            {
                self.marginal_tree
                    .recompress_tree(self.site.get_variant_distance(self.start));
            }
            self.site = self.site.next();
        }
    }

    /// Create a traceback iterator through the Viterbi table generated during a previous call to
    /// `for_each`. The iterator will be initialized with `current_ancestor` for the traceback and
    /// automatically retrieve events from its correct ancestor whenever the actual ancestor is
    /// compressed.
    pub(crate) fn traceback<'s>(
        &'s self,
        current_ancestor: Ancestor,
    ) -> TracebackSequenceIterator<'s, 'o>
    where
        's: 'o,
    {
        TracebackSequenceIterator {
            marginal_tree: &self.marginal_tree,
            start: self.start,
            end: self.end,
            current: self.end,
            current_ancestor,
            iter: self.marginal_tree.viterbi_events[current_ancestor.0]
                .iter()
                .rev()
                .peekable(),
            inner: None,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum SequenceEventKind {
    Start { parent: Ancestor },
    StartFree,
    ChangeParent { new_parent: Ancestor },
    Mutation,
    End,
}

impl PartialOrd<Self> for SequenceEventKind {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SequenceEventKind {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Start { parent: p1 }, Self::Start { parent: p2 }) => p1.cmp(p2),
            (Self::Start { .. }, _) => Ordering::Less,
            (_, Self::Start { .. }) => Ordering::Greater,
            (Self::StartFree, Self::StartFree) => Ordering::Equal,
            (Self::StartFree, _) => Ordering::Less,
            (_, Self::StartFree) => Ordering::Greater,
            (Self::ChangeParent { new_parent: p1 }, Self::ChangeParent { new_parent: p2 }) => {
                p1.cmp(p2)
            }
            (Self::ChangeParent { .. }, _) => Ordering::Less,
            (_, Self::ChangeParent { .. }) => Ordering::Greater,
            (Self::Mutation, Self::Mutation) => Ordering::Equal,
            (Self::Mutation, _) => Ordering::Less,
            (_, Self::Mutation) => Ordering::Greater,
            (Self::End, Self::End) => Ordering::Equal,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SequenceEvent {
    site: VariantIndex,
    node: Ancestor,
    kind: SequenceEventKind,
}

impl PartialOrd<Self> for SequenceEvent {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SequenceEvent {
    fn cmp(&self, other: &Self) -> Ordering {
        self.site
            .cmp(&other.site)
            .then(self.kind.cmp(&other.kind))
            .then(match self.kind {
                // we need to delete the tree from the bottom,
                SequenceEventKind::End => self.node.cmp(&other.node).reverse(),
                // but build it from the top
                _ => self.node.cmp(&other.node),
            })
    }
}

/// Iterator to iterate back through the viterbi table generated by the Viterbi algorithm.
/// It holds the marginal tree that was modified by the Viterbi algorithm and iterates back through
/// some aspects of the tree sequence, but only on the most likely path.
///
/// TODO this needs way more description
pub(super) struct TracebackSequenceIterator<'a, 'o>
where
    'o: 'a,
{
    marginal_tree: &'a MarginalTree<'o>,
    start: VariantIndex,
    end: VariantIndex,
    current: VariantIndex,
    current_ancestor: Ancestor,
    iter: Peekable<Rev<Iter<'o, ViterbiEvent>>>,
    inner: Option<Box<TracebackSequenceIterator<'a, 'o>>>,
}

impl<'o> TracebackSequenceIterator<'o, 'o> {
    pub fn switch_ancestor(&mut self, ancestor: Ancestor) {
        self.current_ancestor = ancestor;

        // find the position BEFORE the current site in the new ancestor,
        // since we do not want to copy the current site
        let new_cursor = self.marginal_tree.viterbi_events[ancestor.0]
            .binary_search_by(|e| e.site.cmp(&self.current));
        let pos = match new_cursor {
            Ok(mut pos) => {
                while pos > 0
                    && self.marginal_tree.viterbi_events[ancestor.0][pos - 1].site == self.current
                {
                    pos -= 1;
                }
                pos
            }
            Err(pos) => pos,
        };

        self.iter = self.marginal_tree.viterbi_events[ancestor.0][..pos]
            .iter()
            .rev()
            .peekable();

        // println!("Switching to ancestor {} at site {}", ancestor.0, self.current);
        // println!("Next site is {:?}", self.iter.peek().map(|e| e.site));

        self.inner = None;
    }
}

impl<'o> Iterator for TracebackSequenceIterator<'o, 'o> {
    type Item = &'o ViterbiEvent;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(ref mut inner) = self.inner {
            let element = inner.next();
            if element.is_none() {
                self.inner = None;
            } else {
                self.current = element.unwrap().site;
                return element;
            }
        }

        let element = self.iter.next();
        if let Some(event) = element {
            if event.site < self.start {
                return None;
            } else {
                self.current = event.site;
            }

            match event.kind {
                ViterbiEventKind::Copy(ancestor) => {
                    let start_site = if let Some(ViterbiEvent { site, kind: _ }) = self.iter.peek()
                    {
                        *site
                    } else {
                        self.start
                    };

                    let end_cursor = self.marginal_tree.viterbi_events[ancestor.0]
                        .binary_search_by(|e| e.site.cmp(&event.site));
                    let start_cursor = self.marginal_tree.viterbi_events[ancestor.0]
                        .binary_search_by(|e| e.site.cmp(&start_site));

                    // find the position AFTER the current site in the new ancestor
                    // since we want to also copy the current site
                    let end = match end_cursor {
                        Ok(mut pos) => {
                            while self.marginal_tree.viterbi_events[ancestor.0].len() > pos
                                && self.marginal_tree.viterbi_events[ancestor.0][pos].site
                                    == event.site
                            {
                                pos += 1;
                            }
                            pos
                        }
                        Err(pos) => pos,
                    };

                    // find the position AFTER the start_cursor,
                    // since we do not want to copy that site
                    let start = match start_cursor {
                        Ok(mut pos) => {
                            while self.marginal_tree.viterbi_events[ancestor.0].len() > pos
                                && self.marginal_tree.viterbi_events[ancestor.0][pos].site
                                    == start_site
                            {
                                pos += 1;
                            }
                            pos
                        }
                        Err(pos) => pos,
                    };

                    // println!("Copying from {} to {} ({} to {}) from {}", start, end, start_site, event.site, ancestor.0);
                    self.inner = Some(Box::new(TracebackSequenceIterator {
                        marginal_tree: self.marginal_tree,
                        start: start_site,
                        end: event.site,
                        current: event.site,
                        current_ancestor: ancestor,
                        iter: self.marginal_tree.viterbi_events[ancestor.0][start..end]
                            .iter()
                            .rev()
                            .peekable(),
                        inner: None,
                    }));
                    self.next()
                }
                ViterbiEventKind::Decompress => self.next(),
                _ => element,
            }
        } else {
            None
        }
    }
}

/// An access structure into the tree sequence at a specific site.
/// This structure is generated and by the [`AncestorIterator`] during iteration.
/// It provides method for tree node compression during the viterbi algorithm.
///
/// It further serves as state for the Viterbi algorithm, as it holds the likelihoods of the
/// nodes in the tree, as well as recombination and mutation sites.
///
/// # Free Nodes
/// The structure might include some nodes that aren't part of the tree,
/// if no edges have been calculated for those nodes yet.
/// As a result, those nodes cannot be recompressed (they have no parents).
/// The marginal tree offers methods to iterate through all uncompressed nodes.
///
/// [`AncestorIterator`]: AncestorIndex
#[derive(Debug)]
pub(crate) struct MarginalTree<'o> {
    /// The actual parents of the tree nodes. `None` for free nodes and the root, parents might be
    /// compressed.
    parents: &'o mut [Option<Ancestor>],

    /// Node children
    children: &'o mut [Vec<Ancestor>],

    /// The most recent uncompressed ancestor node of each node in the marginal tree.
    /// `None` for root and free nodes
    uncompressed_parents: &'o mut [Option<Ancestor>],

    /// Whether the node is currently compressed.
    is_compressed: &'o mut [bool],

    /// The list of active nodes, i.e. nodes that are not compressed and are used in the Viterbi
    /// algorithm.
    active_nodes: &'o mut Vec<Ancestor>,

    /// The likelihood of each node in the tree. Likelihood for compressed nodes is wrong, but
    /// they are ignored during the Viterbi algorithm.
    likelihoods: &'o mut [f64],

    /// Mutations, Recombinations and indirections
    viterbi_events: &'o mut [Vec<ViterbiEvent>],

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: &'o mut [usize],

    /// The number of nodes in the tree. Nodes with a higher index are ignored in the tree.
    limit_nodes: usize,

    /// Start index of the current iteration
    start: VariantIndex,

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
    inv_recompression_threshold: usize,
}

impl<'o> MarginalTree<'o> {
    fn new(
        start: VariantIndex,
        num_nodes: usize,
        limit_nodes: usize,
        ancestor_length: usize,
        parents: &'o mut [Option<Ancestor>],
        children: &'o mut [Vec<Ancestor>],
        uncompressed_parents: &'o mut [Option<Ancestor>],
        is_compressed: &'o mut [bool],
        active_nodes: &'o mut Vec<Ancestor>,
        likelihoods: &'o mut [f64],
        viterbi_events: &'o mut [Vec<ViterbiEvent>],
        last_compressed: &'o mut [usize],
        use_recompression_threshold: bool,
        inv_recompression_threshold: usize,
    ) -> Self {
        // println!("========================================================");
        // println!("Creating new marginal tree: {} nodes, {} limit, {} start", num_nodes, limit_nodes, start);
        debug_assert!(num_nodes > 0, "Tree must have at least one node");

        // re-initialize vectors into the default state where needed
        active_nodes.clear();
        viterbi_events.iter_mut().for_each(|i| i.clear());
        last_compressed.fill(0);

        // the other states are updated whenever nodes are added to the marginal tree

        let mut marginal_tree = Self {
            parents,
            children,
            uncompressed_parents,
            is_compressed,
            active_nodes,
            likelihoods,
            viterbi_events,
            last_compressed,
            limit_nodes,
            start,
            use_recompression_threshold,
            inv_recompression_threshold,
        };

        marginal_tree.add_initial_node(Ancestor(0), false);
        assert!(
            limit_nodes <= num_nodes,
            "Limit nodes must be less than or equal to num nodes"
        );

        marginal_tree
    }

    /// An iterator through the currently uncompressed nodes
    pub(crate) fn nodes(&self) -> impl Iterator<Item = Ancestor> {
        // TODO can we maybe get away without cloning here?
        self.active_nodes.clone().into_iter()
    }

    /// Get the number of valid nodes in the tree.
    /// This includes nodes that cannot be copied from at the current location, but might
    /// influence nodes that can be copied from.
    /// This also includes free nodes that aren't inserted into the tree yet.
    /// This value is interesting mostly for the Markov Chain during the Viterbi algorithm.
    pub(crate) fn num_nodes(&self) -> usize {
        self.likelihoods.len()
    }

    /// Get the likelihood of the given node.
    pub(crate) fn likelihood(&mut self, node: Ancestor) -> &mut f64 {
        debug_assert!(!self.is_compressed[node.0]);
        &mut self.likelihoods[node.0]
    }

    /// Insert a recombination for the given ancestor at the given site (absolute site, not relative
    /// to the candidate)
    pub(crate) fn insert_recombination_event(&mut self, node: Ancestor, site: VariantIndex) {
        self.viterbi_events[node.0].push(ViterbiEvent {
            kind: ViterbiEventKind::Recombination,
            site,
        });
    }

    /// Insert a mutation for the given ancestor at the given site (absolute site, not relative to the
    /// candidate)
    pub(crate) fn insert_mutation_event(&mut self, node: Ancestor, site: VariantIndex) {
        self.viterbi_events[node.0].push(ViterbiEvent {
            kind: ViterbiEventKind::Mutation,
            site,
        });
    }

    /// Return whether the node is currently compressed.
    fn is_compressed(&self, node: Ancestor) -> bool {
        self.is_compressed[node.0]
    }

    /// Compress all active nodes that have the same likelihood as their parent, so they can be
    /// ignored during the Viterbi algorithm.
    ///
    /// # Parameters
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    pub(crate) fn recompress_tree(&mut self, site_index: usize) {
        self.active_nodes.retain(|&node| {
            debug_assert!(self.is_compressed[node.0] == false);

            if let Some(parent) = self.uncompressed_parents[node.0] {
                if self.likelihoods[node.0] == self.likelihoods[parent.0] {
                    self.last_compressed[node.0] = site_index + 1;
                    self.is_compressed[node.0] = true;
                    // insert event that informs traceback that starting from here we are uncompressed
                    self.viterbi_events[node.0].push(ViterbiEvent {
                        kind: ViterbiEventKind::Decompress,
                        site: self.start + site_index,
                    });

                    let mut queue = self.children[node.0].clone();
                    while let Some(child) = queue.pop() {
                        if self.uncompressed_parents[child.0] == Some(node) {
                            self.uncompressed_parents[child.0] = Some(parent);

                            if self.is_compressed[child.0] {
                                // insert event that informs traceback that starting from here we are compressed into node instead into the new parent
                                self.viterbi_events[child.0].push(ViterbiEvent {
                                    kind: ViterbiEventKind::Copy(node),
                                    site: self.start + site_index,
                                });

                                queue.extend(self.children[child.0].iter().copied());
                            }
                        }
                    }
                    false
                } else {
                    true
                }
            } else {
                true
            }
        });
    }

    /// Search the parent of the given node and update the uncompressed_tree_parents array, the is_compressed array
    /// and add the node to active nodes.
    /// Returns the parent.
    fn set_uncompressed(&mut self, node: Ancestor) -> Ancestor {
        let parent = self.uncompressed_parents[node.0];
        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        parent.unwrap()
    }

    /// Ensures a node is decompressed, i.e. it is in the active node list and if it is not,
    /// it is added to the list and its likelihood is set to the likelihood of its parent.
    /// Recombination and mutation sites are copied from the parent.
    /// If you want to decompress a node but not inherit from a parent, use [`set_uncompressed`]
    ///
    /// If the node is already decompressed, this method does nothing.
    ///
    /// # Parameters
    /// - `node`: The node to ensure is decompressed
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    ///
    /// [`set_uncompressed`]: MarginalTree::set_uncompressed
    fn ensure_decompressed(&mut self, node: Ancestor, site_index: usize) {
        if self.is_compressed(node) {
            debug_assert!(!self.active_nodes.contains(&node));

            // search the first uncompressed ancestor node
            let uncompressed_parent = self.set_uncompressed(node);

            self.likelihoods[node.0] = self.likelihoods[uncompressed_parent.0];

            // TODO we dont need last_compressed anymore, since we record compression events anyway?
            // copy the recombination and mutation sites from the parent
            let last_compressed_begin = self.last_compressed[node.0];

            // record event for traceback that starting from here we are compressed into the parent
            if site_index > 0 {
                self.viterbi_events[node.0].push(ViterbiEvent {
                    kind: ViterbiEventKind::Copy(uncompressed_parent),
                    site: self.start + site_index - 1,
                });
            }
        }
    }

    /// Insert a node into the tree that has no parent, i.e. it is a free node or is the root node.
    /// This method is used to initialize the tree sequence and to insert unmatched nodes.
    ///
    /// # Parameters
    /// - `node` the free node to insert
    /// - `during_viterbi` must be true, if the viterbi algorithm is already running
    fn add_initial_node(&mut self, node: Ancestor, during_viterbi: bool) {
        debug_assert!(self.parents[node.0].is_none());

        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        self.likelihoods[node.0] = if during_viterbi { 0.0 } else { 1.0 };
        self.parents[node.0] = None;
        self.children[node.0].clear();
        self.uncompressed_parents[node.0] = None;
    }

    /// Insert a new node into the tree. This is done whenever a new node begins, it is not the
    /// correct method if a node is decompressed because of a change of its parent or a mutation.
    ///
    /// # Parameters
    /// - `parent`: The parent of the new node
    /// - `child`: The new node to insert
    /// - `keep_compressed`: Whether the node should be kept compressed, i.e. it should not be
    /// added to the set of active nodes.
    fn insert_new_node(&mut self, parent: Ancestor, child: Ancestor, keep_compressed: bool) {
        self.parents[child.0] = Some(parent);
        self.children[parent.0].push(child);

        self.children[child.0].clear();
        self.is_compressed[child.0] = true;

        if !keep_compressed {
            // update the compressed parent. If the parent is compressed, its compressed parent
            // will be used, otherwise the parent itself will be used.
            if self.is_compressed(parent) {
                debug_assert!(
                    self.uncompressed_parents[parent.0].is_some(),
                    "{} has no uncompressed parent",
                    parent.0
                );
                self.uncompressed_parents[child.0] = self.uncompressed_parents[parent.0];
            } else {
                self.uncompressed_parents[child.0] = Some(parent);
            }

            self.set_uncompressed(child);
            self.likelihoods[child.0] = 0.0;
        } else {
            // initially, everything is compressed into the root node
            self.uncompressed_parents[child.0] = Some(Ancestor(0));
        }
    }

    /// Update the tree to reflect a change of parent for a node.
    ///
    /// Whenever a node changes parent, it must be decompressed, as its likelihood will now diverge
    /// from its parent (unless the tree is not yet being used for matching, but just prepared).
    /// This method ensures the node is decompressed and updates the tree to reflect the new parent.
    ///
    /// # Parameters
    /// - `node`: The node that has a new parent
    /// - `new_parent`: The new parent of the node
    /// - `keep_compressed`: Whether the node should be kept compressed, i.e. it should not be
    /// decompressed. This is useful if the tree is not yet being used for matching, but just being
    /// prepared (i.e. the Viterbi algorithm is not yet running, but the ancestor iterator updates
    /// the marginal tree by updating it for every event in the queue before the start-site of
    /// the current candidate ancestor).
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    fn update_node_parent(
        &mut self,
        node: Ancestor,
        new_parent: Ancestor,
        keep_compressed: bool,
        site_index: usize,
    ) {
        let mut old_parent = None;
        if !keep_compressed {
            old_parent = self.uncompressed_parents[node.0];
            self.ensure_decompressed(node, site_index);
        }

        // update parent afterwards to ensure it copies from the correct ancestor
        self.children[self.parents[node.0].unwrap().0].retain(|&n| n != node);
        self.parents[node.0] = Some(new_parent);
        self.children[new_parent.0].push(node);

        if !keep_compressed {
            // update the compressed parent. If the parent is compressed, its compressed parent
            // will be used, otherwise the parent itself will be used.
            if self.is_compressed(new_parent) {
                debug_assert!(self.uncompressed_parents[new_parent.0].is_some());
                self.uncompressed_parents[node.0] = self.uncompressed_parents[new_parent.0];
            } else {
                self.uncompressed_parents[node.0] = Some(new_parent);
            }

            self.update_subtree(node, old_parent, Some(node), site_index)
        }
    }

    /// Whenever a node has a mutation from its parent, it must be decompressed, as its likelihood
    /// will no diverge from its parent.
    /// This method ensures the node is decompressed.
    /// This method must not be called when the tree is just being prepared, as the node will be
    /// decompressed, which is unnecessary outside of the Viterbi algorithm.
    ///
    /// # Parameters
    /// - `node`: The node that has a mutation
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    fn update_node_mutation(&mut self, node: Ancestor, site_index: usize) {
        if self.is_compressed(node) {
            self.ensure_decompressed(node, site_index);
            self.update_subtree(
                node,
                self.uncompressed_parents[node.0],
                Some(node),
                site_index,
            );
        }
    }

    /// Update a subtree to reflect a change of compression of its root.
    ///
    /// # Parameters
    /// - `subtree_root`: The root of the subtree to update
    /// - `old_parent`: The previous uncompressed parent for the subtree. The parameter is an option
    /// for implementation convenience, but cannot be `None`
    /// - `new_parent` The new uncompressed parent for the subtree. The parameter is an option for
    /// implementation convenience, but cannot be `None`
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    fn update_subtree(
        &mut self,
        subtree_root: Ancestor,
        old_parent: Option<Ancestor>,
        new_parent: Option<Ancestor>,
        site_index: usize,
    ) {
        debug_assert!(old_parent.is_some());
        debug_assert!(new_parent.is_some());

        let mut queue = self.children[subtree_root.0].clone();
        while let Some(node) = queue.pop() {
            if self.uncompressed_parents[node.0] == old_parent {
                self.uncompressed_parents[node.0] = new_parent;

                if self.is_compressed(node) && site_index > 0 {
                    // record event for traceback that starting from here we are compressed into the old parent
                    self.viterbi_events[node.0].push(ViterbiEvent {
                        kind: ViterbiEventKind::Copy(old_parent.unwrap()),
                        site: self.start + site_index - 1,
                    });

                    queue.extend(self.children[node.0].iter().copied());
                }
            }
        }
    }

    /// Remove a node from the tree. This is done whenever an ancestor ends.
    fn remove_node(&mut self, node: Ancestor) {
        debug_assert!(
            !self.parents.iter().any(|&n| n == Some(node)),
            "Node {} cannot be removed while it has children",
            node.0
        );

        // if node has no parent, it is a free node and no updates are needed
        if self.parents[node.0].is_some() {
            self.children[self.parents[node.0].unwrap().0].retain(|&n| n != node);
            self.parents[node.0] = None;
        }

        self.is_compressed[node.0] = true;
        self.active_nodes.retain(|&n| n != node);
    }

    /// Advance the tree to the given site by updating the tree to reflect all events in the queue
    /// up to (but excluding) the given site.
    ///
    /// # Parameters
    /// - `event_queue`: The queue of events to process, it must be peekable
    /// - `site`: The site to advance to
    /// - `keep_compressed`: Whether the nodes should be kept compressed, i.e. they should not be
    /// added to the set of active nodes. This is useful if the tree is not yet being used for
    /// matching, but just prepared (i.e. the Viterbi algorithm is not yet running, but the
    /// ancestor iterator updates the marginal tree by updating it for every event in the queue
    /// before the start-site of the current candidate ancestor).
    /// - `mutations_only`: Whether only mutations should be processed. This is useful to prepare
    /// the first site of a candidate ancestor, as all nodes that dont have mutations on that site
    /// have the same likelihood anyway, so they dont need to be decompressed
    fn advance_to_site<'b, I: Iterator<Item = &'b SequenceEvent>>(
        &mut self,
        event_queue: &mut Peekable<I>,
        site: VariantIndex,
        keep_compressed: bool,
        mutations_only: bool,
    ) {
        while event_queue.peek().is_some() && event_queue.peek().unwrap().site < site {
            let event = event_queue.next().unwrap();

            // TODO this should be changed into the absolute index, since we no longer need the
            //  relative index for anything
            let site_index = self.start.get_variant_distance(event.site);

            if event.node.0 >= self.limit_nodes {
                continue;
            }

            match event {
                SequenceEvent {
                    site: _,
                    node,
                    kind: SequenceEventKind::Mutation,
                } => {
                    if !keep_compressed {
                        self.update_node_mutation(*node, site_index);
                    }
                }
                SequenceEvent {
                    site: _,
                    node,
                    kind: SequenceEventKind::Start { parent },
                } => {
                    self.insert_new_node(*parent, *node, keep_compressed || mutations_only);
                }
                SequenceEvent {
                    site: _,
                    node,
                    kind: SequenceEventKind::StartFree,
                } => {
                    self.add_initial_node(*node, !(keep_compressed || mutations_only));
                }
                SequenceEvent {
                    site: _,
                    node,
                    kind: SequenceEventKind::ChangeParent { new_parent },
                } => {
                    self.update_node_parent(
                        *node,
                        *new_parent,
                        keep_compressed || mutations_only,
                        site_index,
                    );
                }
                SequenceEvent {
                    site: _,
                    node,
                    kind: SequenceEventKind::End,
                } => {
                    self.remove_node(*node);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::ts::ancestor_iterator::ViterbiEventKind::{Mutation, Recombination};
    use std::hint::black_box;

    use super::*;

    fn find_uncompressed_parent(tree: &MarginalTree, node: Ancestor) -> Option<Ancestor> {
        let mut parent = tree.parents[node.0];

        while let Some(p) = parent {
            if !tree.is_compressed(p) {
                break;
            }
            parent = tree.parents[p.0];
        }

        parent
    }

    fn find_mutation_site(tree: &MarginalTree, site: usize, node: Ancestor) -> bool {
        tree.viterbi_events[node.0]
            .iter()
            .find(|e| e.site.0 == site && e.kind == ViterbiEventKind::Mutation)
            .is_some()
    }

    #[test]
    fn test_empty_range() {
        // test whether iterating over an empty range doesn't crash and calls the closure zero times
        let mut ix = AncestorIndex::new(2, 1, false, 1);

        ix.insert_free_node(
            Ancestor(1),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );

        ix.sites(
            VariantIndex::from_usize(10),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|_| panic!("closure called for empty range"));
    }

    #[test]
    fn test_ancestor_iteration() {
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        ix.insert_free_node(
            Ancestor(1),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );
        ix.insert_free_node(
            Ancestor(2),
            VariantIndex::from_usize(0),
            VariantIndex::from_usize(10),
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
                assert_eq!(site, VariantIndex::from_usize(counter));
                assert_eq!(tree.num_nodes(), 2);
                assert_eq!(tree.nodes().count(), 2);
                assert_eq!(find_uncompressed_parent(tree, Ancestor(0)), None);
                assert_eq!(find_uncompressed_parent(tree, Ancestor(1)), None);
                counter += 1;
            });
    }

    #[test]
    fn test_simple_tree_compression() {
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        // insert edge from first to root node
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(5)],
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(site, VariantIndex::from_usize(counter));
                assert_eq!(tree.num_nodes(), 2);
                assert_eq!(
                    tree.nodes().count(),
                    if counter < 5 { 1 } else { 2 },
                    "wrong number of nodes at site {}",
                    counter
                );
                counter += 1;
            });
    }

    #[test]
    fn test_simple_tree_recompression() {
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        // insert edge from first to root node
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
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
                    tree.nodes().count(),
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
        let mut ix = AncestorIndex::new(3, 9, false, 1);
        let mut counter = 0;

        // insert edges for two nodes
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(9)],
        );
        ix.insert_sequence_node(
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

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(9), 3)
            .for_each(|(site, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(site, VariantIndex::from_usize(counter));
                assert_eq!(tree.num_nodes(), 3);
                assert_eq!(
                    tree.nodes().count(),
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
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            // diverge on site 1, then again on site 5, recompress in between
            &vec![VariantIndex::from_usize(1), VariantIndex::from_usize(5)],
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(_, tree)| {
                // likelihoods stay at 0, so recompression happens immediately
                assert_eq!(
                    tree.nodes().count(),
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
        let mut ix = AncestorIndex::new(3, 10, false, 1);
        let mut counter = 0;

        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );
        ix.insert_sequence_node(
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

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3)
            .for_each(|(_, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(
                    tree.nodes().count(),
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
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        // insert two nodes, one will be ignored
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // check the tree size is always 2
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(_, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(
                    tree.nodes().count(),
                    2,
                    "wrong number of nodes at site {}",
                    counter
                );
                counter += 1usize;
            });
    }

    #[test]
    fn test_subset_iterator() {
        // test whether an incomplete range of sites iterated yields correct trees
        let mut ix = AncestorIndex::new(3, 10, false, 1);
        let mut counter = 2;

        // insert two nodes
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(5)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(1)],
        );

        // check the tree
        ix.sites(VariantIndex::from_usize(2), VariantIndex::from_usize(10), 3)
            .for_each(|(_, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(
                    tree.nodes().count(),
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
        let mut ix = AncestorIndex::new(2, 10, false, 1);
        let mut counter = 0;

        // insert incomplete ancestor
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(5),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // check the tree
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(_, tree)| {
                // fake likelihoods to prevent recompression
                for i in 0..tree.num_nodes() {
                    tree.likelihoods[i] = 0.1 * i as f64
                }

                assert_eq!(
                    tree.nodes().count(),
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

        let mut ix = AncestorIndex::new(3, 10, false, 1);
        let mut counter = 0;

        // insert a child node that serves as a second node
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // insert another child node that copies from both the root and the child node
        ix.insert_sequence_node(
            Ancestor(2),
            &vec![
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
            &vec![VariantIndex::from_usize(4), VariantIndex::from_usize(9)],
        );

        // add mutations and recombinations with the root and first child node, and check whether the second child node
        // copies them
        let mut iterator = ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3);
        iterator.for_each(|(_, tree)| {
            match counter {
                // don't recompress first two nodes
                0 => tree
                    .nodes()
                    .for_each(|n| tree.likelihoods[n.0] = n.0 as f64 * 0.1),

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

        iterator.marginal_tree.viterbi_events[2]
            .iter()
            .filter(|state| {
                if let Mutation | Recombination = state.kind {
                    true
                } else {
                    false
                }
            })
            .for_each(|state| match state.site.0 {
                2 => assert_eq!(
                    state.kind, Mutation,
                    "expected mutation at {}",
                    state.site.0
                ),
                6 => assert_eq!(
                    state.kind, Mutation,
                    "expected mutation at {}",
                    state.site.0
                ),
                3 => assert_eq!(
                    state.kind, Recombination,
                    "expected recombination at {}",
                    state.site.0
                ),
                7 => assert_eq!(
                    state.kind, Recombination,
                    "expected recombination at {}",
                    state.site.0
                ),
                _ => assert_ne!(
                    state.kind, Recombination,
                    "expected no event at site {}",
                    state.site.0
                ),
            });
    }

    #[test]
    fn test_copying_from_parent_on_change_parent() {
        // test whether the iterator copies the recombination and mutation sites from the parent when a child is
        // decompressed

        let mut ix = AncestorIndex::new(3, 10, false, 1);
        let mut counter = 0;

        // insert a child node that serves as a second node
        ix.insert_sequence_node(
            Ancestor(1),
            &vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            &vec![VariantIndex::from_usize(0)],
        );

        // insert another child node that copies from both the root and the child node, but has no mutations
        ix.insert_sequence_node(
            Ancestor(2),
            &vec![
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
            &vec![],
        );

        // add mutations and recombinations with the root, and check whether the second child node copies them
        let mut iterator = ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3);
        iterator.for_each(|(_, tree)| {
            match counter {
                // don't recompress first two nodes
                0 => tree
                    .nodes()
                    .for_each(|n| tree.likelihoods[n.0] = n.0 as f64 * 0.1),

                // insert mutations and recombination for nodes 0 and 1, and later check if 2 inherited them
                2 => tree.insert_mutation_event(Ancestor(0), VariantIndex::from_usize(2)),
                3 => tree.insert_recombination_event(Ancestor(0), VariantIndex::from_usize(3)),
                5 => tree.insert_mutation_event(Ancestor(1), VariantIndex::from_usize(5)), // this should not be copied because the second ancestor is uncompressed at this site

                _ => {}
            }

            counter += 1;
        });

        iterator.marginal_tree.viterbi_events[2]
            .iter()
            .filter(|state| {
                if let Mutation | Recombination = state.kind {
                    true
                } else {
                    false
                }
            })
            .for_each(|state| match state.site.0 {
                2 => assert_eq!(
                    state.kind, Mutation,
                    "expected mutation at {}",
                    state.site.0
                ),
                3 => assert_eq!(
                    state.kind, Recombination,
                    "expected recombination at {}",
                    state.site.0
                ),
                _ => assert_ne!(
                    state.kind, Recombination,
                    "expected no events at site {}",
                    state.site.0
                ),
            });
    }

    #[test]
    fn test_queue_order() {
        // test whether the queue order is sensible in upholding invariants.
        // specifically it shouldn't remove nodes before their children have changed parents or have been removed,
        // and shouldnt insert nodes before their parents have been inserted

        let mut ix = AncestorIndex::new(4, 10, false, 1);

        ix.insert_sequence_node(
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
        ix.insert_sequence_node(
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
        ix.insert_sequence_node(
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
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 4)
            .for_each(|tuple| {
                black_box(tuple);
            });
    }

    #[test]
    fn test_uncompressed_tree_integrity() {
        // test whether the uncompressed tree array always points to the correct parent
        let mut ix = AncestorIndex::new(3, 10, false, 1);
        let mut counter = 0;

        ix.insert_sequence_node(
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
        ix.insert_sequence_node(
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
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3)
            .for_each(|(_, tree)| {
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
