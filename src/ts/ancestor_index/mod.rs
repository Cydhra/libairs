use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::iter::Rev;
use std::slice::Iter;

use crate::ancestors::Ancestor;
use crate::ts::ancestor_index::iter::*;
use crate::ts::partial_sequence::{PartialSequenceEdge, PartialTreeSequence};
use crate::variants::VariantIndex;

pub mod iter;

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

    /// Beginning from here, traceback is impossible from the current ancestor
    Compressed(VariantIndex),
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
        use_recompression_threshold: bool,
        inv_recompression_threshold: usize,
        tree_sequence: &PartialTreeSequence,
    ) -> Self {
        let mut ancestor_index = Self::new(
            max_nodes,
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
        let mut marginal_tree = MarginalTree::new(
            start,
            self.num_nodes,
            limit_nodes,
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

        PartialTreeSequenceIterator::new(marginal_tree, start, site, end, queue)
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
