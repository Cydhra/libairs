use crate::ancestors::Ancestor;
use crate::ts::partial_sequence::{PartialSequenceEdge, PartialTreeSequence};
use crate::variants::VariantIndex;
use indexset::BTreeSet;
use std::cmp::Ordering;

/// The kind of event that occurs while iterating a partial tree sequence.
/// The events are ordered by the site they occur at, the kind of event, and the node they are
/// associated with.
/// The events are used by the [`iterator`] to update a [`MarginalTree`] during the Viterbi
/// algorithm.
///
/// [`iterator`]: super::TreeSequenceState
/// [`MarginalTree`]: super::tree::MarginalTree
#[derive(Clone, Debug, Eq, PartialEq)]
pub(super) enum SequenceEventKind {
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
            // TODO: ChangeParent should logically come after Start, since new nodes must be present
            //  before they are being switched to. Having ChangeParent before Start will fail,
            //  if an edge starts right at the beginning of an incomplete ancestor. However, we
            //  achieve tsinfer fidelity by changing the order those events are processed
            (Self::ChangeParent { .. }, Self::ChangeParent { .. }) => Ordering::Equal,
            (Self::ChangeParent { .. }, _) => Ordering::Less,
            (_, Self::ChangeParent { .. }) => Ordering::Greater,
            (Self::Start { parent: p1 }, Self::Start { parent: p2 }) => p1.cmp(p2),
            (Self::Start { .. }, _) => Ordering::Less,
            (_, Self::Start { .. }) => Ordering::Greater,
            (Self::StartFree, Self::StartFree) => Ordering::Equal,
            (Self::StartFree, _) => Ordering::Less,
            (_, Self::StartFree) => Ordering::Greater,

            (Self::Mutation, Self::Mutation) => Ordering::Equal,
            (Self::Mutation, _) => Ordering::Less,
            (_, Self::Mutation) => Ordering::Greater,
            (Self::End, Self::End) => Ordering::Equal,
        }
    }
}

/// Any kind of event in a partial tree sequence that influences the Viterbi algorithm.
/// The events are ordered by the site they occur at, the kind of event, and the node they are
/// associated with.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(in crate::ts) struct SequenceEvent {
    pub(super) site: VariantIndex,
    pub(super) node: Ancestor,
    pub(super) kind: SequenceEventKind,
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

/// A sorted set of all tree-sequence events (edges and mutations), that a [`ViterbiIterator`]
/// iterates over. Since this is common to all threads, it is not part of the iterator, which is
/// thread-local.
pub(in crate::ts) struct EdgeSequence {
    /// A sorted set containing all tree-sequence events:
    /// - Start of a new tree node
    /// - Change of parent of a tree node
    /// - End of a tree node
    pub(super) edge_index: BTreeSet<SequenceEvent>,

    /// The number of nodes in the tree sequence.
    /// This number is the maximum amount of nodes that have likelihoods during the Viterbi
    /// algorithm.
    /// The actual iterator can be limited to a smaller number of nodes.
    /// Not all nodes included in this number necessarily have edges in the tree sequence.
    pub(super) num_nodes: u32,
}

impl EdgeSequence {
    pub(in crate::ts) fn new() -> Self {
        Self {
            edge_index: BTreeSet::new(),
            num_nodes: 1,
        }
    }

    /// Instantiate a new ancestor index and insert up to `max_nodes` nodes from the given tree
    /// sequence.
    ///
    /// # Parameters
    /// - `tree_sequence`: The tree sequence to insert into the iterator
    pub(crate) fn from_tree_sequence(tree_sequence: &PartialTreeSequence) -> Self {
        let mut edge_sequence = Self::new();

        for (node, (edges, mutations)) in tree_sequence
            .edges
            .iter()
            .zip(tree_sequence.mutations.iter())
            .enumerate()
            .skip(1)
        {
            edge_sequence.insert_sequence_node(Ancestor(node as u32), edges, mutations);
        }

        edge_sequence
    }

    /// Insert a free node into the iterator. This node will be available during the Viterbi
    /// algorithm, but it will not be compressed into the tree (since it hasn't been matched against
    /// the tree sequence yet).
    ///
    /// The nodes must be inserted in order of increasing index.
    ///
    /// A free node can be upgraded to a sequence node by calling [`ViterbiIterator::insert_sequence_node`].
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
}
