use crate::ancestors::{Ancestor, AncestorArray, VariantIndex};
use crate::ts::tree_sequence::{TreeSequence, TreeSequenceEdge, TreeSequenceNode};

/// Internal tree sequence edge used during the Viterbi algorithm.
/// Its interval has not been converted into sequence positions yet, conveniently allowing use
/// during the algorithm for the tree compression.
/// The edge is associated with the child node, hence only the parent is stored.
#[derive(Clone, Eq, PartialEq, Debug)]
pub(crate) struct PartialSequenceEdge {
    start: VariantIndex,
    end: VariantIndex,
    parent: Ancestor,
}

impl PartialSequenceEdge {
    /// Create a new edge for the partial tree sequence
    pub(crate) fn new(start: VariantIndex, end: VariantIndex, parent: Ancestor) -> Self {
        Self { start, end, parent }
    }

    /// Get the (inclusive) start index of the edge
    pub(crate) fn start(&self) -> VariantIndex {
        self.start
    }

    /// Get the (exclusive) end index of the edge
    pub(crate) fn end(&self) -> VariantIndex {
        self.end
    }

    /// Get the parent ancestor index of the edge.
    pub(crate) fn parent(&self) -> Ancestor {
        self.parent
    }
}

/// Newtype over a collection of tree sequence edges with variant index intervals.
/// Can be converted into a tree sequence.
pub(crate) struct PartialTreeSequence {
    edges: Vec<Vec<PartialSequenceEdge>>,
    mutations: Vec<Vec<VariantIndex>>,
}

impl PartialTreeSequence {
    /// Create a new partial tree sequence with a given node capacity. Does not set the edge capacity
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            edges: Vec::with_capacity(capacity),
            mutations: Vec::with_capacity(capacity),
        }
    }

    /// Push a new set of edges to the partial tree sequence. The edges are associated with the ancestor that shares
    /// its index with this set of edges.
    pub(crate) fn push(&mut self, edges: Vec<PartialSequenceEdge>, mutations: Vec<VariantIndex>) {
        self.edges.push(edges);
        self.mutations.push(mutations);
    }

    pub(crate) fn as_tree_sequence(&self, ancestors: AncestorArray) -> TreeSequence {
        TreeSequence {
            0: self
                .edges
                .iter()
                .enumerate()
                .map(|(idx, edges)| {
                    TreeSequenceNode::new(
                        idx,
                        edges
                            .iter()
                            .map(|edge| {
                                let parent = edge.parent();
                                let start = ancestors.variant_index_to_sequence_pos(edge.start());
                                let end = ancestors.variant_index_to_sequence_pos(edge.end());
                                TreeSequenceEdge::new(parent.0, start, end)
                            })
                            .collect(),
                        &self.mutations[idx],
                    )
                })
                .collect(),
            1: ancestors,
        }
    }
}
