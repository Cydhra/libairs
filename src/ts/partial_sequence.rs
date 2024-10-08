use crate::ancestors::{Ancestor, AncestorArray};
use crate::ts::tree_sequence::{TreeSequence, TreeSequenceEdge, TreeSequenceNode};
use crate::ts::Mutation;
use crate::variants::VariantIndex;

/// Internal tree sequence edge used during the Viterbi algorithm.
/// Its interval has not been converted into sequence positions yet, conveniently allowing use
/// during the algorithm for the tree compression.
/// The edge is associated with the child node, hence only the parent is stored.
#[derive(Clone, Eq, PartialEq, Debug, serde::Serialize, serde::Deserialize)]
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
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct PartialTreeSequence {
    pub(super) edges: Vec<Vec<PartialSequenceEdge>>,
    pub(super) mutations: Vec<Vec<VariantIndex>>,
    inner_nodes: usize,
}

impl PartialTreeSequence {
    /// Create a new partial tree sequence with a given node capacity. Does not set the edge capacity
    pub(crate) fn with_capacity(capacity: usize) -> Self {
        Self {
            edges: Vec::with_capacity(capacity),
            mutations: Vec::with_capacity(capacity),
            inner_nodes: 0,
        }
    }

    /// Push a new set of edges to the partial tree sequence. The edges are associated with the ancestor that shares
    /// its index with this set of edges.
    pub(crate) fn push(
        &mut self,
        edges: Vec<PartialSequenceEdge>,
        mutations: Vec<VariantIndex>,
        inner_node: bool,
    ) {
        debug_assert!(
            !inner_node || self.edges.len() == self.inner_nodes,
            "cannot append inner nodes after samples have already been appended"
        );
        if inner_node {
            self.inner_nodes += 1;
        }
        self.edges.push(edges);
        self.mutations.push(mutations);
    }

    pub(crate) fn as_tree_sequence(&self, ancestors: &AncestorArray) -> TreeSequence {
        TreeSequence {
            nodes: self
                .edges
                .iter()
                .enumerate()
                .map(|(idx, edges)| {
                    TreeSequenceNode::new(
                        idx as u32,
                        edges
                            .iter()
                            .map(|edge| {
                                let parent = edge.parent();
                                let start = ancestors.variant_index_to_sequence_pos(edge.start());
                                let end = ancestors.variant_index_to_sequence_pos(edge.end());
                                TreeSequenceEdge::new(parent.0, start, end)
                            })
                            .collect(),
                        self.mutations[idx]
                            .iter()
                            .map(|v| Mutation {
                                variant_index: v.unwrap(),
                                derived_state: ancestors.get_derived_state(*v),
                            })
                            .collect::<Vec<_>>(),
                        self.inner_nodes <= idx,
                    )
                })
                .collect(),
            ancestors: ancestors.clone(),
        }
    }
}
