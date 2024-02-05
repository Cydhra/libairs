use crate::ancestors::{Ancestor, AncestorArray, VariantIndex};
use crate::dna::SequencePosition;
use crate::ts::tree_sequence::{TreeSequence, TreeSequenceInterval, TreeSequenceNode};

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

    pub(crate) fn as_tree_sequence(
        &self,
        sequence_length: SequencePosition,
        variant_positions: &[SequencePosition],
        ancestors: AncestorArray,
    ) -> TreeSequence {
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
                                // TODO hide conversion of variant index to sequence position somewhere in ancestors module
                                let start = Self::convert_variant_index(
                                    edge.start(),
                                    sequence_length,
                                    variant_positions,
                                );
                                let end = Self::convert_variant_index(
                                    edge.end(),
                                    sequence_length,
                                    variant_positions,
                                );
                                TreeSequenceInterval::new(parent.0, start, end)
                            })
                            .collect(),
                        // TODO don't convert to usize here
                        self.mutations[idx].iter().map(|index| index.0).collect(),
                    )
                })
                .collect(),
            1: ancestors,
        }
    }

    fn convert_variant_index(
        index: VariantIndex,
        sequence_length: SequencePosition,
        variant_positions: &[SequencePosition],
    ) -> SequencePosition {
        if let 0 = index.0 {
            SequencePosition::from_usize(0)
        } else if index.0 == variant_positions.len() {
            sequence_length
        } else {
            variant_positions[index.0]
        }
    }
}
