use crate::ts::ancestor_array::{Ancestor, VariantIndex};

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
