/// An interval in an ancestor that is covered by a parent node in the tree sequence.
/// The interval is defined by the start and (exclusive) end position of the interval and the index of the
/// parent node.
#[derive(Debug, Clone)]
pub struct TreeSequenceInterval {
    pub(crate) parent: usize,
    pub(crate) start: usize,
    pub(crate) end: usize,
}

impl TreeSequenceInterval {
    pub fn new(parent: usize, start: usize, end: usize) -> Self {
        Self { parent, start, end }
    }
}

/// A node in the tree sequence. The node is defined by the index of the ancestor sequence it
/// represents and a list of intervals that define what parent nodes cover the ancestor sequence.
#[derive(Debug, Clone)]
pub struct TreeSequenceNode {
    // todo hide fields
    pub(crate) ancestor_index: usize,
    pub(crate) node_intervals: Vec<TreeSequenceInterval>,
}

impl TreeSequenceNode {
    pub fn new(ancestor_index: usize, node_intervals: Vec<TreeSequenceInterval>) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals,
        }
    }

    pub fn empty(ancestor_index: usize) -> Self {
        TreeSequenceNode {
            ancestor_index,
            node_intervals: Vec::new(),
        }
    }
}