use indexset::BTreeSet;

use crate::ts::ancestor_array::{Ancestor, VariantIndex};
use crate::ts::partial_sequence::PartialSequenceEdge;

/// A helper structure for the Viterbi algorithm that helps iterating through the sites, updating
/// the marginal tree at the currently visited site and helping with tree compression.
/// The iterator only produces indices into the ancestor sequence structure, it does not hold the
/// actual ancestors.
///
/// Whenever new tree edges are calculated, they can be inserted into the iterator so the tree
/// compression improves for future iterations.
pub(crate) struct AncestorIterator {
    edge_index: BTreeSet<usize>,
}

impl AncestorIterator {
    /// Create a new empty ancestor iterator
    pub(crate) fn new() -> Self {
        Self {
            edge_index: BTreeSet::new(),
        }
    }

    /// Insert tree sequence edges into the partial tree sequence, which will be exploited by the
    /// iterator to compress trees.
    pub(crate) fn insert_sequence_edge(&mut self, edge: PartialSequenceEdge) {
        todo!()
    }

    /// Iterate through the variant sites of the ancestral tree sequence.
    /// The iterator will start at the provided start point and iterate through all sites
    /// until the given end index (exclusive).
    /// It provides the marginal tree for each site.
    /// If multiple sites have the same marginal tree, it is provided multiple times.
    pub(crate) fn sites(
        &self,
        start: VariantIndex,
        end: VariantIndex,
    ) -> impl Iterator<Item=(VariantIndex, MarginalTree)> {
        todo!()
    }
}

/// A node in the marginal tree of the tree sequence iterator.
#[derive(Copy, Clone, Debug, Hash, Eq, PartialEq)]
pub(crate) struct TreeNode(usize);

/// An access structure into the tree sequence at a specific site.
/// This structure is generated and by the [`AncestorIterator`] during iteration.
/// It provides method for tree node compression during the viterbi algorithm.
///
/// # Free Nodes
/// The structure might include some nodes that aren't part of the tree,
/// if no edges have been calculated for those nodes yet.
/// As a result, those nodes cannot be recompressed (they have no parents).
/// The marginal tree offers methods to iterate through all uncompressed nodes.
///
/// [`AncestorIterator`]: AncestorIterator
pub(crate) struct MarginalTree {
    actual_parents: Vec<Option<TreeNode>>,
    compressed_tree_parents: Vec<Option<TreeNode>>,
}

impl MarginalTree {
    /// An iterator through the currently uncompressed nodes
    pub(crate) fn nodes(&self) -> impl Iterator<Item=TreeNode> {
        todo!()
    }

    /// Recompress a node into its parent node. The method assumes that the node has the same
    /// likelihood as its parent, otherwise recompression is a logical error in Viterbi.
    pub(crate) fn recompress(&mut self, node: TreeNode) {
        todo!()
    }

    /// Get the uncompressed parent of the given tree node.
    /// Returns `None` if the given node is the tree root, or a [free node].
    ///
    /// [free node]: MarginalTree
    pub(crate) fn parent(&self, node: TreeNode) -> Option<TreeNode> {
        todo!()
    }
}
