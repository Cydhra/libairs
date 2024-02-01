use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::iter::Peekable;

use crate::ancestors::{Ancestor, VariantIndex};
use crate::ts::partial_sequence::PartialSequenceEdge;

pub(crate) type Site<'a> = (VariantIndex, &'a mut MarginalTree);

/// A helper structure for the Viterbi algorithm that helps iterating through the sites, updating
/// the marginal tree at the currently visited site and helping with tree compression.
/// The iterator only produces indices into the ancestor sequence structure, it does not hold the
/// actual ancestors.
///
/// Whenever new tree edges are calculated, they can be inserted into the iterator so the tree
/// compression improves for future iterations.
pub(crate) struct AncestorIndex {
    edge_index: BTreeSet<SequenceEvent>,
    num_nodes: usize, // number of nodes in the tree sequence. All nodes with an index >= this
                      // value are not part of the tree sequence yet, and are thus free nodes.
                      // The root node is always part of the tree, so this value is at least 1.
}

impl AncestorIndex {
    /// Create a new empty ancestor iterator
    pub(crate) fn new() -> Self {
        Self {
            edge_index: BTreeSet::new(),
            num_nodes: 1,
        }
    }

    /// Insert a tree sequence node into the partial tree sequence, which will be exploited by the
    /// iterator to compress trees.
    /// The nodes must be inserted in order of increasing index.
    ///
    /// # Parameters
    /// - `tree_node`: The tree node handle that will be inserted into the tree sequence
    /// - `edges`: The edges that are associated with the tree node
    /// - `mutations`: The mutations that are associated with the tree node
    ///
    /// # Panics
    /// Panics if the given tree node index is not the next index in the sequence.
    pub(crate) fn insert_sequence_node(
        &mut self,
        tree_node: Ancestor,
        edges: Vec<PartialSequenceEdge>,
        mutations: Vec<VariantIndex>,
    ) {
        self.num_nodes += 1;
        assert_eq!(
            tree_node.0,
            self.num_nodes - 1,
            "Tree nodes must be inserted in order"
        );

        if !edges.is_empty() {
            self.edge_index.insert(SequenceEvent {
                site: edges.last().unwrap().end(),
                node: tree_node,
                kind: SequenceEventKind::End,
            });

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

        for mutation in mutations {
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
    /// It provides the marginal tree for each site.
    /// If multiple sites have the same marginal tree, it is provided multiple times.
    ///
    /// # Parameters
    /// - `start`: The index of the first site to iterate through
    /// - `end`: The index of the last site to iterate through (exclusive)
    /// - `limit_nodes`: The number of nodes to consider for the marginal tree. All nodes with a
    /// index equal or greater than this value will be ignored.
    pub(crate) fn sites(
        &self,
        start: VariantIndex,
        end: VariantIndex,
        limit_nodes: usize,
    ) -> PartialTreeSequenceIterator<impl Iterator<Item = &'_ SequenceEvent>> {
        let mut marginal_tree = MarginalTree::new(
            start,
            self.num_nodes,
            limit_nodes,
            end.get_variant_distance(start),
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
/// Not an actual implementation of IterMut, because it borrows the same marginal tree for each site
/// mutably, so we need to make sure the reference cannot escape the iterator.
pub(crate) struct PartialTreeSequenceIterator<'a, I: Iterator<Item = &'a SequenceEvent>> {
    marginal_tree: MarginalTree,
    start: VariantIndex,
    site: VariantIndex,
    end: VariantIndex,
    queue: Peekable<I>,
}

impl<'a, I: Iterator<Item = &'a SequenceEvent>> PartialTreeSequenceIterator<'a, I> {
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
            self.marginal_tree
                .recompress_tree(self.site.get_variant_distance(self.start));
            self.site = self.site.next();
        }
    }

    /// Get the flag whether the given node requires recombination at the given index.
    /// The index is a zero-based integer beginning at the first site of the current candidate, meaning it is not a
    /// variant site index.
    pub(crate) fn recombination_site(&self, node: Ancestor, index: usize) -> bool {
        self.marginal_tree.recombination_sites[node.0][index]
    }

    /// Get the flag whether the given node requires mutation at the given site
    /// The index is a zero-based integer beginning at the first site of the current candidate, meaning it is not a
    /// variant site index.
    pub(crate) fn mutation_site(&self, node: Ancestor, index: usize) -> bool {
        self.marginal_tree.mutation_sites[node.0][index]
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum SequenceEventKind {
    End,
    Start { parent: Ancestor },
    ChangeParent { new_parent: Ancestor },
    Mutation,
}

impl PartialOrd<Self> for SequenceEventKind {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SequenceEventKind {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::End, Self::End) => Ordering::Equal,
            // we can afford to make this the smallest element, because it cannot be at the same site
            // as any other event regarding the same node. And it has advantages to make this value
            // a sentinel rather than one with parameters
            (Self::End, _) => Ordering::Greater,
            (_, Self::End) => Ordering::Less,
            (Self::Start { parent: p1 }, Self::Start { parent: p2 }) => p1.cmp(p2),
            (Self::Start { .. }, _) => Ordering::Less,
            (_, Self::Start { .. }) => Ordering::Greater,
            (Self::ChangeParent { new_parent: p1 }, Self::ChangeParent { new_parent: p2 }) => {
                p1.cmp(p2)
            }
            (Self::ChangeParent { .. }, _) => Ordering::Less,
            (_, Self::ChangeParent { .. }) => Ordering::Greater,
            (Self::Mutation, Self::Mutation) => Ordering::Equal,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct SequenceEvent {
    site: VariantIndex,
    node: Ancestor,
    kind: SequenceEventKind,
}

impl SequenceEvent {
    /// Returns a sentinel event that is smaller or equal to the smallest event at the given site.
    /// This is useful for defining ranges in the BTreeSet.
    fn sentinel(pos: VariantIndex) -> Self {
        Self {
            site: pos,
            node: Ancestor(0),
            kind: SequenceEventKind::End,
        }
    }
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
            .then(self.node.cmp(&other.node))
            .then(self.kind.cmp(&other.kind))
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
#[derive(Clone, Debug)]
pub(crate) struct MarginalTree {
    /// The actual parents of the tree nodes. None for free nodes and the root, parents might be
    /// compressed.
    actual_parents: Vec<Option<Ancestor>>,

    /// Uncompressed ancestors of uncompressed tree nodes. None for free nodes and the root, wrong
    /// for compressed nodes.
    uncompressed_tree_parents: Vec<Option<Ancestor>>,

    /// Whether the node is currently compressed.
    is_compressed: Vec<bool>,

    /// The list of active nodes, i.e. nodes that are not compressed and are used in the Viterbi
    /// algorithm.
    active_nodes: Vec<Ancestor>,

    /// The likelihood of each node in the tree. Likelihood for compressed nodes is wrong, but
    /// they are ignored during the Viterbi algorithm.
    likelihoods: Vec<f64>,

    /// Recombination sites array
    recombination_sites: Vec<Vec<bool>>,

    /// Mutation sites array
    mutation_sites: Vec<Vec<bool>>,

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: Vec<usize>,

    /// The number of nodes in the tree. Nodes with a higher index are ignored in the tree.
    limit_nodes: usize,

    /// Start index of the current iteration
    start: VariantIndex,
}

impl MarginalTree {
    fn new(
        start: VariantIndex,
        num_nodes: usize,
        limit_nodes: usize,
        sequence_length: usize,
    ) -> Self {
        debug_assert!(num_nodes > 0, "Tree must have at least one node");

        let mut marginal_tree = Self {
            actual_parents: vec![None; limit_nodes],
            uncompressed_tree_parents: vec![None; limit_nodes],
            is_compressed: vec![true; limit_nodes],
            active_nodes: Vec::new(),
            likelihoods: vec![-1.0f64; limit_nodes],
            recombination_sites: vec![vec![false; sequence_length]; limit_nodes],
            mutation_sites: vec![vec![false; sequence_length]; limit_nodes],
            last_compressed: vec![0; limit_nodes],
            limit_nodes,
            start,
        };

        marginal_tree.add_initial_node(Ancestor(0));
        (num_nodes..limit_nodes)
            .map(Ancestor)
            .for_each(|n| marginal_tree.add_initial_node(n));

        marginal_tree
    }

    /// An iterator through the currently uncompressed nodes
    pub(crate) fn nodes(&self) -> impl Iterator<Item = Ancestor> {
        // TODO can we maybe get away without cloning here?
        self.active_nodes.clone().into_iter()
    }

    /// Get the parent of the given tree node ignoring compressed nodes.
    /// Returns `None` if the given node is the tree root, or a [free node].
    ///
    /// [free node]: MarginalTree
    pub(crate) fn parent(&self, node: Ancestor) -> Option<Ancestor> {
        // compressed parents arent tracked for compressed nodes
        debug_assert!(
            self.is_compressed[node.0] == false,
            "Cannot get parent of compressed node"
        );

        self.uncompressed_tree_parents[node.0]
    }

    /// Get the number of active nodes in the tree.
    /// This includes nodes that cannot be copied from at the current location, but might
    /// influence nodes that can be copied from.
    /// This also includes free nodes that aren't matched against the tree yet.
    /// This value is interesting mostly for the Markov Chain during the Viterbi algorithm.
    pub(crate) fn num_nodes(&self) -> usize {
        self.likelihoods.len()
    }

    /// Get the likelihood of the given node.
    pub(crate) fn likelihood(&mut self, node: Ancestor) -> &mut f64 {
        debug_assert!(!self.is_compressed[node.0]);
        &mut self.likelihoods[node.0]
    }

    /// Get access to the flag whether the given node requires recombination at the given index.
    /// The index is a zero-based integer beginning at the first site of the current candidate, meaning it is not a
    /// variant site index.
    pub(crate) fn recombination_site(&mut self, node: Ancestor, index: usize) -> &mut bool {
        debug_assert!(!self.is_compressed[node.0]);
        &mut self.recombination_sites[node.0][index]
    }

    /// Get access to the flag whether the given node requires mutation at the given site
    /// The index is a zero-based integer beginning at the first site of the current candidate, meaning it is not a
    /// variant site index.
    pub(crate) fn mutation_site(&mut self, node: Ancestor, index: usize) -> &mut bool {
        debug_assert!(!self.is_compressed[node.0]);
        &mut self.mutation_sites[node.0][index]
    }

    fn find_uncompressed_parent(&self, node: Ancestor) -> Option<Ancestor> {
        // If the node isn't compressed, the parent is recorded in the uncompressed tree parents array
        debug_assert!(
            self.is_compressed[node.0],
            "Uncompressed parent requested for uncompressed node {}",
            node.0
        );

        let mut parent = self.actual_parents[node.0];

        while let Some(p) = parent {
            if !self.is_compressed(p) {
                break;
            }
            parent = self.actual_parents[p.0];
        }
        debug_assert!(
            parent.is_some(),
            "No uncompressed ancestor found for node {}",
            node.0
        );
        parent
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
        let mut recompress = Vec::new();
        self.active_nodes.iter().for_each(|&node| {
            debug_assert!(self.is_compressed[node.0] == false);

            if let Some(parent) = self.parent(node) {
                if self.likelihoods[node.0] == self.likelihoods[parent.0] {
                    recompress.push(node);
                }
            }
        });

        recompress.into_iter().for_each(|node| {
            self.recompress_node(
                node,
                self.parent(node)
                    .expect("node has no parent but was marked for recompression"),
                site_index,
            );
        });
    }

    /// Recompress a node into its parent node. The method assumes that the node has the same
    /// likelihood as its parent, otherwise recompression is a logical error in Viterbi.
    ///
    /// # Parameters
    /// - `node`: The node to recompress
    /// - `parent`: The parent of the node
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    fn recompress_node(&mut self, node: Ancestor, parent: Ancestor, site_index: usize) {
        debug_assert!(self.likelihoods[node.0] == self.likelihoods[parent.0]);
        debug_assert!(self.is_compressed[node.0] == false);
        debug_assert!(self.is_compressed[parent.0] == false);

        self.uncompressed_tree_parents.iter_mut().for_each(|entry| {
            if *entry == Some(node) {
                *entry = Some(parent);
            }
        });
        self.active_nodes.retain(|&n| n != node);
        self.last_compressed[node.0] = site_index;
        self.is_compressed[node.0] = true;
    }

    /// Search the parent of the given node and update the uncompressed_tree_parents array, the is_compressed array
    /// and add the node to active nodes.
    /// Returns the parent.
    fn set_uncompressed(&mut self, node: Ancestor) -> Ancestor {
        let parent = self.find_uncompressed_parent(node);
        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        self.uncompressed_tree_parents[node.0] = parent;
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
            let parent = self.set_uncompressed(node);

            self.likelihoods[node.0] = self.likelihoods[parent.0];

            /// copy the recombination and mutation sites from the parent
            let last_compressed_begin = self.last_compressed[node.0];

            Self::copy_parent_sites(
                &mut self.recombination_sites,
                parent,
                node,
                last_compressed_begin,
                site_index,
            );
            Self::copy_parent_sites(
                &mut self.mutation_sites,
                parent,
                node,
                last_compressed_begin,
                site_index,
            );
        }
    }

    /// Helper function to copy slices from two-dimensional arrays from a parent node to a child node. The code is a
    /// bit awkward because we need to convince the borrow checker that the slices don't intersect.
    fn copy_parent_sites(
        sites: &mut Vec<Vec<bool>>,
        parent: Ancestor,
        node: Ancestor,
        last_compressed_begin: usize,
        site_index: usize,
    ) {
        assert_ne!(parent, node);

        if node.0 < parent.0 {
            let (a, b) = sites.split_at_mut(parent.0);
            a[node.0][last_compressed_begin..site_index]
                .copy_from_slice(&b[0][last_compressed_begin..site_index]);
        } else {
            let (a, b) = sites.split_at_mut(node.0);
            b[0][last_compressed_begin..site_index]
                .copy_from_slice(&a[parent.0][last_compressed_begin..site_index]);
        }
    }

    /// Insert a node into the tree that has no parent, i.e. it is a free node or is the root node.
    /// This method is used to initialize the tree sequence, and should not be called during the
    /// Viterbi algorithm.
    fn add_initial_node(&mut self, node: Ancestor) {
        debug_assert!(self.actual_parents[node.0].is_none());
        debug_assert!(self.uncompressed_tree_parents[node.0].is_none());

        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        self.likelihoods[node.0] = 1.0;
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
        self.actual_parents[child.0] = Some(parent);

        if !keep_compressed {
            self.set_uncompressed(child);
            self.likelihoods[child.0] = 0.0;
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
        self.actual_parents[node.0] = Some(new_parent);

        if !keep_compressed {
            self.ensure_decompressed(node, site_index);
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
        self.ensure_decompressed(node, site_index)
    }

    /// Remove a node from the tree. This is done whenever an ancestor ends.
    fn remove_node(&mut self, node: Ancestor) {
        self.actual_parents[node.0] = None;
        self.uncompressed_tree_parents[node.0] = None;
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
    use super::*;

    #[test]
    fn test_ancestor_iteration() {
        let ix = AncestorIndex::new();
        let mut counter = 0;

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
                assert_eq!(site, VariantIndex::from_usize(counter));
                assert_eq!(tree.num_nodes(), 2);
                assert_eq!(tree.nodes().count(), 2);
                assert_eq!(tree.parent(Ancestor(0)), None);
                assert_eq!(tree.parent(Ancestor(1)), None);
                counter += 1;
            });
    }

    #[test]
    fn test_simple_tree_compression() {
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert edge from first to root node
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(5)],
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert edge from first to root node
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert edges for two nodes
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(9)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            vec![
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
            vec![VariantIndex::from_usize(10)],
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            // diverge on site 1, then again on site 5, recompress in between
            vec![VariantIndex::from_usize(1), VariantIndex::from_usize(5)],
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            vec![
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
            vec![VariantIndex::from_usize(5)],
        );

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3)
            .for_each(|(site, tree)| {
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert two nodes, one will be ignored
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
        );

        // check the tree size is always 2
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
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
        let mut ix = AncestorIndex::new();
        let mut counter = 2;

        // insert two nodes
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(5)],
        );
        ix.insert_sequence_node(
            Ancestor(2),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(1)],
        );

        // check the tree
        ix.sites(VariantIndex::from_usize(2), VariantIndex::from_usize(10), 3)
            .for_each(|(site, tree)| {
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
        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert incomplete ancestor
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(5),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
        );

        // check the tree
        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
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
    fn test_copying_from_parent() {
        // test whether the iterator copies the recombination and mutation sites from the parent when a child is
        // decompressed

        let mut ix = AncestorIndex::new();
        let mut counter = 0;

        // insert a child node that serves as a second node
        ix.insert_sequence_node(
            Ancestor(1),
            vec![PartialSequenceEdge::new(
                VariantIndex::from_usize(0),
                VariantIndex::from_usize(10),
                Ancestor(0),
            )],
            vec![VariantIndex::from_usize(0)],
        );

        // insert another child node that copies from both the root and the child node
        ix.insert_sequence_node(
            Ancestor(2),
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
        );

        // add mutations and recombinations with the root and first child node, and check whether the second child node
        // copies them
        let mut iterator = ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 3);
        iterator.for_each(|(site, tree)| {
            match counter {
                // don't recompress first two nodes
                0 => tree
                    .nodes()
                    .for_each(|n| tree.likelihoods[n.0] = n.0 as f64 * 0.1),

                // insert mutations and recombination for nodes 0 and 1, and later check if 2 inherited them
                2 => *tree.mutation_site(Ancestor(0), 2) = true,
                3 => *tree.recombination_site(Ancestor(0), 3) = true,
                4 => *tree.mutation_site(Ancestor(0), 4) = true, // this shouldnt be copied, because the second node is uncompressed at this site
                6 => *tree.mutation_site(Ancestor(1), 6) = true,
                7 => *tree.recombination_site(Ancestor(1), 7) = true,
                _ => {}
            }

            counter += 1;
        });

        iterator.marginal_tree.mutation_sites[2]
            .iter()
            .enumerate()
            .for_each(|(i, state)| match i {
                2 => assert!(*state),
                6 => assert!(*state),
                _ => assert!(!*state),
            });

        iterator.marginal_tree.recombination_sites[2]
            .iter()
            .enumerate()
            .for_each(|(i, state)| match i {
                3 => assert!(*state),
                7 => assert!(*state),
                _ => assert!(!*state),
            });
    }
}
