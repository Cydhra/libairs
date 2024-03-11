use std::cmp::Ordering;
use std::collections::BTreeSet;
use std::iter::Peekable;

use crate::ancestors::Ancestor;
use crate::ts::partial_sequence::PartialSequenceEdge;
use crate::variants::VariantIndex;

pub(crate) type Site<'a, 'o> = (VariantIndex, &'a mut MarginalTree<'o>);

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

    /// Recombination sites array
    recombination_sites: Vec<Vec<bool>>,

    /// Mutation sites array
    mutation_sites: Vec<Vec<bool>>,

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: Vec<usize>,

    /// pre-allocated stack for the DFS in subtree updates
    stack: Vec<Ancestor>,
}

impl AncestorIndex {
    /// Create a new empty ancestor iterator
    pub(crate) fn new(max_nodes: usize, variant_count: usize) -> Self {
        Self {
            edge_index: BTreeSet::new(),
            num_nodes: 1,
            parents: vec![None; max_nodes],
            children: vec![Vec::new(); max_nodes],
            uncompressed_parents: vec![None; max_nodes],
            is_compressed: vec![true; max_nodes],
            active_nodes: Vec::new(),
            likelihoods: vec![-1.0f64; max_nodes],
            recombination_sites: vec![vec![false; variant_count]; max_nodes],
            mutation_sites: vec![vec![false; variant_count]; max_nodes],
            last_compressed: vec![0; max_nodes],
            stack: Vec::new(),
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
        edges: &[PartialSequenceEdge],
        mutations: &[VariantIndex],
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
    /// It provides the marginal tree for each site.
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
            &mut self.recombination_sites[0..limit_nodes],
            &mut self.mutation_sites[0..limit_nodes],
            &mut self.last_compressed[0..limit_nodes],
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
    Start { parent: Ancestor },
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

    /// Recombination sites array
    recombination_sites: &'o mut [Vec<bool>],

    /// Mutation sites array
    mutation_sites: &'o mut [Vec<bool>],

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    last_compressed: &'o mut [usize],

    /// The number of nodes in the tree. Nodes with a higher index are ignored in the tree.
    limit_nodes: usize,

    /// Start index of the current iteration
    start: VariantIndex,

    /// Length of the sequence processed by the current tree
    current_sequence_length: usize,
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
        recombination_sites: &'o mut [Vec<bool>],
        mutation_sites: &'o mut [Vec<bool>],
        last_compressed: &'o mut [usize],
    ) -> Self {
        debug_assert!(num_nodes > 0, "Tree must have at least one node");

        // re-initialize vectors into the default state where needed
        active_nodes.clear();
        recombination_sites
            .iter_mut()
            .for_each(|i| i[0..ancestor_length].fill(false));
        mutation_sites
            .iter_mut()
            .for_each(|i| i[0..ancestor_length].fill(false));
        last_compressed.fill(0);

        // the other states are updated whenever nodes are added to the marginal tree

        let mut marginal_tree = Self {
            parents,
            children,
            uncompressed_parents,
            is_compressed,
            active_nodes,
            likelihoods,
            recombination_sites,
            mutation_sites,
            last_compressed,
            limit_nodes,
            start,
            current_sequence_length: ancestor_length,
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
        let mut parent = self.parents[node.0];

        while let Some(p) = parent {
            if !self.is_compressed(p) {
                break;
            }
            parent = self.parents[p.0];
        }

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
        // Todo do this in a single pass loop
        self.active_nodes.iter().for_each(|&node| {
            debug_assert!(self.is_compressed[node.0] == false);

            if let Some(parent) = self.uncompressed_parents[node.0] {
                if self.likelihoods[node.0] == self.likelihoods[parent.0] {
                    recompress.push(node);
                }
            }
        });

        recompress.into_iter().for_each(|node| {
            self.recompress_node(node, self.uncompressed_parents[node.0].unwrap(), site_index);
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
        debug_assert!(node.0 != 0, "Root node cannot be compressed");

        self.active_nodes.retain(|&n| n != node);
        // store the first site where the node is compressed, so the next one
        self.last_compressed[node.0] = site_index + 1;
        self.is_compressed[node.0] = true;

        self.update_subtree(node, Some(node), Some(parent));
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

            // copy the recombination and mutation sites from the parent
            let last_compressed_begin = self.last_compressed[node.0];

            Self::copy_parent_sites(
                &mut self.recombination_sites,
                uncompressed_parent,
                node,
                last_compressed_begin,
                site_index,
            );
            Self::copy_parent_sites(
                &mut self.mutation_sites,
                uncompressed_parent,
                node,
                last_compressed_begin,
                site_index,
            );
        }
    }

    /// Helper function to copy slices from two-dimensional arrays from a parent node to a child node. The code is a
    /// bit awkward because we need to convince the borrow checker that the slices don't intersect.
    fn copy_parent_sites(
        sites: &mut [Vec<bool>],
        parent: Ancestor,
        node: Ancestor,
        last_compressed_begin: usize,
        site_index: usize,
    ) {
        debug_assert!(parent != node);
        debug_assert!(parent < node, "Parent must be older than child");

        let (parent_sites, child_sites) = sites.split_at_mut(node.0);
        child_sites[0][last_compressed_begin..site_index]
            .copy_from_slice(&parent_sites[parent.0][last_compressed_begin..site_index]);
    }

    /// Insert a node into the tree that has no parent, i.e. it is a free node or is the root node.
    /// This method is used to initialize the tree sequence, and should not be called during the
    /// Viterbi algorithm.
    fn add_initial_node(&mut self, node: Ancestor) {
        debug_assert!(self.parents[node.0].is_none());

        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        self.likelihoods[node.0] = 1.0;
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
                self.uncompressed_parents[node.0] = self.uncompressed_parents[new_parent.0];
            } else {
                self.uncompressed_parents[node.0] = Some(new_parent);
            }

            self.update_subtree(node, old_parent, Some(node))
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
            self.update_subtree(node, self.uncompressed_parents[node.0], Some(node));
        }
    }

    /// Update a subtree to reflect a change of compression of its root.
    ///
    /// `old_parent` is an option for convenience, but cannot be `None`
    /// `new_parent` is an option for convenience, but cannot be `None`
    fn update_subtree(
        &mut self,
        subtree_root: Ancestor,
        old_parent: Option<Ancestor>,
        new_parent: Option<Ancestor>,
    ) {
        debug_assert!(old_parent.is_some());
        debug_assert!(new_parent.is_some());

        let mut queue = self.children[subtree_root.0].clone();
        while let Some(node) = queue.pop() {
            if self.uncompressed_parents[node.0] == old_parent {
                self.uncompressed_parents[node.0] = new_parent;

                queue.extend(self.children[node.0].iter().copied());
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

        self.children[self.parents[node.0].unwrap().0].retain(|&n| n != node);
        self.parents[node.0] = None;
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
    use std::hint::black_box;

    #[test]
    fn test_empty_range() {
        // test whether iterating over an empty range doesn't crash and calls the closure zero times
        let mut ix = AncestorIndex::new(2, 1);
        ix.sites(
            VariantIndex::from_usize(10),
            VariantIndex::from_usize(10),
            2,
        )
        .for_each(|_| panic!("closure called for empty range"));
    }

    #[test]
    fn test_ancestor_iteration() {
        let mut ix = AncestorIndex::new(2, 10);
        let mut counter = 0;

        ix.sites(VariantIndex::from_usize(0), VariantIndex::from_usize(10), 2)
            .for_each(|(site, tree)| {
                assert_eq!(site, VariantIndex::from_usize(counter));
                assert_eq!(tree.num_nodes(), 2);
                assert_eq!(tree.nodes().count(), 2);
                assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                assert_eq!(tree.find_uncompressed_parent(Ancestor(1)), None);
                counter += 1;
            });
    }

    #[test]
    fn test_simple_tree_compression() {
        let mut ix = AncestorIndex::new(2, 10);
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
        let mut ix = AncestorIndex::new(2, 10);
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
        let mut ix = AncestorIndex::new(3, 9);
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
        let mut ix = AncestorIndex::new(2, 10);
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
        let mut ix = AncestorIndex::new(3, 10);
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
        let mut ix = AncestorIndex::new(2, 10);
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
        let mut ix = AncestorIndex::new(3, 10);
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
        let mut ix = AncestorIndex::new(2, 10);
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

        let mut ix = AncestorIndex::new(3, 10);
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
                2 => *tree.mutation_site(Ancestor(0), 2) = true,
                3 => *tree.recombination_site(Ancestor(0), 3) = true,
                4 => {
                    tree.likelihoods[2] = 0.1; // set likelihood of third node to that of second to trigger recompression
                    *tree.mutation_site(Ancestor(0), 4) = true // this shouldnt be copied, because the second node is uncompressed at this site
                }
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
                2 => assert!(*state, "mutation site {} should be true", i),
                6 => assert!(*state, "mutation site {} should be true", i),
                _ => assert!(!*state, "mutation site {} should be false", i),
            });

        iterator.marginal_tree.recombination_sites[2]
            .iter()
            .enumerate()
            .for_each(|(i, state)| match i {
                3 => assert!(*state, "recombination site {} should be true", i),
                7 => assert!(*state, "recombination site {} should be true", i),
                _ => assert!(!*state, "recombination site {} should be false", i),
            });
    }

    #[test]
    fn test_copying_from_parent_on_change_parent() {
        // test whether the iterator copies the recombination and mutation sites from the parent when a child is
        // decompressed

        let mut ix = AncestorIndex::new(3, 10);
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
                2 => *tree.mutation_site(Ancestor(0), 2) = true,
                3 => *tree.recombination_site(Ancestor(0), 3) = true,
                5 => *tree.mutation_site(Ancestor(1), 5) = true, // this should not be copied because the second ancestor is uncompressed at this site

                _ => {}
            }

            counter += 1;
        });

        iterator.marginal_tree.mutation_sites[2]
            .iter()
            .enumerate()
            .for_each(|(i, state)| match i {
                2 => assert!(*state, "mutation site {} should be true", i),
                _ => assert!(!*state, "mutation site {} should be false", i),
            });

        iterator.marginal_tree.recombination_sites[2]
            .iter()
            .enumerate()
            .for_each(|(i, state)| match i {
                3 => assert!(*state, "recombination site {} should be true", i),
                _ => assert!(!*state, "recombination site {} should be false", i),
            });
    }

    #[test]
    fn test_queue_order() {
        // test whether the queue order is sensible in upholding invariants.
        // specifically it shouldn't remove nodes before their children have changed parents or have been removed,
        // and shouldnt insert nodes before their parents have been inserted

        let mut ix = AncestorIndex::new(4, 10);

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
        let mut ix = AncestorIndex::new(3, 10);
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
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(1)),
                            Some(Ancestor(0))
                        );
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(2)),
                            Some(Ancestor(1))
                        );
                    }
                    1 => {
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                        // ancestor 1 is compressed
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(2)),
                            Some(Ancestor(0))
                        );
                    }
                    2 => {
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(1)),
                            Some(Ancestor(0))
                        );
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(2)),
                            Some(Ancestor(1))
                        );
                    }
                    4 => {
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(1)),
                            Some(Ancestor(0))
                        );
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(2)),
                            Some(Ancestor(1))
                        );

                        // prevent recompression, to test if Ancestor(2) still recombines to Ancestor(0)
                        *tree.likelihood(Ancestor(1)) = 0.1;
                    }
                    5 => {
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(1)),
                            Some(Ancestor(0))
                        );
                        assert_eq!(
                            tree.find_uncompressed_parent(Ancestor(2)),
                            Some(Ancestor(0))
                        );

                        // recompress everything
                        *tree.likelihood(Ancestor(1)) = 0.0;
                    }
                    _ => {
                        assert_eq!(tree.find_uncompressed_parent(Ancestor(0)), None);
                    }
                }
                counter += 1usize;
            });
    }
}
