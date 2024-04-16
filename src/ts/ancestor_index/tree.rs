use crate::ancestors::Ancestor;
use crate::ts::ancestor_index::sequence::{SequenceEvent, SequenceEventKind};
use crate::ts::ancestor_index::{ViterbiEvent, ViterbiEventKind};
use crate::variants::VariantIndex;
use std::iter::Peekable;

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
/// [`AncestorIterator`]: ViterbiIterator
#[derive(Debug)]
pub(crate) struct MarginalTree<'o> {
    /// The actual parents of the tree nodes. `None` for free nodes and the root, parents might be
    /// compressed.
    pub(super) parents: &'o mut [Option<Ancestor>],

    /// Node children
    pub(super) children: &'o mut [Vec<Ancestor>],

    /// The most recent uncompressed ancestor node of each node in the marginal tree.
    /// `None` for root and free nodes
    pub(super) uncompressed_parents: &'o mut [Option<Ancestor>],

    /// Whether the node is currently compressed.
    pub(super) is_compressed: &'o mut [bool],

    /// The list of active nodes, i.e. nodes that are not compressed and are used in the Viterbi
    /// algorithm.
    pub(super) active_nodes: &'o mut Vec<Ancestor>,

    /// The likelihood of each node in the tree. Likelihood for compressed nodes is wrong, but
    /// they are ignored during the Viterbi algorithm.
    pub(super) likelihoods: &'o mut [f64],

    /// Mutations, Recombinations and indirections
    pub(super) viterbi_events: &'o mut [Vec<ViterbiEvent>],

    /// The last site index (relative to current ancestor, not a variant index) where the node was compressed.
    /// Starting from that index, the mutation and recombination sites of that node are invalid
    pub(super) last_compressed: &'o mut [VariantIndex],

    /// The number of nodes in the tree. Nodes with a higher index are ignored in the tree.
    pub(super) limit_nodes: usize,

    /// Start index of the current iteration
    pub(super) start: VariantIndex,

    /// Whether to use the recompression threshold to avoid recompressing when only a few nodes are
    /// active. If this is false, the algorithm will recompress in every iteration, regardless of
    /// efficacy.
    pub(super) use_recompression_threshold: bool,

    /// Inverse recompression threshold.
    /// A threshold of 100 means that recompression of the marginal tree is attempted when
    /// more than 1% of the nodes are active,
    /// a threshold of 50 means that recompression is attempted when more than 2% of the nodes are active.
    /// A threshold of 1 means that recompression is attempted when more than 100% of the nodes are active,
    /// i.e. never.
    pub(super) inv_recompression_threshold: u16,
}

impl<'o> MarginalTree<'o> {
    pub(super) fn new(
        start: VariantIndex,
        num_nodes: usize,
        limit_nodes: usize,
        parents: &'o mut [Option<Ancestor>],
        children: &'o mut [Vec<Ancestor>],
        uncompressed_parents: &'o mut [Option<Ancestor>],
        is_compressed: &'o mut [bool],
        active_nodes: &'o mut Vec<Ancestor>,
        likelihoods: &'o mut [f64],
        viterbi_events: &'o mut [Vec<ViterbiEvent>],
        last_compressed: &'o mut [VariantIndex],
        use_recompression_threshold: bool,
        inv_recompression_threshold: u16,
    ) -> Self {
        debug_assert!(num_nodes > 0, "Tree must have at least one node");

        // re-initialize vectors into the default state where needed
        active_nodes.clear();
        viterbi_events.iter_mut().for_each(|i| i.clear());
        last_compressed.fill(VariantIndex(0));

        // the other states are updated whenever nodes are added to the marginal tree

        let mut marginal_tree = Self {
            parents,
            children,
            uncompressed_parents,
            is_compressed,
            active_nodes,
            likelihoods,
            viterbi_events,
            last_compressed,
            limit_nodes,
            start,
            use_recompression_threshold,
            inv_recompression_threshold,
        };

        marginal_tree.add_initial_node(Ancestor(0), false);
        assert!(
            limit_nodes <= num_nodes,
            "Limit nodes must be less than or equal to num nodes"
        );

        marginal_tree
    }

    /// An iterator through the currently uncompressed nodes
    pub fn nodes(&self) -> impl Iterator<Item = Ancestor> {
        // TODO can we maybe get away without cloning here?
        self.active_nodes.clone().into_iter()
    }

    /// Get the number of valid nodes in the tree.
    /// This includes nodes that cannot be copied from at the current location, but might
    /// influence nodes that can be copied from.
    /// This also includes free nodes that aren't inserted into the tree yet.
    /// This value is interesting mostly for the Markov Chain during the Viterbi algorithm.
    pub fn num_nodes(&self) -> usize {
        self.likelihoods.len()
    }

    /// Get the likelihood of the given node.
    pub fn likelihood(&mut self, node: Ancestor) -> &mut f64 {
        debug_assert!(!self.is_compressed[node.0]);
        &mut self.likelihoods[node.0]
    }

    /// Insert a recombination for the given ancestor at the given site (absolute site, not relative
    /// to the candidate)
    pub fn insert_recombination_event(&mut self, node: Ancestor, site: VariantIndex) {
        self.viterbi_events[node.0].push(ViterbiEvent {
            kind: ViterbiEventKind::Recombination,
            site,
        });
    }

    /// Insert a mutation for the given ancestor at the given site (absolute site, not relative to the
    /// candidate)
    pub fn insert_mutation_event(&mut self, node: Ancestor, site: VariantIndex) {
        self.viterbi_events[node.0].push(ViterbiEvent {
            kind: ViterbiEventKind::Mutation,
            site,
        });
    }

    /// Return whether the node is currently compressed.
    pub(super) fn is_compressed(&self, node: Ancestor) -> bool {
        self.is_compressed[node.0]
    }

    /// Compress all active nodes that have the same likelihood as their parent, so they can be
    /// ignored during the Viterbi algorithm.
    ///
    /// # Parameters
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
    pub fn recompress_tree(&mut self, site_index: VariantIndex) {
        self.active_nodes.retain(|&node| {
            debug_assert!(self.is_compressed[node.0] == false);

            if let Some(parent) = self.uncompressed_parents[node.0] {
                if self.likelihoods[node.0] == self.likelihoods[parent.0] {
                    self.last_compressed[node.0] = site_index + 1;
                    self.is_compressed[node.0] = true;

                    let mut queue = self.children[node.0].clone();
                    while let Some(child) = queue.pop() {
                        if self.uncompressed_parents[child.0] == Some(node) {
                            self.uncompressed_parents[child.0] = Some(parent);

                            if self.is_compressed[child.0] {
                                queue.extend(self.children[child.0].iter().copied());
                            }
                        }
                    }
                    false
                } else {
                    true
                }
            } else {
                true
            }
        });
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
    fn ensure_decompressed(&mut self, node: Ancestor, site_index: VariantIndex) {
        if self.is_compressed(node) {
            debug_assert!(!self.active_nodes.contains(&node));

            // search the first uncompressed ancestor node
            let uncompressed_parent = self.set_uncompressed(node);

            self.likelihoods[node.0] = self.likelihoods[uncompressed_parent.0];

            // copy the recombination and mutation sites from the parent
            let last_compressed_begin = self.last_compressed[node.0];

            // record event for traceback that starting from here we are compressed into the parent
            if site_index > self.start {
                self.viterbi_events[node.0].push(ViterbiEvent {
                    kind: ViterbiEventKind::Compressed(last_compressed_begin),
                    site: site_index - 1,
                });
            }
        }
    }

    /// Insert a node into the tree that has no parent, i.e. it is a free node or is the root node.
    /// This method is used to initialize the tree sequence and to insert unmatched nodes.
    ///
    /// # Parameters
    /// - `node` the free node to insert
    /// - `during_viterbi` must be true, if the viterbi algorithm is already running
    fn add_initial_node(&mut self, node: Ancestor, during_viterbi: bool) {
        debug_assert!(self.parents[node.0].is_none());

        self.active_nodes.push(node);
        self.is_compressed[node.0] = false;
        self.likelihoods[node.0] = if during_viterbi { 0.0 } else { 1.0 };
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
                debug_assert!(
                    self.uncompressed_parents[parent.0].is_some(),
                    "{} has no uncompressed parent",
                    parent.0
                );
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
        site_index: VariantIndex,
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
                debug_assert!(self.uncompressed_parents[new_parent.0].is_some());
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
    fn update_node_mutation(&mut self, node: Ancestor, site_index: VariantIndex) {
        if self.is_compressed(node) {
            self.ensure_decompressed(node, site_index);
            self.update_subtree(node, self.uncompressed_parents[node.0], Some(node));
        }
    }

    /// Update a subtree to reflect a change of compression of its root.
    ///
    /// # Parameters
    /// - `subtree_root`: The root of the subtree to update
    /// - `old_parent`: The previous uncompressed parent for the subtree. The parameter is an option
    /// for implementation convenience, but cannot be `None`
    /// - `new_parent` The new uncompressed parent for the subtree. The parameter is an option for
    /// implementation convenience, but cannot be `None`
    /// - `site_index`: The index of the current site relative to the beginning of the candidate.
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

                if self.is_compressed(node) {
                    queue.extend(self.children[node.0].iter().copied());
                }
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

        // if node has no parent, it is a free node and no updates are needed
        if self.parents[node.0].is_some() {
            self.children[self.parents[node.0].unwrap().0].retain(|&n| n != node);
            self.parents[node.0] = None;
        }

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
    pub(super) fn advance_to_site<'b, I: Iterator<Item = &'b SequenceEvent>>(
        &mut self,
        event_queue: &mut Peekable<I>,
        site: VariantIndex,
        keep_compressed: bool,
        mutations_only: bool,
    ) {
        while event_queue.peek().is_some() && event_queue.peek().unwrap().site < site {
            let event = event_queue.next().unwrap();

            let site_index = event.site;

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
                    kind: SequenceEventKind::StartFree,
                } => {
                    self.add_initial_node(*node, !(keep_compressed || mutations_only));
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
