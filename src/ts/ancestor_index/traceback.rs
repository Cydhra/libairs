use crate::ancestors::Ancestor;
use crate::ts::ancestor_index::tree::MarginalTree;
use crate::ts::ancestor_index::{ViterbiEvent, ViterbiEventKind};
use crate::ts::partial_sequence::PartialTreeSequence;
use crate::variants::VariantIndex;
use std::num::NonZeroUsize;

/// Iterator to iterate back through the viterbi table generated by the Viterbi algorithm.
/// It holds the marginal tree that was modified by the Viterbi algorithm and iterates back through
/// some aspects of the tree sequence, but only on the most likely path.
///
/// TODO this needs way more description
#[derive(Clone)]
pub(in crate::ts) struct TracebackSequenceIterator<'a, 'o>
where
    'o: 'a,
{
    marginal_tree: &'a MarginalTree<'o>,
    partial_tree_sequence: &'a PartialTreeSequence,

    /// start of the interval that is being iterated over. Since we iterate backwards, this marks the end of iteration.
    start: VariantIndex,

    /// end of the interval that is being iterated over. Since we iterate backwards, this marks the start of iteration.
    end: VariantIndex,

    /// site that has last been yielded by the iterator
    current: VariantIndex,

    /// the ancestor that is currently being iterated over
    current_ancestor: Ancestor,

    /// index into the linked list of viterbi events for the next event yielded by the iterator
    next_event: usize,

    /// how long the current_ancessor is compressed
    compressed_until: Option<VariantIndex>,

    /// index of the parent edge that is being used to decompress the current ancestor
    current_parent_edge_index: Option<usize>,

    /// iterator over the viterbi events of the parent ancestor in case of compression
    inner: Option<Box<TracebackSequenceIterator<'a, 'o>>>,
}

impl<'o> TracebackSequenceIterator<'o, 'o> {
    /// Create a new traceback iterator that iterates over the most likely path in the given interval.
    /// The iterator will start at the end of the interval and iterate backwards to the start of the interval.
    /// The iterator will automatically switch to the most recent uncompressed parent if the
    /// `current_ancestor` is compressed.
    pub(super) fn new(
        marginal_tree: &'o MarginalTree,
        partial_tree_sequence: &'o PartialTreeSequence,
        current_ancestor: Ancestor,
        start: VariantIndex,
        end: VariantIndex,
        next_event: usize,
    ) -> Self {
        TracebackSequenceIterator {
            marginal_tree,
            partial_tree_sequence,
            start,
            end,
            current: end,
            current_ancestor,
            next_event,
            compressed_until: None,
            current_parent_edge_index: None,
            inner: None,
        }
    }

    /// Create a new `TracebackSequenceIterator` that iterates over the interval from `start` to `end` in the
    /// `marginal_tree` and `partial_tree_sequence` of the parent iterator.
    fn new_inner(
        parent: &Self,
        ancestor: Ancestor,
        start: VariantIndex,
        end: VariantIndex,
    ) -> Self {
        // find the exclusive end of the interval
        let (mut end_cursor, pred) = parent.search_site(ancestor, end);

        let mut new_iter = Self {
            marginal_tree: parent.marginal_tree,
            partial_tree_sequence: parent.partial_tree_sequence,
            start,
            end,
            current: end,
            current_ancestor: ancestor,
            next_event: end_cursor,
            compressed_until: None,
            current_parent_edge_index: None,
            inner: None,
        };

        // check if the current ancestor is compressed at the end site. If it is, the event directly
        // after the end site must be the compression event
        if pred.is_some() {
            if let ViterbiEventKind::Compressed =
                parent.marginal_tree.linked_viterbi_events[pred.unwrap().get()].kind
            {
                assert!(matches!(
                    parent.marginal_tree.linked_viterbi_events[end_cursor].kind,
                    ViterbiEventKind::Uncompressed
                ));
                let compressed_until = parent.marginal_tree.linked_viterbi_events[end_cursor].site;

                // move end_cursor past the uncompression event
                end_cursor = parent.marginal_tree.linked_viterbi_events[end_cursor]
                    .prev
                    .map(|p| p.get())
                    .unwrap_or(0);

                // if the compressed interval reaches into our interval
                if compressed_until <= end {
                    let parent_index = new_iter.search_parent_index(end);
                    let parent_edge =
                        &new_iter.partial_tree_sequence.edges[ancestor.0 as usize][parent_index];

                    // recursively create a new inner iterator that iterates over the compressed interval
                    let inner = Self::new_inner(
                        &new_iter,
                        parent_edge.parent(),
                        start.max(compressed_until),
                        end,
                    );

                    debug_assert!(inner.current <= new_iter.current);
                    new_iter.inner = Some(Box::new(inner));
                    new_iter.compressed_until = Some(compressed_until);
                    new_iter.current_parent_edge_index = Some(parent_index);
                }
            }
        }

        new_iter
    }

    pub fn switch_ancestor(&mut self, ancestor: Ancestor) {
        self.current_ancestor = ancestor;

        // find the position BEFORE the current site in the new ancestor,
        // since we do not want to copy the current site
        self.next_event = self.search_beyond_site(ancestor, self.current);

        self.destroy_inner();
    }

    /// Return the next event from the current ancestor and decrease the cursor
    fn next_event<'mutable>(&'mutable mut self) -> Option<&'o ViterbiEvent> {
        if self.next_event > 0 {
            let e = &self.marginal_tree.linked_viterbi_events[self.next_event];
            self.next_event = e.prev.map(|p| p.get()).unwrap_or(0);

            Some(e)
        } else {
            None
        }
    }

    /// Search the event index of the given site in the given ancestor.
    /// The index will be of the first site that has a site index less than or equal to the given site.
    /// Used for inner iterators to find the start of the interval, as well as for switching ancestors.
    ///
    /// Returns the index of the first event that has a site index less than or equal to the given site,
    /// and its immediate predecessor, if it exists.
    fn search_site(&self, ancestor: Ancestor, site: VariantIndex) -> (usize, Option<NonZeroUsize>) {
        let mut pred = None;
        let mut current = self.marginal_tree.last_event_index[ancestor.0 as usize];
        while self.marginal_tree.linked_viterbi_events[current].site > site {
            pred = NonZeroUsize::new(current);
            current = self.marginal_tree.linked_viterbi_events[current]
                .prev
                .map(|p| p.get())
                .unwrap_or(0);
        }
        (current, pred)
    }

    /// Search the event index of the given site in the given ancestor.
    /// The index will be of the first site that has a site index less than the given site, meaning
    /// if sites have the given site index, they will be right of the returned index.
    fn search_beyond_site(&self, ancestor: Ancestor, site: VariantIndex) -> usize {
        let mut current = self.marginal_tree.last_event_index[ancestor.0 as usize];
        while current > 0 && self.marginal_tree.linked_viterbi_events[current].site >= site {
            current = self.marginal_tree.linked_viterbi_events[current]
                .prev
                .map(|p| p.get())
                .unwrap_or(0);
        }
        current
    }

    /// Searches the parent of the current ancestor for the given site
    fn search_parent_index(&self, site: VariantIndex) -> usize {
        self.partial_tree_sequence.edges[self.current_ancestor.0 as usize]
            .binary_search_by(|edge| edge.start().cmp(&site))
            .unwrap_or_else(|next_element_index| {
                debug_assert!(next_element_index > 0);
                next_element_index - 1
            })
    }

    /// Used by next() to end the iterator
    fn end_iteration(&mut self) {
        self.current = self.end;
        self.inner = None;
        self.current_parent_edge_index = None;
        self.compressed_until = None;
    }

    /// end iteration of inner iterator and continue to iterate the outer iterator
    fn destroy_inner(&mut self) {
        self.inner = None;
        self.compressed_until = None;
        self.current_parent_edge_index = None;
    }
}

impl<'o> Iterator for TracebackSequenceIterator<'o, 'o> {
    type Item = &'o ViterbiEvent;

    fn next(&mut self) -> Option<Self::Item> {
        // if we have an inner iterator, we need to check if it has any elements left and return those
        if let Some(ref mut inner) = self.inner {
            let element = inner.next();
            debug_assert!(element.is_none() || element.unwrap().site <= self.current);

            if element.is_none() {
                // if the inner iterator is empty, we arrived at the site it ends
                self.current = self.inner.as_ref().unwrap().start;

                // if we are at the start of the interval, we end the iteration
                if self.current == self.start {
                    self.end_iteration();
                    return None;
                }

                debug_assert!(self.compressed_until.is_some());
                // if we haven't reached an uncompressed interval, we need to switch to the next parent
                if self.compressed_until.unwrap() < self.current {
                    assert!(self.current_parent_edge_index.unwrap() > 0, "last inner iterator already reached the first parent, but self.current has not reached start yet");
                    let next_parent_index = self.current_parent_edge_index.unwrap() - 1;

                    let next_parent_edge = &self.partial_tree_sequence.edges
                        [self.current_ancestor.0 as usize][next_parent_index];

                    // we should not be in a position where the current site is already beyond the
                    // end of the next edge, since that would mean the last inner iterator yielded
                    // wrong elements
                    debug_assert!(next_parent_edge.end() < self.current);

                    // if we reached the end of interval, we end the iterator.
                    // this happens if self.current has not reached self.start yet, but the last
                    // inner iterator has no more elements between self.current and self.start,
                    // yet the parent is defined until self.start
                    if next_parent_edge.end() < self.start {
                        self.end_iteration();
                        return None;
                    }

                    // if the next parent reaches into the interval of interest, update the inner
                    // iterator, otherwise destroy it
                    if next_parent_edge.end() > self.compressed_until.unwrap() {
                        let new_inner = Self::new_inner(
                            self,
                            next_parent_edge.parent(),
                            next_parent_edge.start().max(self.compressed_until.unwrap()),
                            self.current,
                        );

                        debug_assert!(new_inner.current <= self.current);
                        self.inner = Some(Box::new(new_inner));
                        self.current_parent_edge_index = Some(next_parent_index);

                        return self.next();
                    } else {
                        self.destroy_inner();
                    }
                } else {
                    self.destroy_inner();
                }
            } else {
                self.current = element.unwrap().site;
                return element;
            }
        }

        // otherwise, we look at the current sequence and return the next element
        let element = self.next_event();
        debug_assert!(
            element.is_none() || element.unwrap().site <= self.current,
            "current: {}, element: {:?}",
            self.current,
            element
        );
        if let Some(event) = element {
            if event.site < self.start {
                self.end_iteration();
                return None;
            } else {
                self.current = event.site;
            }

            match event.kind {
                ViterbiEventKind::Compressed => {
                    let compressed_until =
                        self.marginal_tree.linked_viterbi_events[event.prev.unwrap().get()].site;

                    // search for the parent of the compressed node
                    let current_parent_edge_index = self.search_parent_index(event.site);
                    let parent = self.partial_tree_sequence.edges[self.current_ancestor.0 as usize]
                        [current_parent_edge_index]
                        .parent();

                    let inner =
                        Self::new_inner(self, parent, self.start.max(compressed_until), event.site);
                    debug_assert!(inner.current <= self.current);

                    self.inner = Some(Box::new(inner));
                    self.compressed_until = Some(compressed_until);
                    self.current_parent_edge_index = Some(current_parent_edge_index);
                    self.next()
                }
                ViterbiEventKind::Sentinel => {
                    panic!("Sentinel event in traceback iterator")
                }
                _ => element,
            }
        } else {
            self.end_iteration();
            None
        }
    }
}
