use crate::ancestors::AncestralSequence;
use crate::ancestors::{Ancestor, AncestorArray};
use crate::ts::ancestor_iterator::AncestorIndex;
use crate::ts::partial_sequence::PartialSequenceEdge;

/// A matcher runs the viterbi algorithm for a set of sequences.
/// It will generate a tree sequence from an array of ancestral sequences and can then match
/// samples against the generated tree sequence.
///
/// It exposes two main methods for matching sequences:
/// - [`Matcher::match_ancestors`]: Match the given ancestors into a partial tree sequence
/// - [`Matcher::match_samples`]: Match a set of samples against the partial tree sequence
pub struct Matcher {
    ancestors: AncestorArray,
    ancestor_iterator: AncestorIndex,
}

impl Matcher {
    /// Create a new matcher for the given ancestral sequences
    pub fn new(ancestors: AncestorArray) -> Self {
        let ancestor_iterator = AncestorIndex::new();
        Self {
            ancestors,
            ancestor_iterator,
        }
    }

    /// Find a copying path for the given sequence given the ancestor iterator as a partial
    /// tree sequence.
    fn find_copy_path(&self, candidate: &AncestralSequence) -> Vec<PartialSequenceEdge> {
        todo!()
    }

    /// Generate a tree sequence from the given ancestor array.
    /// This will modify the internal state to represent the tree sequence for all ancestors
    /// within the array.
    pub fn match_ancestors(&mut self) {
        todo!()
    }

    /// Insert a set of samples into an ancestral tree sequence. The samples will be matched
    /// against the existing sequence, but not against each other.
    pub fn match_samples(&mut self, samples: &[AncestralSequence]) {
        todo!()
    }

    /// Finalize the tree sequence.
    pub fn get_tree_sequence(&self) -> ! {
        todo!()
    }
}
