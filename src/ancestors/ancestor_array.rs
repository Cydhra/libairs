use crate::ancestors::{Ancestor, AncestralSequence};
use std::ops::{Deref, Index};

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
pub struct AncestorArray {
    // TODO figure out if transposing the ancestors improves cache locality
    ancestors: Vec<AncestralSequence>,
}
impl AncestorArray {
    pub(crate) fn from(ancestors: Vec<AncestralSequence>) -> Self {
        Self { ancestors }
    }

    /// Get the number of ancestors in the array
    pub fn len(&self) -> usize {
        self.ancestors.len()
    }
}

impl Index<Ancestor> for AncestorArray {
    type Output = AncestralSequence;

    fn index(&self, index: Ancestor) -> &Self::Output {
        &self.ancestors[index.0]
    }
}

impl Deref for AncestorArray {
    type Target = Vec<AncestralSequence>;

    fn deref(&self) -> &Self::Target {
        &self.ancestors
    }
}
