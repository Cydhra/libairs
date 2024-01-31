use crate::ancestors::{Ancestor, AncestralSequence};
use std::ops::Index;

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
pub(crate) struct AncestorArray {
    // TODO figure out if transposing the ancestors improves cache locality
    ancestors: Vec<AncestralSequence>,
    num_variants: usize,
}
impl AncestorArray {
    pub(crate) fn from(ancestors: Vec<AncestralSequence>, num_variants: usize) -> Self {
        Self {
            ancestors,
            num_variants,
        }
    }
}

impl Index<Ancestor> for AncestorArray {
    type Output = AncestralSequence;

    fn index(&self, index: Ancestor) -> &Self::Output {
        &self.ancestors[index.0]
    }
}
