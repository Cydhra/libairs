use std::ops::Index;
use crate::ancestors::AncestralSequence;

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
pub(crate) struct AncestorArray {
    // TODO figure out if transposing the ancestors improves cache locality
    ancestors: Vec<AncestralSequence>,
}

impl AncestorArray {
    pub(crate) fn from(ancestors: Vec<AncestralSequence>) -> Self {
        Self {
            ancestors
        }
    }
}

impl Index<Ancestor> for AncestorArray {
    type Output = AncestralSequence;

    fn index(&self, index: Ancestor) -> &Self::Output {
        &self.ancestors[index.0]
    }
}

/// An index into the [`AncestorArray`]
///
/// [`AncestorArray`]: AncestorArray
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) struct VariantIndex(usize);

/// An index into the ancestor array which uniquely identifies an ancestor
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub(crate) struct Ancestor(usize);