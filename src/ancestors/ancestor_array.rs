use crate::ancestors::AncestralSequence;
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

/// An index into the [`AncestorArray`]
///
/// [`AncestorArray`]: AncestorArray
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) struct VariantIndex(usize);

impl VariantIndex {
    #[cfg(test)]
    pub fn from_usize(index: usize) -> Self {
        Self(index)
    }

    /// Get the next variant index after this one
    pub(crate) fn next(&self) -> Self {
        Self(self.0 + 1)
    }
}

/// An index into the ancestor array which uniquely identifies an ancestor
#[derive(Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub(crate) struct Ancestor(pub(crate) usize);
