use std::fmt::{Display, Formatter};

mod data;
mod site;

pub use data::VariantData;
pub use site::VariantSite;

/// A position in a DNA sequence. This newtype ensures that sequence positions and variant indices (indices into the
/// variant site vector) aren't mixed up.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct SequencePosition(usize); // TODO do not make this public

impl SequencePosition {
    /// Create a new sequence position from a usize.
    pub fn from_usize(position: usize) -> Self {
        Self(position)
    }

    #[inline]
    pub fn from_vec(positions: Vec<usize>) -> Vec<Self> {
        positions.into_iter().map(|p| Self::from_usize(p)).collect()
    }

    /// Get the underlying usize value of the sequence position.
    pub fn unwrap(&self) -> usize {
        self.0
    }
}

impl Display for SequencePosition {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

/// An index into the [`VariantData`]. The new-type guarantees that variant indices aren't mixed up
/// with [`SequencePosition`].
#[derive(Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq)]
pub(crate) struct VariantIndex(pub(crate) usize);

impl VariantIndex {
    /// Create a new variant index from a raw index. Only for testing purposes, actual module code shouldn't work
    /// with raw values
    #[cfg(test)]
    pub(crate) fn from_usize(index: usize) -> Self {
        Self(index)
    }

    /// Get the next variant index after this one
    pub(crate) fn next(&self) -> Self {
        Self(self.0 + 1)
    }

    /// Calculate the distance in variants between this index and another. Does not return the distance in sequence bases.
    pub(crate) fn get_variant_distance(&self, other: VariantIndex) -> usize {
        if self.0 > other.0 {
            self.0 - other.0
        } else {
            other.0 - self.0
        }
    }

    /// Get the underlying usize value of the variant index.
    pub fn unwrap(&self) -> usize {
        self.0
    }
}

impl Display for VariantIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}
