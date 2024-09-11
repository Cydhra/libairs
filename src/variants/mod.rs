use std::fmt::{Display, Formatter};
use std::ops::{Add, Sub};

mod builder;
mod data;
mod sequence;
mod site;

pub use builder::VariantDataBuilder;
pub use data::{SampleData, VariantData};
pub use sequence::VariantSequence;
pub use site::VariantSite;

/// A position in a DNA sequence. This newtype ensures that sequence positions and variant indices (indices into the
/// variant site vector) aren't mixed up.
#[derive(
    Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, serde::Serialize, serde::Deserialize,
)]
pub struct SequencePosition(u32);

impl SequencePosition {
    /// Create a new sequence position from a usize.
    pub fn from_usize(position: usize) -> Self {
        Self(position as u32)
    }

    #[inline]
    pub fn from_vec(positions: Vec<usize>) -> Vec<Self> {
        positions.into_iter().map(|p| Self::from_usize(p)).collect()
    }

    /// Get the underlying usize value of the sequence position.
    pub fn unwrap(&self) -> u32 {
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
#[derive(
    Copy, Clone, Debug, Ord, PartialOrd, Eq, PartialEq, Hash, serde::Serialize, serde::Deserialize,
)]
pub(crate) struct VariantIndex(pub(crate) u32);

impl VariantIndex {
    /// Create a new variant index from a raw index. Only for testing purposes, actual module code shouldn't work
    /// with raw values
    pub(crate) fn from_usize(index: usize) -> Self {
        Self(index as u32)
    }

    /// Get the next variant index after this one
    pub(crate) fn next(&self) -> Self {
        Self(self.0 + 1)
    }

    /// Calculate the distance in variants between this index and another. Does not return the distance in sequence bases.
    pub(crate) fn get_variant_distance(&self, other: VariantIndex) -> u32 {
        if self.0 > other.0 {
            self.0 - other.0
        } else {
            other.0 - self.0
        }
    }

    /// Get the underlying usize value of the variant index.
    pub fn unwrap(&self) -> u32 {
        self.0
    }
}

impl Display for VariantIndex {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}

impl Add for VariantIndex {
    type Output = VariantIndex;

    fn add(self, rhs: Self) -> Self::Output {
        VariantIndex(self.0 + rhs.0)
    }
}

impl Add<u32> for VariantIndex {
    type Output = VariantIndex;

    fn add(self, rhs: u32) -> Self::Output {
        VariantIndex(self.0 + rhs)
    }
}

impl Sub for VariantIndex {
    type Output = VariantIndex;

    fn sub(self, rhs: Self) -> Self::Output {
        VariantIndex(self.0 - rhs.0)
    }
}

impl Sub<VariantIndex> for u32 {
    type Output = VariantIndex;

    fn sub(self, rhs: VariantIndex) -> Self::Output {
        VariantIndex(self - rhs.0)
    }
}

impl Sub<u32> for VariantIndex {
    type Output = VariantIndex;

    fn sub(self, rhs: u32) -> Self::Output {
        VariantIndex(self.0 - rhs)
    }
}
