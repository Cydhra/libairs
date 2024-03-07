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
