use crate::variants::VariantIndex;
use std::ops::{Deref, DerefMut, Index, IndexMut, Range};

pub type MutationState = u8;

/// A genomic sequence expressed only through variant calls at fixed mutation sites.
/// The data only makes sense with respect to a reference genome which defines the ancestral state,
/// and a [`super::VariantData`] instance that defines which index in the sequence corresponds to
/// which genomic site.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct VariantSequence(Vec<MutationState>);

impl VariantSequence {
    /// Create a new empty [`VariantSequence`]
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new [`VariantSequence`] with the given length, filled entirely with the ancestral
    /// state (i.e. with zeros).
    pub fn from_ancestral_state(len: usize) -> Self {
        Self(vec![0; len])
    }

    /// Create a new [`VariantSequence`] from a vector of mutation states
    pub fn from_vec(vec: Vec<MutationState>) -> Self {
        Self(vec)
    }
}

/// Decay into the underlying [`Vec`]
impl Deref for VariantSequence {
    type Target = Vec<MutationState>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for VariantSequence {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl Index<usize> for VariantSequence {
    type Output = MutationState;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for VariantSequence {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

/// Index a VariantSequence by a [`VariantIndex`] to get the state at the site.
impl Index<VariantIndex> for VariantSequence {
    type Output = MutationState;

    fn index(&self, index: VariantIndex) -> &Self::Output {
        &self.0[index.0]
    }
}

impl IndexMut<VariantIndex> for VariantSequence {
    fn index_mut(&mut self, index: VariantIndex) -> &mut Self::Output {
        &mut self.0[index.0]
    }
}

impl Index<Range<VariantIndex>> for VariantSequence {
    type Output = [MutationState];

    fn index(&self, index: Range<VariantIndex>) -> &Self::Output {
        &self.0[index.start.0..index.end.0]
    }
}

impl PartialEq<Vec<u8>> for VariantSequence {
    fn eq(&self, other: &Vec<u8>) -> bool {
        self.0.eq(other)
    }

    fn ne(&self, other: &Vec<u8>) -> bool {
        self.0.ne(other)
    }
}
