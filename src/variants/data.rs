use super::{SequencePosition, VariantIndex, VariantSite};
use std::ops::Index;

/// Holds variant data and associated metadata.
/// This data is used as input for both ancestor generation and matching sample data against the
/// partial tree sequence.
pub struct VariantData {
    sites: Vec<VariantSite>,
    positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
}

impl VariantData {
    fn new() -> Self {
        todo!("implement")
    }

    /// Get the sequence length of the genome this variant data is about. It is not the length of the
    /// variant site vector, but the genome length.
    ///
    /// # Returns
    /// A [`SequencePosition`] containing the genome length
    fn get_sequence_length(&self) -> SequencePosition {
        self.sequence_length
    }

    /// Convert a variant index to a sequence position
    pub(crate) fn variant_index_to_sequence_pos(&self, index: VariantIndex) -> SequencePosition {
        if index.0 == 0 {
            SequencePosition::from_usize(0)
        } else if index.0 == self.positions.len() {
            self.sequence_length
        } else {
            self.positions[index.0]
        }
    }
}

/// Index variant data by [`VariantIndex`].
impl Index<VariantIndex> for VariantData {
    type Output = VariantSite;

    fn index(&self, index: VariantIndex) -> &Self::Output {
        &self.sites[index.0]
    }
}
