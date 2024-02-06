use crate::ancestors::{Ancestor, AncestralSequence, VariantIndex};
use crate::dna::SequencePosition;
use std::ops::{Deref, Index};

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
#[derive(Debug, Clone)]
pub struct AncestorArray {
    // TODO figure out if transposing the ancestors improves cache locality
    ancestors: Vec<AncestralSequence>,

    /// Maps variant indices to sequence positions
    variant_positions: Vec<SequencePosition>,

    /// The length of the sequence
    sequence_length: SequencePosition,
}
impl AncestorArray {
    pub(crate) fn new(
        ancestors: Vec<AncestralSequence>,
        variant_positions: Vec<SequencePosition>,
        sequence_length: SequencePosition,
    ) -> Self {
        Self {
            ancestors,
            variant_positions,
            sequence_length,
        }
    }

    /// Get the number of ancestors in the array
    pub fn len(&self) -> usize {
        self.ancestors.len()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (Ancestor, &AncestralSequence)> {
        self.ancestors
            .iter()
            .enumerate()
            .map(|(i, a)| (Ancestor(i), a))
    }

    /// Convert a variant index to a sequence position
    pub(crate) fn variant_index_to_sequence_pos(&self, index: VariantIndex) -> SequencePosition {
        if index.0 == 0 {
            SequencePosition::from_usize(0)
        } else if index.0 == self.variant_positions.len() {
            self.sequence_length
        } else {
            self.variant_positions[index.0]
        }
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
