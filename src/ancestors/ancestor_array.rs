use crate::ancestors::{Ancestor, AncestralSequence};
use crate::variants::SequencePosition;
use crate::variants::VariantIndex;
use std::fs::File;
use std::io;
use std::io::Write;
use std::ops::{Deref, Index};
use std::path::Path;

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
#[derive(Debug, Clone)]
pub struct AncestorArray {
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

    /// Export the ancestor array in a TSV format that can be read by the test suite
    pub fn export_ancestors(&self, path: &Path) -> io::Result<()> {
        let mut node_file = path.to_path_buf();
        node_file.push("ancestors.tsv");
        let mut writer = File::create(node_file)?;

        writer.write_fmt(format_args!("start\tend\tage\tfocal_sites\tstate\n"))?;
        for ancestor in self.ancestors.iter() {
            ancestor.export(&mut writer)?;
        }

        Ok(())
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
