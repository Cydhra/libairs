use std::fs::File;
use std::io;
use std::io::Write;
use std::ops::{Deref, Index};
use std::path::Path;

use crate::ancestors::{Ancestor, AncestralSequence};
use crate::variants::{SampleData, VariantIndex};
use crate::variants::{SequencePosition, VariantData};

/// This is a helper struct for the Viterbi algorithm that manages the ancestral sequences.
#[derive(Clone, Debug)]
pub struct AncestorArray {
    ancestors: Vec<AncestralSequence>,

    /// The variant data that the ancestral sequences are based on
    variant_data: VariantData,
}
impl AncestorArray {
    pub(crate) fn new(ancestors: Vec<AncestralSequence>, variant_data: VariantData) -> Self {
        Self {
            ancestors,
            variant_data,
        }
    }

    /// Get the number of ancestors in the array
    pub fn len(&self) -> usize {
        self.ancestors.len()
    }

    /// Get the number of variants that the ancestors in the array are based on
    pub fn get_num_variants(&self) -> usize {
        self.variant_data.len()
    }

    pub(crate) fn iter(&self) -> impl Iterator<Item = (Ancestor, &AncestralSequence)> {
        self.ancestors
            .iter()
            .enumerate()
            .map(|(i, a)| (Ancestor(i), a))
    }

    /// Convert a variant index to a sequence position
    pub(crate) fn variant_index_to_sequence_pos(&self, index: VariantIndex) -> SequencePosition {
        self.variant_data.variant_index_to_sequence_pos(index)
    }

    /// Obtain the DNA sample data that makes up the variant data
    pub fn generate_sample_data(&self) -> SampleData {
        self.variant_data.into_samples()
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
