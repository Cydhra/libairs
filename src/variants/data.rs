use super::{SequencePosition, VariantIndex, VariantSequence, VariantSite};
use std::ops::{Index, Range};

/// Holds variant data and associated metadata.
/// This data is used as input for both ancestor generation and matching sample data against the
/// partial tree sequence.
pub struct VariantData {
    sites: Vec<VariantSite>,
    positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
    num_samples: usize,
}

impl VariantData {
    /// Construct a new VariantData instance
    pub(super) fn new(
        sites: Vec<VariantSite>,
        positions: Vec<SequencePosition>,
        sequence_length: SequencePosition,
        num_samples: usize,
    ) -> Self {
        Self {
            sites,
            positions,
            sequence_length,
            num_samples,
        }
    }

    /// Convert the variant data into sample data by transposing the variant sites into
    /// [`VariantSequence`]s.
    /// Data is copied, so the variant data object can still be used after this method is called.
    pub fn into_samples(&self) -> SampleData {
        // transpose sites into sample sequences
        let mut samples =
            vec![VariantSequence::from_ancestral_state(self.sites.len()); self.sites.len()];
        for (i, site) in self.sites.iter().enumerate() {
            for (j, state) in site.genotypes.iter().enumerate() {
                samples[j][VariantIndex(i)] = *state;
            }
        }

        SampleData {
            samples,
            positions: self.positions.clone(),
            sequence_length: self.sequence_length,
        }
    }

    /// Iterate through the [`VariantSite`]s in this instance
    pub fn iter<'s>(
        &'s self,
    ) -> impl Iterator<Item = &VariantSite> + ExactSizeIterator + DoubleEndedIterator + 's {
        self.sites.iter()
    }

    /// Conveniently iterate through the [`VariantSite`]s in this instance, yielding the
    /// [`VariantIndex`] of the site as well.
    pub(crate) fn iter_with_index<'s>(
        &'s self,
    ) -> impl Iterator<Item = (VariantIndex, &VariantSite)> + ExactSizeIterator + DoubleEndedIterator + 's
    {
        self.sites
            .iter()
            .enumerate()
            .map(|(i, s)| (VariantIndex(i), s))
    }

    /// Get the array that maps variant indices to sequence positions.
    /// TODO: maybe it's unnecessary to expose this here
    pub fn variant_positions(&self) -> &[SequencePosition] {
        &self.positions
    }

    /// Get the sequence length of the genome this variant data is about. It is not the length of the
    /// variant site vector, but the genome length.
    ///
    /// # Returns
    /// A [`SequencePosition`] containing the genome length
    pub fn get_sequence_length(&self) -> SequencePosition {
        self.sequence_length
    }

    /// Get the number of samples that make up each variant site
    pub fn get_num_samples(&self) -> usize {
        self.num_samples
    }

    /// Get the number of variant sites in the collection
    pub fn len(&self) -> usize {
        self.sites.len()
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

impl Index<Range<VariantIndex>> for VariantData {
    type Output = [VariantSite];

    fn index(&self, index: Range<VariantIndex>) -> &Self::Output {
        &self.sites[index.start.0..index.end.0]
    }
}

/// A collection of DNA sequence samples generated from variant data.
pub struct SampleData {
    samples: Vec<VariantSequence>,
    positions: Vec<SequencePosition>,
    sequence_length: SequencePosition,
}

impl SampleData {
    /// Get the sequence length of the genome this sample data is about. It is not the length of the
    /// variant site vector, but the genome length.
    ///
    /// # Returns
    /// A [`SequencePosition`] containing the genome length
    pub fn get_sequence_length(&self) -> SequencePosition {
        self.sequence_length
    }

    /// Iterate over the sample sequences in this instance
    pub fn iter<'s>(&'s self) -> impl Iterator + 's {
        self.samples.iter()
    }

    /// Get the number of samples in the collection
    pub fn len(&self) -> usize {
        self.samples.len()
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

/// Turn the sample data into an iterator over its sample sequences
impl IntoIterator for SampleData {
    type Item = VariantSequence;
    type IntoIter = <Vec<Self::Item> as IntoIterator>::IntoIter;

    fn into_iter(self) -> Self::IntoIter {
        self.samples.into_iter()
    }
}
