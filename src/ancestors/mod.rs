use std::fmt::{Debug, Display, Formatter};
use std::io;
use std::io::Write;
use std::ops::{Deref, Index, IndexMut};

pub use ancestor_array::AncestorArray;
pub use generator::AncestorGenerator;

use crate::variants::{VariantIndex, VariantSequence};

mod ancestor_array;
mod generator;

const ANCESTRAL_STATE: u8 = 0;
const DERIVED_STATE: u8 = 1;

/// A DNA sequence expressed through a bit vector where each bit defines whether the DNA sequence at
/// the given site has the ancestral state, or the derived state. This only works if we do not accept
/// more than two variants (ancestral and one derived state) per site.
#[derive(Clone, serde::Serialize, serde::Deserialize)]
pub struct AncestralSequence {
    state: VariantSequence,
    focal_sites: Vec<VariantIndex>,
    /// start of valid data in the state vector, inclusive
    start: VariantIndex,
    /// end of valid data in the state vector, exclusive
    end: VariantIndex,
    age: f64,
}

impl AncestralSequence {
    fn from_ancestral_state(len: usize, age: f64) -> Self {
        AncestralSequence {
            state: VariantSequence::from_ancestral_state(len),
            focal_sites: Vec::new(),
            start: VariantIndex(0),
            end: VariantIndex(len),
            age,
        }
    }

    /// Create a new ancestral sequence with state copied from another [`VariantSequence`].
    fn copy_from(
        variant_sequence: &VariantSequence,
        start: VariantIndex,
        end: VariantIndex,
        focal_sites: Vec<VariantIndex>,
        age: f64,
    ) -> Self {
        AncestralSequence {
            state: VariantSequence::from_vec(variant_sequence[start..end].to_vec()),
            focal_sites,
            start,
            end,
            age,
        }
    }

    /// Create an ancestral sequence that contains sample data from an actual input DNA sequence.
    /// The age will be 0, as this sequence is not inferred from any other sequence.
    /// The start and end indices are set to the full length of the sequence.
    pub(crate) fn new_sample(sample_sequence: VariantSequence) -> Self {
        let sequence_len = sample_sequence.len();
        AncestralSequence {
            state: sample_sequence,
            focal_sites: Vec::new(),
            start: VariantIndex(0),
            end: VariantIndex(sequence_len),
            age: 0.0,
        }
    }

    /// Get the haplotype sequence for the ancestral sequence. This sequence starts at the first
    /// known site and ends at the last known site. Where this sequence is located in the genome
    /// is defined by [`start`] and [`end`] respectively.
    /// This means that the indices of the haplotype sequence do not correspond to the indices of
    /// the genome.
    ///
    /// [`start`]: Self::start
    /// [`end`]: Self::end
    pub fn haplotype(&self) -> &[u8] {
        &self.state
    }

    /// Get an enumerated iterator over the haplotype sequence. This sequence starts at the first
    /// known site and ends at the last known site. Each returned element also contains the index
    /// of the site in the genome (regarding the genome's variant site vector, the actual genome
    /// position will differ from this).
    pub(crate) fn site_iter(
        &self,
    ) -> impl Iterator<Item = (VariantIndex, &'_ u8)> + DoubleEndedIterator + '_ {
        self.state
            .iter()
            .zip(self.start.0..self.end.0)
            .map(|(bit, idx)| (VariantIndex(idx), bit))
    }

    /// Get the position of the first known site in the genome (inclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub(crate) fn start(&self) -> VariantIndex {
        self.start
    }

    /// Get the position of the last known site in the genome (exclusive), regarding the genome's variant site
    /// vector (so it doesn't necessarily correspond to the actual position in the genome).
    pub(crate) fn end(&self) -> VariantIndex {
        self.end
    }

    /// Get the inferred relative age of the ancestral sequence. This is derived from the inferred
    /// relative age of the focal sites that were used to infer this ancestor.
    pub fn relative_age(&self) -> f64 {
        self.age
    }

    /// Get the length of the ancestral sequence. Only known sites are considered, so the length
    /// might be shorter than the length of the underlying DNA sequence.
    pub fn len(&self) -> usize {
        self.end.get_variant_distance(self.start)
    }

    /// Dump the ancestral sequence into a text file for the testing suite.
    pub fn export(&self, writer: &mut dyn Write) -> io::Result<()> {
        writer.write_fmt(format_args!(
            "{start}\t{end}\t{age}\t{focal_sites:?}\t",
            start = self.start.0,
            end = self.end.0,
            age = self.age,
            focal_sites = self.focal_sites.iter().map(|s| s.0).collect::<Vec<_>>(),
        ))?;

        for b in self.state.deref().iter() {
            writer.write_fmt(format_args!("{}", b))?;
        }

        writer.write_fmt(format_args!("\n"))?;
        Ok(())
    }
}

impl Debug for AncestralSequence {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.write_str("AncestralSequence { ")?;
        f.write_fmt(format_args!("focal_sites={:?},\t", self.focal_sites))?;

        f.write_str("genotype=[ ")?;
        let mut iter = self.state.iter();
        let mut idx = 0;
        while let Some(b) = iter.next() {
            if idx < self.start.0 || idx >= self.end.0 {
                f.write_str("-")?;
            } else {
                f.write_fmt(format_args!("{}", b))?;
            }

            if iter.len() > 0 {
                f.write_str(", ")?;
            }

            idx += 1;
        }
        f.write_str(" ] }")?;
        Ok(())
    }
}

impl Index<VariantIndex> for AncestralSequence {
    type Output = u8;

    fn index(&self, index: VariantIndex) -> &Self::Output {
        &self.state[index - self.start]
    }
}

impl IndexMut<VariantIndex> for AncestralSequence {
    fn index_mut(&mut self, index: VariantIndex) -> &mut Self::Output {
        &mut self.state[index - self.start]
    }
}

/// An index into the ancestor array which uniquely identifies an ancestor
#[derive(
    Copy, Clone, Debug, Eq, PartialEq, Ord, PartialOrd, Hash, serde::Serialize, serde::Deserialize,
)]
pub(crate) struct Ancestor(pub(crate) usize);

impl Display for Ancestor {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        Display::fmt(&self.0, f)
    }
}
